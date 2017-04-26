import sys
import re
import os
import tempfile
import subprocess
import shutil
import multiprocessing
import logging
import csv

import Bio
import Bio.SeqIO
import Bio.SeqIO.QualityIO
import Bio.Alphabet.IUPAC
import Bio.Blast
import Bio.Blast.NCBIXML

import Bio.File
import Bio.SeqIO.FastaIO

import paftol.tools


logger = logging.getLogger(__name__)

keepTmp = False


def isSane(filename):
    """Check whether a file name is sane, in the sense that it does not contain any "funny" characters"""
    if filename == '':
        return False
    funnyCharRe = re.compile('\t/ ;,$#')
    m = funnyCharRe.search(filename)
    if m is not None:
        return False
    if filename[0] == '-':
        return False
    return True


def cmpExonerateResultByQueryAlignmentStart(e1, e2):
    """Comparator function for sorting C{ExonerateResult}s by query alignment start.

@param e1: first exonerate result
@type e1: C{ExonerateResult}
@param e2: second exonerate result
@type e2: C{ExonerateResult}
@return: one of -1, 0 or 1
@rtype: C{int}
"""
    if e1.queryAlignmentStart < e2.queryAlignmentStart:
        return -1
    elif e1.queryAlignmentStart > e2.queryAlignmentStart:
        return 1
    return 0


# FIXME: use abc for this class?
class HybseqAnalyser(object):
    """Base class for Hybseq analysers.
    
Instances of this class take a FASTA file of target PAFTOL gene sequences
and FASTQ files (one or two, for single / paired end, respectively),
and provide methods for running analyses to reconstruct sequences of
the target genes.
"""

    def __init__(self, targetsSourcePath, forwardFastq, reverseFastq=None, workdirTgz=None, workDirname='paftoolstmp'):
        self.targetsSourcePath = targetsSourcePath
        self.forwardFastq = forwardFastq
        self.reverseFastq = reverseFastq
        self.workdirTgz = workdirTgz
        self.workDirname = workDirname
        self.tmpDirname = None
        # parameters for ensuring file names don't clash, e.g. because paftolGene / organism name is same as targets basename etc.
        self.targetsFname = 'targets.fasta'
        self.geneFnamePattern = 'gene-%s.fasta'
        self.allowInvalidBases = False
    
    def __str__(self):
        return 'HybseqAnalyser(targetsSourcePath=%s, forwardFastq=%s, reverseFastq=%s)' % (repr(self.targetsSourcePath), repr(self.forwardFastq), repr(self.reverseFastq))
        
    def checkTargets(self):
        # FIXME: merge with __init__()? parsing is redundant with HybpiperAnalyser.initPaftolTargetDicts too
        for targetSr in Bio.SeqIO.parse(self.targetsSourcePath, 'fasta', alphabet = Bio.Alphabet.IUPAC.ambiguous_dna):
            if not self.allowInvalidBases:
                setDiff = set(str(targetSr.seq).lower()) - set('acgt')
                if len(setDiff) != 0:
                    raise StandardError('target %s: illegal base(s) %s' % (targetSr.id, ', '.join(setDiff)))

    def isPaired(self):
        return self.reverseFastq is not None
    
    def analyse(self):
        raise StandardError('not implemented in this "abstract" base class')
    
    def setupTmpdir(self):
        if self.tmpDirname is not None:
            raise StandardError('illegal state: already have generated working directory %s' % self.tmpDirname)
        self.tmpDirname = tempfile.mkdtemp(prefix=self.workDirname)
        os.mkdir(self.makeWorkDirname())

    def cleanupTmpdir(self):
        if self.tmpDirname is not None:
            if keepTmp:
                logger.warning('not removing temporary directory %s', self.tmpDirname)
            else:
                shutil.rmtree(self.tmpDirname)
            self.tmpDirname = None

    def makeWorkDirname(self):
        if self.tmpDirname is None:
            raise StandardError('illegal state: no temporary directory and hence no working directory')
        # logger.debug('tmpDirname = %s, workDirname = %s', self.tmpDirname, self.workDirname)
        return os.path.join(self.tmpDirname, self.workDirname)
    
    def makeTargetsFname(self, absolutePath=False):
        if absolutePath:
            return os.path.join(self.makeWorkDirname(), self.targetsFname)
        else:
            return self.targetsFname

    def makeGeneFname(self, geneName, absolutePath=False):
        geneFname = self.geneFnamePattern % geneName
        if absolutePath:
            return os.path.join(self.makeWorkDirname(), geneFname)
        else:
            return geneFname
        
    def makeTgz(self):
        if self.workdirTgz is not None:
            if self.tmpDirname is None:
                raise StandardError('illegal state: no temporary directory generated')
            tmpTgz = os.path.join(self.tmpDirname, '%s.tgz' % self.workDirname)
            tgzArgv = ['tar', '-zcf', tmpTgz, self.workDirname]
            tgzProcess = subprocess.Popen(tgzArgv, cwd = self.tmpDirname)
            tgzReturncode = tgzProcess.wait()
            if tgzReturncode != 0:
                raise StandardError('process "%s" returned %d' % (' '.join(tgzArgv), tgzReturncode))
            # FIXME: clumsy to first create tgz in temp dir and then
            # moving it to final destination, compute absolute path to
            # final destination and use that directly?
            shutil.move(os.path.join(self.tmpDirname, tmpTgz), self.workdirTgz)


class SamAlignment(object):
    """Class to represent a SAM record.
This class follows the naming and definitions of the SAMv1 spec. It is incomplete
to provide fields required for Hyb-Seq analysis only.

@ivar qname: SAM query name (C{QNAME}), read or read pair ID
@type qname: C{str}
@ivar rname: SAM reference name (C{RNAME})
@type rname: C{str}
@ivar mapq: SAM mapping quality (C{MAPQ})
@type mapq: C{int}
@ivar cigar: SAM CIGAR string (unexpanded) (C{CIGAR})
@type mapq: C{str}
@ivar seq: SAM query (read) sequence (C{SEQ})
@type seq: C{str}
"""
    
    cigarElementRe = re.compile('([0-9]+)([MIDNSHP=X])')
    
    def __init__(self, samLine):
        if samLine[-1] == '\n':
            samLine = samLine[:-1]
        w = samLine.split('\t')
        self.qname = w[0]
        self.rname = w[2]
        self.mapq = int(w[4])
        self.cigar = w[5]
        self.seq = w[9]
    
    def expandedCigar(self):
        if self.cigar is None:
            return None
        e = ''
        c = self.cigar
        while c != '':
            m = self.cigarElementRe.match(c)
            if m is None:
                raise StandardError('malformed CIGAR "%s" (stuck at "%s")' % (self.cigar, c))
            e = e + (m.group(2) * int(m.group(1)))
            c = c[len(m.group()):]
        return e
    
    def numCigarMatches(self):
        e = self.expandedCigar()
        if e is None:
            return None
        if e.count('=') > 0:
            logger.warning('found sequence match ("=") characters, unimplemented')
        if e.count('X') > 0:
            logger.warning('found sequence mismatch ("X") characters, unimplemented')
        return e.count('M')


class PaftolTarget(object):
    """Represent a PAFTOL target, specific to an organism (i.e. species, specimen etc.).
    
The main content of instances of this class is a C{SeqRecord}
containing the sequence of the gene in the organism, thus
facilitating handling of multiple genes and multiple organisms.

@ivar organism: the organism
@type organism: C{Organism}
@ivar paftolGene: the PAFTOL gene 
@type paftolGene: C{PaftolGene}
@ivar seqRecord: the sequence of this gene in this organism
@type seqRecord: C{Bio.SeqRecord.SeqRecord}
"""
    def __init__(self, organism, paftolGene, seqRecord):
        self.organism = organism
        self.paftolGene = paftolGene
        self.seqRecord = seqRecord
        self.samAlignmentList = []
        if paftolGene.name in organism.paftolTargetDict or organism.name in paftolGene.paftolTargetDict:
            raise StandardError('duplicate organism/gene: organism = %s, gene = %s, seqId = %s' % (organism.name, paftolGene.name, seqRecord.id))
        organism.paftolTargetDict[paftolGene.name] = self
        paftolGene.paftolTargetDict[organism.name] = self
        
    def addSamAlignment(self, samAlignment):
        self.samAlignmentList.append(samAlignment)
            
    def mapqSum(self):
        if len(self.samAlignmentList) == 0:
            return None
        return sum([a.mapq for a in self.samAlignmentList])
    
    def qnameSet(self):
        # FIXME: may have to trim away "/1", "/2"?
        return set([a.qname for a in self.samAlignmentList])
    
    @staticmethod
    def makeCsvDictWriter(csvfile):
        csvFieldnames = ['organism', 'gene', 'seqLength', 'numSamAlignments']
        csvDictWriter = csv.DictWriter(csvfile, csvFieldnames)
        csvDictWriter.writeheader()
        return csvDictWriter
    
    def writeCsvRow(self, csvDictWriter):
        d = {}
        d['organism'] = self.organism.name
        d['gene'] = self.paftolGene.name
        d['seqLength'] = len(self.seqRecord)
        d['numSamAlignments'] = len(self.samAlignmentList)
        logger.debug('writing CSV row: %s, %s', self.organism.name, self.paftolGene.name)
        csvDictWriter.writerow(d)


class Organism(object):
    """Represent an organism (in the GenBank / NCBI sense of the term).
    
@ivar name: this organism's name
@type name: C{str}
@ivar paftolTargetDict: dictionary of genes in this organism
@type paftolTargetDict: C{dict} of C{PaftolTarget} instances with PAFTOL gene names as keys
"""
    
    def __init__(self, name):
        self.name = name
        self.paftolTargetDict = {}


class PaftolGene(object):
    """Represent a PAFTOL gene.
    
This class does not represent genes in terms of intron / exon models
and other features. Its main purpose is to contain a collection of
PAFTOL targets, i.e. sequences found for this gene in various
organisms.

@ivar name: the name of this PAFTOL gene
@type name: C{str}
@ivar paftolTargetDict: dictionary of organisms with this PAFTOL gene
@type paftolTargetDict: C{dict} of C{PaftolTarget} instances with organism names as keys

    """
    
    def __init__(self, name):
        self.name = name
        self.paftolTargetDict = {}

    def qnameSet(self):
        s = set()
        for paftolTarget in self.paftolTargetDict.values():
            s = s | paftolTarget.qnameSet()
        return s
    

class PaftolTargetSet(object):

    paftolTargetRe = re.compile('([^-]+)-([^-]+)')
    
    def __init__(self):
        self.paftolGeneDict = {}
        self.organismDict = {}
        
    def makeFastaId(organismName, geneName):
        return '%s-%s' % (organismName, geneName)
        
    def extractOrganismAndGeneNames(self, s):
        m = self.paftolTargetRe.match(s)
        if m is not None:
            organismName = m.group(1)
            geneName = m.group(2)
        else:
            organismName = 'unknown'
            geneName = s
        return organismName, geneName
    
    def readFasta(self, fastaHandle):
        self.paftolGeneDict = {}
        self.organismDict = {}
        for sr in Bio.SeqIO.parse(fastaHandle, 'fasta'):
            organismName, geneName = self.extractOrganismAndGeneNames(sr.id)
            if not isSane(organismName):
                raise StandardError('bad organism name: %s' % organismName)
            if not isSane(geneName):
                raise StandardError('bad gene name: %s' % geneName)
            if organismName not in self.organismDict:
                self.organismDict[organismName] = Organism(organismName)
            if geneName not in self.paftolGeneDict:
                self.paftolGeneDict[geneName] = PaftolGene(geneName)
            paftolTarget = PaftolTarget(self.organismDict[organismName], self.paftolGeneDict[geneName], sr)
            
    def getSeqRecordList(self):
        srList = []
        for organism in self.organismDict.values():
            for paftolTarget in organism.paftolTargetDict.values():
                srList.append(paftolTarget.seqRecord)
        return srList
            
    def writeFasta(self, fastaHandle):
        srList = self.getSeqRecordList()
        sys.stderr.write('writeFasta: writing %d sequences\n' % len(srList))
        Bio.SeqIO.write(srList, fastaHandle, 'fasta')
        
    def addSamAlignment(self, samAlignment):
        organismName, geneName = self.extractOrganismAndGeneNames(samAlignment.rname)
        if organismName not in self.organismDict:
            raise StandardError('unknown organism: %s' % organismName)
        if geneName not in self.paftolGeneDict:
            raise StandardError('unknown gene: %s' % geneName)
        if geneName not in self.organismDict[organismName].paftolTargetDict:
            raise StandardError('no entry for gene %s in organism %s' % (geneName, organismName))
        paftolTarget = self.organismDict[organismName].paftolTargetDict[geneName]
        paftolTarget.addSamAlignment(samAlignment)
    
    
class ReferenceGene(object):
    
    def __init__(self, geneId, referenceGenome, seqRecord, geneFeature, mrnaFeature=None, cdsFeature=None):
        self.geneId = geneId
        self.referenceGenome = referenceGenome
        self.seqRecord = seqRecord
        self.geneFeature = geneFeature
        self.mrnaFeature = mrnaFeature
        self.cdsFeature = cdsFeature
        
    def getSequenceId(self):
        return self.seqRecord.id.split('.')[0]
        
    def containsHsp(self, hspAccession, hsp):
        if self.getSequenceId() != hspAccession:
            return False
        return self.geneFeature.location.start <= hsp.sbjct_start and self.geneFeature.location.end >= hsp.sbjct_end
                
    
class ReferenceGenome(object):
    """Represent a reference genome, provided via FASTA and GenBank files (possibly both).

@ivar name: the name of this reference genome
@type name: C{str}
@ivar fastaFname: name of FASTA file containing the sequences of this genome
@type fastaFname: C{str}
@ivar genbankFname: name of GenBank file containing the sequences of this genome
@type genbankFname: C{str}
"""
    
    def __init__(self, name, fastaFname, genbankFname):
        self.name = name
        self.fastaFname = fastaFname
        self.genbankFname = genbankFname
        self.geneList = None
        
    def scanGenesAth(self):
        if self.genbankFname is None:
            raise StandardError('no GenBank file name, cannot scan genes (ath method)')
        mrnaFeatureDict = {}
        cdsFeatureDict = {}
        self.geneList = []
        geneDict = {}
        with open(self.genbankFname, 'r') as f:
            for seqRecord in Bio.SeqIO.parse(f, 'genbank'):
                for seqFeature in seqRecord.features:
                    if seqFeature.type == 'gene':
                        # CHECKME: just presuming that locus_tag qualifier will always be present and have exactly one value
                        geneId = seqFeature.qualifiers['locus_tag'][0]
                        if geneId in geneDict:
                            raise StandardError('duplicate gene id: %s' % geneId)
                        gene = ReferenceGene(geneId, self, seqRecord, seqFeature)
                        self.geneList.append(gene)
                        geneDict[geneId] = gene
        # somewhat clumsy to re-scan GenBank file for additional features...
        with open(self.genbankFname, 'r') as f:
            for seqRecord in Bio.SeqIO.parse(f, 'genbank'):
                for seqFeature in seqRecord.features:
                    if seqFeature.type == 'mRNA':
                        geneId = seqFeature.qualifiers['locus_tag'][0]
                        if geneId in geneDict:
                            gene = geneDict[geneId]
                            if gene.mrnaFeature is not None:
                                sys.stderr.write('gene %s: duplicate mRNA feature, ignoring\n' % geneId)
                            else:
                                gene.mrnaFeature = seqFeature
                    elif seqFeature.type == 'CDS':
                        geneId = seqFeature.qualifiers['locus_tag'][0]
                        if geneId in geneDict:
                            gene = geneDict[geneId]
                            if gene.cdsFeature is not None:
                                sys.stderr.write('gene %s: duplicate CDS feature, ignoring\n' % geneId)
                            else:
                                gene.cdsFeature = seqFeature

    def scanGenes(self, scanMethod):
        """Populate C{self.geneList} by scanning an appropriate file.

Currently, the only method supported is C{ath}, which is designed to
work with the Arabidopsis thaliana genome (specifically the TAIR10
release). Currently, specifying any other method will raise an
exception. In the future, more genomes with different annotation
conventions may be added.
        
@param scanMethod: the method to use for scanning genes
@type scanMethod: C{str}
        """
        if scanMethod == 'ath':
            self.scanGenesAth()
        else:
            raise StandardError('unknown gene scan method: %s' % scanMethod)
        
    def findGenesByHsp(self, hspAccession, hsp):
        """Find genes that contain a given HSP.
"""
        geneList = []
        for gene in self.geneList:
            if gene.containsHsp(hspAccession, hsp):
                geneList.append(gene)
        return geneList

    def blastTargetSet(self, paftolTargetSet):
        blastnArgv = ['blastn', '-db', self.fastaFname, '-outfmt', '5']
        # blastnArgv = ['tee', 'tee.txt']
        logger.debug('%s', ' '.join(blastnArgv))
        sys.stderr.flush()
        sys.stdout.flush()
        blastnProcess = subprocess.Popen(blastnArgv, stdin=subprocess.PIPE, stdout = subprocess.PIPE)
        subprocess.call(['lsof', '-p', '%d' % os.getpid()])
        # blastnProcess.stdin.flush()
        pid = os.fork()
        if pid == 0:
            # reload(Bio.SeqIO)
            blastnProcess.stdout.close()
            # paftolTargetSet.writeFasta(sys.stderr)
            # srList = paftolTargetSet.getSeqRecordList()
            # sys.stderr.write('target set has %d seqRecords\n' % len(srList))
            sr = Bio.SeqRecord.SeqRecord(Bio.Seq.Seq('A'), id = 'srDummy', description = '')
            # sq = Bio.Seq.Seq('A')
            # s = str(sr.seq)
            # Bio.SeqIO.write([sr], sys.stderr, 'fasta')
            # sys.stderr.write(sr.format('fasta'))
            # for i in xrange(333):
            #     sys.stderr.write('>dummy\n')
            #     for j in xrange(10):
            #         sys.stderr.write('%s\n' % ('A' * 60))
            # with Bio.File.as_handle(sys.stderr) as h:
                # sys.stderr.write('biopython version: %s\n' % Bio.__version__)
                # sys.stderr.write('handle of sys.stderr: %s\n' % str(h))
                # w = Bio.SeqIO.FastaIO.FastaWriter(h)
                # sys.stderr.write('writer: %s\n' % str(w))
                # w.write_file([sr])
                # w.write_header()
                # w.write_records([sr])
                # w.handle.write('>someseq\n')
                # s = str(sr.seq)
                # w.handle.write(str(sr.seq) + '\n')
                # w.handle.write(sr.format('fasta'))
            # x = sr.format('fasta')
            paftolTargetSet.writeFasta(blastnProcess.stdin)
            # for sr in paftolTargetSet.getSeqRecordList():
            #     blastnProcess.stdin.write(sr.format('fasta'))
            blastnProcess.stdin.close()
            os._exit(0)
        blastnProcess.stdin.close()
        targetIdToGeneDict = {}
        for blastRecord in Bio.Blast.NCBIXML.parse(blastnProcess.stdout):
            targetId = blastRecord.query
            if targetId in targetIdToGeneDict:
                raise StandardError('duplicate BLAST record for target %s' % targetId)
            geneList = []
            for blastAlignment in blastRecord.alignments:
                for hsp in blastAlignment.hsps:
                    for gene in self.findGenesByHsp(blastAlignment.accession, hsp):
                        if gene not in geneList:
                            geneList.append(gene)
            targetIdToGeneDict[targetId] = geneList
        blastnProcess.stdout.close()
        wPid, wExit = os.waitpid(pid, 0)
        if pid != wPid:
            raise StandardError('wait returned pid %s (expected %d)' % (wPid, pid))
        if wExit != 0:
            raise StandardError('wait on forked process returned %d' % wExit)
        blastnReturncode = blastnProcess.wait()
        if blastnReturncode != 0:
            raise StandardError('blastn process exited with %d' % blastnReturncode)
        return targetIdToGeneDict


class HybpiperAnalyser(HybseqAnalyser):
    """L{HybseqAnalyser} subclass that implements an analysis process
close to the HybPiper pipeline.

Some parameters to BWA and SPAdes can be controlled via instance
variables as documented below. Defaults of these parameters correspond
to the defaults provided by BWA and SPAdes, respectively (at the time
of developing this).

@ivar bwaNumThreads: BWA number of threads (C{-t} option)
@type bwaNumThreads: C{int}
@ivar bwaMinSeedLength: BWA minimum seed length (C{-k} option)
@type bwaMinSeedLength: C{int}
@ivar bwaScoreThreshold: BWA score threshold for recording reads as mapped (C{-T} option)
@type bwaScoreThreshold: C{int}
@ivar bwaReseedTrigger: BWA re-seed trigger (C{-r} option)
@type bwaScoreThreshold: C{float}
@ivar spadesCovCutoff: SPAdes coverage cutoff (C{--cov-cutoff} option)
@type spadesCovCutoff: C{int}
@ivar spadesKvalList: SPAdes oligomer length value list (C{-k} option)
@type spadesKvalList: C{list} of C{int}, or C{None}
"""
    
    def __init__(self, targetsSourcePath, forwardFastq, reverseFastq=None, workdirTgz=None, workDirname='pafpipertmp'):
        super(HybpiperAnalyser, self).__init__(targetsSourcePath, forwardFastq, reverseFastq, workdirTgz, workDirname)
        self.bwaNumThreads = 1
        self.bwaMinSeedLength = 19
        self.bwaScoreThreshold = 30
        self.bwaReseedTrigger = 1.5
        self.spadesCovCutoff = 8
        self.spadesKvalList = None
        self.statsCsvFilename = None
        self.exoneratePercentIdentityThreshold = 65.0
        self.initPaftolTargetDicts()
        
    def initPaftolTargetDicts(self):
        if self.targetsSourcePath is None:
            raise StandardError('illegal state: cannot init organism and gene dicts with targetsSourcePath = None')
        self.paftolTargetSet = PaftolTargetSet()
        self.paftolTargetSet.readFasta(self.targetsSourcePath)
        logger.info('%s organisms, %s genes' % (len(self.paftolTargetSet.organismDict), len(self.paftolTargetSet.paftolGeneDict)))
        self.representativePaftolTargetDict = None

    def setup(self):
        logger.debug('setting up')
        if self.targetsSourcePath is None:
            raise StandardError('illegal state: cannot set up with targetsSourcePath = None')
        self.setupTmpdir()
        shutil.copy(self.targetsSourcePath, self.makeTargetsFname(True))
            
    def cleanup(self):
        self.cleanupTmpdir()
    
    def bwaIndexReference(self, referenceFname):
        bwaIndexArgv = ['bwa', 'index', referenceFname]
        logger.debug('%s', ' '.join(bwaIndexArgv))
        subprocess.check_call(bwaIndexArgv)
        
    def makeBwaMemArgv(self, referenceFname, fastqArgs):
        return ['bwa', 'mem', '-M', '-k', '%d' % self.bwaMinSeedLength, '-T', '%d' % self.bwaScoreThreshold, '-r', '%f' % self.bwaReseedTrigger, '-t', '%d' % self.bwaNumThreads, referenceFname] + fastqArgs
        
    def mapReadsBwa(self):
        """Map reads to gene sequences (from multiple organisms possibly).
"""
        logger.debug('mapping reads to gene sequences')
        self.bwaIndexReference(self.makeTargetsFname(True))
        fastqArgs = [os.path.join(os.getcwd(), self.forwardFastq)]
        if self.reverseFastq is not None:
            fastqArgs.append(os.path.join(os.getcwd(), self.reverseFastq))
        # bwa parameters for tweaking considerations: -k, -r, -T
        bwaArgv = self.makeBwaMemArgv(self.makeTargetsFname(), fastqArgs)
        samtoolsArgv = ['samtools', 'view', '-h', '-S', '-F', '4', '-']
        logger.debug('%s', ' '.join(bwaArgv))
        bwaProcess = subprocess.Popen(bwaArgv, stdout=subprocess.PIPE, cwd = self.makeWorkDirname())
        logger.debug('%s', ' '.join(samtoolsArgv))
        samtoolsProcess = subprocess.Popen(samtoolsArgv, stdin=bwaProcess.stdout.fileno(), stdout=subprocess.PIPE, cwd = self.makeWorkDirname())
        samLine = samtoolsProcess.stdout.readline()
        while samLine != '':
            # logger.debug(samLine)
            if samLine[0] != '@':
                samAlignment = SamAlignment(samLine)
                self.paftolTargetSet.addSamAlignment(samAlignment)
            samLine = samtoolsProcess.stdout.readline()
        bwaProcess.stdout.close()
        samtoolsProcess.stdout.close()
        bwaReturncode = bwaProcess.wait()
        samtoolsReturncode = samtoolsProcess.wait()
        if bwaReturncode != 0:
            raise StandardError('process "%s" returned %d' % (' '.join(bwaArgv), bwaReturncode))
        if samtoolsReturncode != 0:
            raise StandardError('process "%s" returned %d' % (' '.join(samtoolsArgv), samtoolsReturncode))
    
    def setRepresentativeGenes(self):
        """Roughly equivalent to "distribute targets" in HybPiper."""
        self.representativePaftolTargetDict = {}
        for geneName in self.paftolTargetSet.paftolGeneDict:
            representativePaftolTarget = None
            maxMapqSum = None
            for organismName in self.paftolTargetSet.paftolGeneDict[geneName].paftolTargetDict:
                mapqSum = self.paftolTargetSet.paftolGeneDict[geneName].paftolTargetDict[organismName].mapqSum()
                if representativePaftolTarget is None or (mapqSum is not None and mapqSum > maxMapqSum):
                    representativePaftolTarget = self.paftolTargetSet.paftolGeneDict[geneName].paftolTargetDict[organismName]
                    maxMapqSum = mapqSum
            self.representativePaftolTargetDict[geneName] = representativePaftolTarget
            if representativePaftolTarget is None:
                logger.debug('represenative for %s: none', geneName)
            else:
                logger.debug('representative for %s: %s', representativePaftolTarget.paftolGene.name, representativePaftolTarget.organism.name)
    
    def distributeSingle(self):
        fForward = open(self.forwardFastq, 'r')
        fqiForward = Bio.SeqIO.QualityIO.FastqGeneralIterator(fForward)
        for fwdReadTitle, fwdReadSeq, fwdReadQual in fqiForward:
            readName = fwdReadTitle.split()[0]
            for paftolGene in self.paftolTargetSet.paftolGeneDict.values():
                if readName in paftolGene.qnameSet():
                    f = open(self.makeGeneFname(paftolGene.name, True), 'a')
                    logger.debug('appending to %s', f.name)
                    f.write('>%s\n%s\n' % (fwdReadTitle, fwdReadSeq))
                    f.close()
        fForward.close()
    
    def distributePaired(self):
        # FIXME: consider try...finally to ensure files are closed
        fForward = open(self.forwardFastq, 'r')
        fqiForward = Bio.SeqIO.QualityIO.FastqGeneralIterator(fForward)
        fReverse = open(self.reverseFastq, 'r')
        fqiReverse = Bio.SeqIO.QualityIO.FastqGeneralIterator(fReverse)
        for fwdReadTitle, fwdReadSeq, fwdReadQual in fqiForward:
            readName = fwdReadTitle.split()[0]
            # FIXME: premature end of reverse fastq will trigger
            # StopIteration and premature end of forward will leave
            # rest of reverse ignored
            revReadTitle, revReadSeq, revReadQual = fqiReverse.next()
            if readName != revReadTitle.split()[0]:
                raise StandardError('paired read files %s / %s out of sync at read %s / %s' % (self.forwardFastq, self.reverseFastq, fwdReadTitle, revReadTitle))
            for paftolGene in self.paftolTargetSet.paftolGeneDict.values():
                if readName in paftolGene.qnameSet():
                    f = open(self.makeGeneFname(paftolGene.name, True), 'a')
                    f.write('>%s\n%s\n' % (fwdReadTitle, fwdReadSeq))
                    f.write('>%s\n%s\n' % (revReadTitle, revReadSeq))
                    f.close()
        # FIXME: check for dangling stuff in reverse: should trigger
        # an exception:
        # revReadTitle, revReadSeq, revReadQual = fqiReverse.next()
        fForward.close()
        fReverse.close()
    
    def distribute(self):
        if self.isPaired():
            self.distributePaired()
        else:
            self.distributeSingle()
            
    def assembleSpadesParallel(self):
        """OBSOLETE -- Run SPAdes assemblies using GNU parallel, as the
original HybPiper implementation does.

Replaced by L{assembleGeneSpades} and no longer maintained / functional.
"""
        # consider --fg to ensure wait for all parallel processes?
        # is --eta really of any use here?
        # FIXME: hard-coded fasta pattern '{}_interleaved.fasta' for parallel
        if self.isPaired():
            spadesInputArgs = ['--12', '{}_interleaved.fasta']
        else:
            spadesInputArgs = ['-s', '{}.fasta']
        parallelSpadesArgv = ['parallel', 'spades.py', '--only-assembler', '--threads', '1', '--cov-cutoff', '%d' % self.spadesCovCutoff]
        if self.spadesKvalList is not None:
            parallelSpadesArgv.extend(['-k', ','.join(['%d' % k for k in self.spadesKvalList])])
        parallelSpadesArgv.extend(spadesInputArgs)
        parallelSpadesArgv.extend(['-o', '{}_spades'])
        # time parallel --eta spades.py --only-assembler --threads 1 --cov-cutoff 8 --12 {}/{}_interleaved.fasta -o {}/{}_spades :::: spades_genelist.txt > spades.log
        logger.debug('%s', ' '.join(parallelSpadesArgv))
        parallelSpadesProcess = subprocess.Popen(parallelSpadesArgv, stdin=subprocess.PIPE, cwd = self.makeWorkDirname())
        pid = os.fork()
        if pid == 0:
            for geneName in self.geneNameSet:
                parallelSpadesProcess.stdin.write('%s\n' % geneName)
            parallelSpadesProcess.stdin.close()
            os._exit(0)
        parallelSpadesProcess.stdin.close()
        wPid, wExit = os.waitpid(pid, 0)
        if pid != wPid:
            raise StandardError('wait returned pid %s (expected %d)' % (wPid, pid))
        if wExit != 0:
            raise StandardError('wait on forked process returned %d' % wExit)
        parallelSpadesReturncode = parallelSpadesProcess.wait()
        if parallelSpadesReturncode != 0:
            raise StandardError('parallel spades process exited with %d' % parallelSpadesReturncode)
        
    def makeGeneDirname(self, geneName):
        return 'spades-%s' % geneName
    
    def makeGeneDirPath(self, geneName):
        return os.path.join(self.makeWorkDirname(), self.makeGeneDirname(geneName))
            
    def assembleGeneSpades(self, geneName):
        # FIXME: should return file with contigs / scaffolds upon success, None otherwise
        # consider --fg to ensure wait for all parallel processes?
        # is --eta really of any use here?
        # FIXME: hard-coded fasta pattern '{}_interleaved.fasta' for parallel
        geneFname = self.makeGeneFname(geneName)
        if self.isPaired():
            spadesInputArgs = ['--12', geneFname]
        else:
            spadesInputArgs = ['-s', geneFname]
        if not os.path.exists(os.path.join(self.makeWorkDirname(), geneFname)):
            logger.debug('gene fasta file %s does not exist (no reads?)', geneFname)
            return None
        spadesArgv = ['spades.py', '--only-assembler', '--threads', '1', '--cov-cutoff', '%d' % self.spadesCovCutoff]
        if self.spadesKvalList is not None:
            spadesArgv.extend(['-k', ','.join(['%d' % k for k in self.spadesKvalList])])
        spadesArgv.extend(spadesInputArgs)
        spadesArgv.extend(['-o', self.makeGeneDirname(geneName)])
        logger.debug('%s', ' '.join(spadesArgv))
        spadesProcess = subprocess.Popen(spadesArgv, cwd = self.makeWorkDirname())
        spadesReturncode = spadesProcess.wait()
        if spadesReturncode != 0:
            # raise StandardError('spades process "%s" exited with %d' % (' '.join(spadesArgv), spadesReturncode))
            logger.warning('spades process "%s" exited with %d', ' '.join(spadesArgv), spadesReturncode)
        spadesContigFname = os.path.join(self.makeGeneDirPath(geneName), 'contigs.fasta')
        # logger.debug('spadesContigFname: %s', spadesContigFname)
        if os.path.exists(spadesContigFname):
            spadesContigList = list(Bio.SeqIO.parse(spadesContigFname, 'fasta'))
            # logger.debug('spadesContigFname: %s, %d contigs', spadesContigFname, len(spadesContigList))
        else:
            spadesContigList = None
            # logger.debug('spadesContigFname: %s, no contigs', spadesContigFname)
        return spadesContigList
    
    def translateGene(self, geneDna):
        # FIXME: add support for gene specific translation table setting
        l = len(geneDna) - (len(geneDna) % 3)
        if l < len(geneDna):
            logger.warning('gene %s: length %d is not an integer multiple of 3 -- not a CDS?', geneDna.id, len(geneDna))
        geneProtein = Bio.SeqRecord.SeqRecord(geneDna.seq[:l].translate(), id='%s-pep' % geneDna.id, description='%s, translated' % geneDna.description)
        return geneProtein
    
    def filterByPercentIdentity(self, exonerateResultList):
        return [e for e in exonerateResultList if e.percentIdentity >= self.exoneratePercentIdentityThreshold]
    
    def filterByContainment(self, exonerateResultList):
        
        def isContainedWithTiebreak(exonerateResult, other):
            if not other.containsQueryAlignmentRange(exonerateResult):
                return False
            if not exonerateResult.containsQueryAlignmentRange(other):
                return True
            # prefer shorter target alignment length (fewer gaps)
            if exonerateResult.targetAlignmentLength < other.targetAlignmentLength:
                return False
            elif exonerateResult.targetAlignmentLength > other.targetAlignmentLength:
                return True
            # subsequent tie breaking is arbitrary and intended to yield consistent results only
            # FIXME: resolving tie by arbitrarily preferring target start position
            if exonerateResult.targetAlignmentStart < other.targetAlignmentStart:
                return False
            elif exonerateResult.targetAlignmentStart > other.targetAlignmentStart:
                return True
            # FIXME: resolving tie using contig id, consider using more meaningful criteria but be mindful of biases...???
            if exonerateResult.targetId is None:
                raise StandardError('cannot break tie when exonerateResult.targetId is None')
            if other.targetId is None:
                raise StandardError('cannot break tie when other.targetId is None')
            if exonerateResult.targetId < other.targetId:
                return False
            elif other.targetId < exonerateResult.targetId:
                return True
            raise StandardError('cannot break tie: exonerateResult = %s, other = %s' % (str(exonerateResult), str(other)))
            
        nonContainedExonerateResultList = []
        for exonerateResult in exonerateResultList:
            isContained = False
            for other in exonerateResultList:
                isContained = isContained or ((exonerateResult is not other) and isContainedWithTiebreak(exonerateResult, other))
            if not isContained:
                nonContainedExonerateResultList.append(exonerateResult)
        return nonContainedExonerateResultList
 
    # query:   gattacatgactcga
    # contig1: gattacatga
    # contig2:      ca--actcga
    # trim contig2 because it has (more) gaps in the overlap region??
    # compute consensus -- along overlapping regions, or along entire query?
    def filterByOverlap(self, exonerateResultList):
        logger.warning('scanning for overlaps but not resolving them, pending development of concept')
        nonOverlappingExonerateResultList = []
        for exonerateResult in exonerateResultList:
            for other in exonerateResultList:
                if exonerateResult is not other:
                    if exonerateResult.overlapsQueryAlignmentRange(other):
                        logger.warning('overlap found, but not resolved: %s, %s', str(exonerateResult), str(other))
            nonOverlappingExonerateResultList.append(exonerateResult)
        return nonOverlappingExonerateResultList
    
    def filterExonerateResultList(self, geneName, exonerateResultList):
        logger.debug('gene %s: %d exonerate results', geneName, len(exonerateResultList))
        exonerateResultList = self.filterByPercentIdentity(exonerateResultList)
        logger.debug('gene %s: %d sufficiently close exonerate results', geneName, len(exonerateResultList))
        exonerateResultList = self.filterByContainment(exonerateResultList)
        logger.debug('gene %s: %d non-contained exonerate results', geneName, len(exonerateResultList))
        return exonerateResultList
    
    def reconstructCds(self, geneName):
        logger.debug('reconstructing CDS for gene %s', geneName)
        if self.representativePaftolTargetDict is None:
            raise StandardError('illegal state: no represesentative genes')
        if self.representativePaftolTargetDict[geneName] is None:
            raise StandardError('no representative for gene %s' % geneName)
        os.mkdir(self.makeGeneDirPath(geneName))
        contigList = self.assembleGeneSpades(geneName)
        if contigList is None:
            logger.warning('gene %s: no spades contigs', geneName)
            return None
        if len(contigList) == 0:
            logger.warning('gene %s: empty contig list', geneName)
            return None
        logger.debug('gene %s: %d spades contigs', geneName, len(contigList))
        geneProtein = self.translateGene(self.representativePaftolTargetDict[geneName].seqRecord)
        aminoAcidSet = set(Bio.Alphabet.IUPAC.protein.letters.lower())
        # allow stop translation
        aminoAcidSet.add('*')
        setDiff = set(str(geneProtein.seq).lower()) - aminoAcidSet
        if len(setDiff) > 0:
            logger.warning('gene %s: invalid amino acids %s' % (geneName, ', '.join(setDiff)))
            return None
        contigFname = os.path.join(self.makeGeneDirPath(geneName), '%s-contigs.fasta' % geneName)
        Bio.SeqIO.write(contigList, contigFname, 'fasta')
        exonerateRunner = paftol.tools.ExonerateRunner()
        exonerateResultList = exonerateRunner.parse(geneProtein, contigFname, 'protein2genome', len(contigList))
        logger.debug('%d contigs, %d exonerate results', len(contigList), len(exonerateResultList))
        if len(exonerateResultList) == 0:
            logger.warning('gene %s: no exonerate results from %d contigs', geneName, len(contigList))
        exonerateResultList.sort(cmpExonerateResultByQueryAlignmentStart)
        for exonerateResult in exonerateResultList:
            if exonerateResult.targetStrand == '-':
                exonerateResult.reverseComplementTarget()
        logger.warning('provisional filtering and supercontig construction, handling of overlapping contigs not finalised')
        filteredExonerateResultList = self.filterExonerateResultList(geneName, exonerateResultList)
        if len(filteredExonerateResultList) == 0:
            logger.warning('gene %s: no exonerate results left after filtering', geneName)
            return None
        supercontig = Bio.SeqRecord.SeqRecord(Bio.Seq.Seq(''.join([str(e.targetCdsSeq.seq) for e in filteredExonerateResultList])), id='%s_supercontig' % geneName)
        if len(supercontig) == 0:
            logger.warning('gene %s: empty supercontig', geneName)
            return None
        supercontigFname = os.path.join(self.makeGeneDirPath(geneName), '%s-supercontig.fasta' % geneName)
        Bio.SeqIO.write([supercontig], supercontigFname, 'fasta')
        supercontigErList = exonerateRunner.parse(geneProtein, supercontigFname, 'protein2genome', len(contigList))
        if len(supercontigErList) == 0:
            logger.warning('gene %s: no exonerate results from supercontig', geneName)
            return None
        # not filtering for percent identity to gene again, as that is already done
        if self.reverseFastq is not None:
            readsSpec = '%s, %s' % (self.forwardFastq, self.reverseFastq)
        else :
            readsSpec = self.forwardFastq
        splicedSupercontig = Bio.SeqRecord.SeqRecord(Bio.Seq.Seq(''.join([str(e.targetCdsSeq.seq) for e in supercontigErList])), id=geneName, description='reconstructed CDS computed by paftol.HybpiperAnalyser, targets: %s, reads: %s' % (self.targetsSourcePath, readsSpec))
        return splicedSupercontig
    
    # ideas for hybrid / consensus sequence for (multiple) re-mapping
    # reference CDS:     atgtac------catacagaagagacgtga
    # reconstructed CDS:    cactcatttcat---gga
    # "consensus"        atgCACTCAATTCAT   GGAgagacgtga
    # principe: Where reconstructed symbol is available, use that in preference.
    #   * gap in reference: use symbols from reconstructed (must be non-gap if pairwise alignment)
    #   * gap in reconstructed: skip symbols from reference
    #   * ends / portions with no alignment to reconstructed: fill in from reference
    # Problem: avoid non-homologous alignment portions (e.g. around borders of reconstructed)?

    def analyse(self):
        self.checkTargets()
        try:
            self.setup()
            self.mapReadsBwa()
            self.distribute()
            self.setRepresentativeGenes()
            reconstructedCdsDict = {}
            for geneName in self.paftolTargetSet.paftolGeneDict:
                reconstructedCdsDict[geneName] = self.reconstructCds(geneName)
            if self.statsCsvFilename is not None:
                csvFile = open(self.statsCsvFilename, 'w')
                csvDictWriter = PaftolTarget.makeCsvDictWriter(csvFile)
                for organismName in self.paftolTargetSet.organismDict:
                    for geneName in self.paftolTargetSet.organismDict[organismName].paftolTargetDict:
                        self.paftolTargetSet.organismDict[organismName].paftolTargetDict[geneName].writeCsvRow(csvDictWriter)
                csvFile.close()
            return reconstructedCdsDict
        finally:
            self.makeTgz()
            self.cleanup()
