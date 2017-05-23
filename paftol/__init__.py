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
    funnyCharRe = re.compile('[\t/ ;,$#]')
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


class DataFrame(object):
    
    def __init__(self, columnHeaderList):
        self.columnHeaderList = columnHeaderList[:]
        self.rowDictList = []
        
    def addRow(self, rowDict):
        if set(rowDict.keys()) != set(self.columnHeaderList):
            raise StandardError, 'key set %s is not compatible with column headers %s' % (', '.join([str(k) for k in rowDict.keys()]), ', '.join(self.columnHeaderList))
        self.rowDictList.append(rowDict)
        
    def writeCsv(self, f):
        csvDictWriter = csv.DictWriter(f, self.columnHeaderList)
        csvDictWriter.writeheader()
        for rowDict in self.rowDictList:
            csvDictWriter.writerow(rowDict)
            

class BwaParams(object):
    """Hold parameters for C{bwa} and provide argument vectors on that basis.
    
@ivar numThreads: BWA number of threads (C{-t} option)
@type numThreads: C{int}
@ivar minSeedLength: BWA minimum seed length (C{-k} option)
@type minSeedLength: C{int}
@ivar scoreThreshold: BWA score threshold for recording reads as mapped (C{-T} option)
@type scoreThreshold: C{int}
@ivar reseedTrigger: BWA re-seed trigger (C{-r} option)
@type reseedTrigger: C{float}
"""
    def __init__(self, numThreads=None, minSeedLength=None, scoreThreshold=None, reseedTrigger=None):
        self.numThreads = numThreads
        self.minSeedLength = minSeedLength
        self.scoreThreshold = scoreThreshold
        self.reseedTrigger = reseedTrigger
        
    def indexReferenceArgv(self, referenceFname):
        return ['bwa', 'index', referenceFname]
        
    def mappingMemArgv(self, referenceFname, forwardReadsFname, reverseReadsFname=None):
        argv = ['bwa', 'mem', '-M']
        if self.minSeedLength is not None:
            argv.extend(['-k', '%d' % self.minSeedLength])
        if self.reseedTrigger is not None:
            argv.extend(['-r', '%f' % self.reseedTrigger])
        if self.scoreThreshold is not None:
            argv.extend(['-T', '%d' % self.scoreThreshold])
        if self.numThreads is not None:
            argv.extend(['-t', '%d' % self.numThreads])
        argv.append(referenceFname)
        argv.append(forwardReadsFname)
        if reverseReadsFname is not None:
            argv.append(reverseReadsFname)
        return argv
    
    def referenceIndexArgv(self, referenceFname):
        return ['bwa', 'index', referenceFname]


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
@ivar flag: SAM flag (C{FLAG})
@type flag: C{int}
@ivar pos: SAM mapping position (C{POS})
@type pos: C{int}
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
        self.flag = int(w[1])
        self.rname = w[2]
        self.pos = int(w[3])
        self.mapq = int(w[4])
        self.cigar = w[5]
        self.seq = w[9]
        
    def isMapped(self):
        return self.flag & 4 == 0
        
    def getMatchLength(self):
        e = self.expandedCigar()
        return e.count('M') + e.count('D')

    def getEndpos(self):
        return self.pos + self.getMatchLength()
    
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
    
    csvFieldNames = ['organism', 'gene', 'seqLength', 'numSamAlignments']
    
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
    
    def numSamAlignments(self):
        return len(self.samAlignmentList)
    
    def csvRowDict(self):
        d = {}
        d['organism'] = self.organism.name
        d['gene'] = self.paftolGene.name
        d['seqLength'] = len(self.seqRecord)
        d['numSamAlignments'] = self.numSamAlignments()
        return d


class Organism(object):
    """Represent an organism (in the GenBank / NCBI sense of the term).
    
@ivar name: this organism's name
@type name: C{str}
@ivar paftolTargetDict: dictionary of genes in this organism
@type paftolTargetDict: C{dict} of C{PaftolTarget} instances with PAFTOL gene names as keys
"""
    
    csvFieldNames = ['organism', 'numGenes', 'numSamAlignments']
    
    def __init__(self, name):
        self.name = name
        self.paftolTargetDict = {}
        
    def numSamAlignments(self):
        return sum([len(t.samAlignmentList) for t in self.paftolTargetDict.values()])

    def csvRowDict(self):
        d = {}
        d['organism'] = self.name
        d['numGenes'] = len(self.paftolTargetDict)
        d['numSamAlignments'] = self.numSamAlignments()
        return d

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
    
    csvFieldNames = ['gene', 'numOrganisms', 'meanSeqLength', 'numSamAlignments']
    
    def __init__(self, name):
        self.name = name
        self.paftolTargetDict = {}

    def qnameSet(self):
        s = set()
        for paftolTarget in self.paftolTargetDict.values():
            s = s | paftolTarget.qnameSet()
        return s
    
    def meanSequenceLength(self):
        if len(self.paftolTargetDict) == 0:
            return None
        else:
            return float(sum([len(t.seqRecord) for t in self.paftolTargetDict.values()])) / float(len(self.paftolTargetDict))

    def numSamAlignments(self):
        return sum([len(t.samAlignmentList) for t in self.paftolTargetDict.values()])

    def csvRowDict(self):
        d = {}
        d['gene'] = self.name
        d['numOrganisms'] = len(self.paftolTargetDict)
        d['meanSeqLength'] = self.meanSequenceLength()
        d['numSamAlignments'] = self.numSamAlignments()
        return d
    

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
            
    def meanTargetLength(self, geneName):
        if geneName not in self.paftolGeneDict:
            raise StandardError, 'gene %s not contained in this target set'
        return self.paftolGeneDict[geneName].meanSequenceLength()
            
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
        
    def targetStats(self):
        dataFrame = DataFrame(PaftolTarget.csvFieldNames)
        for organism in self.organismDict.values():
            for paftolTarget in organism.paftolTargetDict.values():
                dataFrame.addRow(paftolTarget.csvRowDict())
        return dataFrame
    
    def geneStats(self):
        dataFrame = DataFrame(PaftolGene.csvFieldNames)
        for paftolGene in self.paftolGeneDict.values():
            dataFrame.addRow(paftolGene.csvRowDict())
        return dataFrame
    
    def organismStats(self):
        dataFrame = DataFrame(Organism.csvFieldNames)
        for organism in self.organismDict.values():
            dataFrame.addRow(organism.csvRowDict())
        return dataFrame
    
    def numSamAlignments(self):
        n = 0
        for organism in self.organismDict.values():
            for paftolTarget in organism.paftolTargetDict.values():
                n = n + paftolTarget.numSamAlignments()
        return n

    
class ReferenceGene(object):
    
    def __init__(self, geneId, referenceGenome, seqRecord, geneFeature, mrnaFeatureList=None, cdsFeatureList=None):
        self.geneId = geneId
        self.referenceGenome = referenceGenome
        self.seqRecord = seqRecord
        self.geneFeature = geneFeature
        self.mrnaFeatureList = [] if mrnaFeatureList is None else mrnaFeatureList[:]
        self.cdsFeatureList = [] if cdsFeatureList is None else cdsFeatureList[:]
        
    def getSequenceId(self):
        return self.seqRecord.id.split('.')[0]
        
    def containsHsp(self, hspAccession, hsp):
        if self.getSequenceId() != hspAccession:
            return False
        return self.geneFeature.location.start <= hsp.sbjct_start and self.geneFeature.location.end >= hsp.sbjct_end
    
    def getLength(self):
        return abs(self.geneFeature.location.end - self.geneFeature.location.start)
    
    def containsSamAlignment(self, samAlignment):
        if self.getSequenceId() != samAlignment.rname:
            return False
        return self.geneFeature.location.start <= samAlignment.pos and self.geneFeature.location.end >= samAlignment.getEndpos()
                
    def getGeneName(self):
        if 'name' in self.geneFeature.qualifiers:
            return self.geneFeature.qualifiers['name'][0]
        else:
            return None
                
    def getGeneNote(self):
        if 'note' in self.geneFeature.qualifiers:
            return self.geneFeature.qualifiers['note'][0]
        else:
            return None
        
    def getMrnaProduct(self):
        # CHECKME: returning 'product' qualifier value from feature with that qualifier -- may be more thorough to check that all are the same?
        for mrnaFeature in self.mrnaFeatureList:
            if 'product' in mrnaFeature.qualifiers:
                return mrnaFeature.qualifiers['product'][0]
        return None
        
    def getCdsProduct(self):
        # CHECKME: returning 'product' qualifier value from feature with that qualifier -- may be more thorough to check that all are the same?
        for cdsFeature in self.cdsFeatureList:
            if 'product' in cdsFeature.qualifiers:
                return cdsFeature.qualifiers['product'][0]
        return None


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
        self.genomeLength = None
        
    def scanGenesAth(self):
        if self.genbankFname is None:
            raise StandardError('no GenBank file name, cannot scan genes (ath method)')
        mrnaFeatureDict = {}
        cdsFeatureDict = {}
        self.geneList = []
        self.genomeLength = 0
        geneDict = {}
        with open(self.genbankFname, 'r') as f:
            for seqRecord in Bio.SeqIO.parse(f, 'genbank'):
                self.genomeLength = self.genomeLength + len(seqRecord)
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
                            geneDict[geneId].mrnaFeatureList.append(seqFeature)
                    elif seqFeature.type == 'CDS':
                        geneId = seqFeature.qualifiers['locus_tag'][0]
                        if geneId in geneDict:
                            geneDict[geneId].cdsFeatureList.append(seqFeature)

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
        # sys.stderr.flush()
        # sys.stdout.flush()
        blastnProcess = subprocess.Popen(blastnArgv, stdin=subprocess.PIPE, stdout = subprocess.PIPE)
        # subprocess.call(['lsof', '-p', '%d' % os.getpid()])
        # blastnProcess.stdin.flush()
        pid = os.fork()
        if pid == 0:
            # reload(Bio.SeqIO)
            blastnProcess.stdout.close()
            # paftolTargetSet.writeFasta(sys.stderr)
            # srList = paftolTargetSet.getSeqRecordList()
            # sys.stderr.write('target set has %d seqRecords\n' % len(srList))
            # sr = Bio.SeqRecord.SeqRecord(Bio.Seq.Seq('A'), id = 'srDummy', description = '')
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
            # paftolTargetSet.writeFasta(blastnProcess.stdin)
            # FIXME: generating FASTA string and writing that manually to work around unresolved broken pipe issue
            for sr in paftolTargetSet.getSeqRecordList():
                blastnProcess.stdin.write(sr.format('fasta'))
            blastnProcess.stdin.close()
            os._exit(0)
        blastnProcess.stdin.close()
        # dict solely serves to check for duplicate BLAST records
        targetIdToGeneDict = {}
        targetGeneTable = DataFrame(['targetId', 'geneId', 'geneName', 'geneNote', 'mrnaProduct', 'cdsProduct'])
        for blastRecord in Bio.Blast.NCBIXML.parse(blastnProcess.stdout):
            targetId = blastRecord.query
            if targetId in targetIdToGeneDict:
                raise StandardError('duplicate BLAST record for target %s' % targetId)
            geneList = []
            for blastAlignment in blastRecord.alignments:
                for hsp in blastAlignment.hsps:
                    for gene in self.findGenesByHsp(blastAlignment.accession, hsp):
                        if gene not in geneList:
                            dfRow = {}
                            dfRow['targetId'] = targetId
                            dfRow['geneId'] = gene.geneId
                            dfRow['geneName'] = gene.getGeneName()
                            dfRow['geneNote'] = gene.getGeneNote()
                            dfRow['mrnaProduct'] = gene.getMrnaProduct()
                            dfRow['cdsProduct'] = gene.getCdsProduct()
                            targetGeneTable.addRow(dfRow)
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
        return targetGeneTable
    
    def findGeneIdForSamAlignment(self, samAlignment):
        # FIXME: clumsy linear search
        for gene in self.geneList:
            if gene.containsSamAlignment(samAlignment):
                return gene.geneId

    def mapReadsStatsBwaMem(self, forwardReadsFname, reverseReadsFname=None, bwaParams=None):
        if bwaParams is None:
            bwaParams = BwaParams()
        bwaArgv = bwaParams.mappingMemArgv(self.fastaFname, forwardReadsFname, reverseReadsFname)
        logger.debug('%s', ' '.join(bwaArgv))
        bwaProcess = subprocess.Popen(bwaArgv, stdout=subprocess.PIPE)
        dfRowDict = {}
        intergenicLength = self.genomeLength
        for gene in self.geneList:
            geneLength = gene.getLength()
            dfRowDict[gene.geneId] = {'geneId': gene.geneId, 'geneLength': geneLength, 'numHits': 0}
            if intergenicLength is not None:
                intergenicLength = intergenicLength - geneLength
        intergenicId = 'intergenic'
        dfRowDict[intergenicId] = {'geneId': intergenicId, 'geneLength': intergenicLength, 'numHits': 0}
        unmappedId = 'unmapped'
        dfRowDict[unmappedId] = {'geneId': unmappedId, 'geneLength': None, 'numHits': 0}
        rawmapTable = DataFrame(['qname', 'rname', 'pos'])
        samLine = bwaProcess.stdout.readline()
        while samLine != '':
            # logger.debug(samLine)
            if samLine[0] != '@':
                samAlignment = SamAlignment(samLine)
                if samAlignment.isMapped():
                    rawmapTable.addRow({'qname': samAlignment.qname, 'rname': samAlignment.rname, 'pos': samAlignment.pos})
                    geneId = self.findGeneIdForSamAlignment(samAlignment)
                    if geneId is None:
                        geneId = intergenicId
                else:
                    geneId = unmappedId
                dfRowDict[geneId]['numHits'] = dfRowDict[geneId]['numHits'] + 1
            samLine = bwaProcess.stdout.readline()
        bwaProcess.stdout.close()
        bwaReturncode = bwaProcess.wait()
        if bwaReturncode != 0:
            raise StandardError('process "%s" returned %d' % (' '.join(bwaArgv), bwaReturncode))
        statsTable = DataFrame(['geneId', 'geneLength', 'numHits'])
        for gene in self.geneList:
            statsTable.addRow(dfRowDict[gene.geneId])
        statsTable.addRow(dfRowDict[intergenicId])
        statsTable.addRow(dfRowDict[unmappedId])
        return statsTable, rawmapTable

        
class HybpiperAnalyser(HybseqAnalyser):
    """L{HybseqAnalyser} subclass that implements an analysis process
close to the HybPiper pipeline.

Some parameters to BWA and SPAdes can be controlled via instance
variables as documented below. Defaults of these parameters correspond
to the defaults provided by BWA and SPAdes, respectively (at the time
of developing this).

@ivar spadesCovCutoff: SPAdes coverage cutoff (C{--cov-cutoff} option)
@type spadesCovCutoff: C{int}
@ivar spadesKvalList: SPAdes oligomer length value list (C{-k} option)
@type spadesKvalList: C{list} of C{int}, or C{None}
"""
    
    def __init__(self, targetsSourcePath, forwardFastq, reverseFastq=None, workdirTgz=None, workDirname='pafpipertmp', bwaParams=None):
        super(HybpiperAnalyser, self).__init__(targetsSourcePath, forwardFastq, reverseFastq, workdirTgz, workDirname)
        if bwaParams is None:
            self.bwaParams = BwaParams()
        else:
            self.bwaParams = bwaParams
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
        bwaIndexArgv = self.bwaParams.indexReferenceArgv(referenceFname)
        logger.debug('%s', ' '.join(bwaIndexArgv))
        subprocess.check_call(bwaIndexArgv)
        
    def mapReadsBwa(self):
        """Map reads to gene sequences (from multiple organisms possibly).
"""
        logger.debug('mapping reads to gene sequences')
        referenceFname = self.makeTargetsFname(True)
        self.bwaIndexReference(referenceFname)
        forwardReadsFname = os.path.join(os.getcwd(), self.forwardFastq)
        if self.reverseFastq is None:
            reverseReadsFname = None
        else:
            reverseReadsFname = os.path.join(os.getcwd(), self.reverseFastq)
        bwaArgv = self.bwaParams.mappingMemArgv(referenceFname, forwardReadsFname, reverseReadsFname)
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
                tStats = self.paftolTargetSets.targetStats()
                with open(self.statsCsvFilename, 'w') as csvFile:
                    tStats.writeCsv(csvFile)
                csvFile.close()
            return reconstructedCdsDict
        finally:
            self.makeTgz()
            self.cleanup()
