import sys
import re
import os
import tempfile
import subprocess
import shutil
import multiprocessing
import logging

import Bio
import Bio.SeqIO
import Bio.SeqIO.QualityIO
import Bio.Alphabet.IUPAC

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
    
Instances of this class take a FASTA file of target locus sequences
and FASTQ files (one or two, for single / paired end, respectively),
and provide methods for running analyses to reconstruct sequences of
the target loci.
"""

    def __init__(self, targetsSourcePath, forwardFastq, reverseFastq=None, workdirTgz=None, workDirname='paftoolstmp'):
        self.targetsSourcePath = targetsSourcePath
        self.forwardFastq = forwardFastq
        self.reverseFastq = reverseFastq
        self.workdirTgz = workdirTgz
        self.workDirname = workDirname
        self.tmpDirname = None
        # parameters for ensuring file names don't clash, e.g. because locus / organism name is same as targets basename etc.
        self.targetsFname = 'targets.fasta'
        self.locusFnamePattern = 'locus-%s.fasta'
        self.allowInvalidBases = False
    
    def __str__(self):
        return 'HybseqAnalyser(targetsSourcePath=%s, forwardFastq=%s, reverseFastq=%s)' % (repr(self.targetsSourcePath), repr(self.forwardFastq), repr(self.reverseFastq))
        
    def checkTargets(self):
        # FIXME: merge with __init__()?
        for targetSr in Bio.SeqIO.parse(self.targetsSourcePath, 'fasta', alphabet = Bio.Alphabet.IUPAC.ambiguous_dna):
            if not self.allowInvalidBases:
                setDiff = set(str(targetSr.seq).lower()) - set('acgt')
                if len(setDiff) != 0:
                    raise StandardError, 'target %s: illegal base(s) %s' % (targetSr.id, ', '.join(setDiff))

    def isPaired(self):
        return self.reverseFastq is not None
    
    def analyse(self):
        raise StandardError, 'not implemented in this "abstract" base class'
    
    def setupTmpdir(self):
        if self.tmpDirname is not None:
            raise StandardError, 'illegal state: already have generated working directory %s' % self.tmpDirname
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
            raise StandardError, 'illegal state: no temporary directory and hence no working directory'
        # logger.debug('tmpDirname = %s, workDirname = %s', self.tmpDirname, self.workDirname)
        return os.path.join(self.tmpDirname, self.workDirname)
    
    def makeTargetsFname(self, absolutePath=False):
        if absolutePath:
            return os.path.join(self.makeWorkDirname(), self.targetsFname)
        else:
            return self.targetsFname

    def makeLocusFname(self, locusName, absolutePath=False):
        locusFname = self.locusFnamePattern % locusName
        if absolutePath:
            return os.path.join(self.makeWorkDirname(), locusFname)
        else:
            return locusFname
        
    def makeTgz(self):
        if self.workdirTgz is not None:
            if self.tmpDirname is None:
                raise StandardError, 'illegal state: no temporary directory generated'
            tmpTgz = os.path.join(self.tmpDirname, '%s.tgz' % self.workDirname)
            tgzArgv = ['tar', '-zcf', tmpTgz, self.workDirname]
            tgzProcess = subprocess.Popen(tgzArgv, cwd = self.tmpDirname)
            tgzReturncode = tgzProcess.wait()
            if tgzReturncode != 0:
                raise StandardError, 'process "%s" returned %d' % (' '.join(tgzArgv), tgzReturncode)
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
    
    def __init__(self, samLine):
        if samLine[-1] == '\n':
            samLine = samLine[:-1]
        w = samLine.split('\t')
        self.qname = w[0]
        self.rname = w[2]
        self.mapq = int(w[4])
        self.cigar = w[5]
        self.seq = w[9]


class OrganismLocus(object):
    """Represent a locus in an organism.
    
The main content of instances of this class is a C{SeqRecord}
containing the sequence of the locus in the organism, thus
facilitating handling of multiple loci and multiple organisms.

@ivar organism: the organism
@type organism: C{Organism}
@ivar locus: the locus 
@type locus: C{Locus}
@ivar seqRecord: the sequence of this organism at this locus
@type seqRecord: C{Bio.SeqRecord.SeqRecord}
"""
    
    def __init__(self, organism, locus, seqRecord):
        self.organism = organism
        self.locus = locus
        self.seqRecord = seqRecord
        self.samAlignmentList = []
        if locus.name in organism.organismLocusDict or organism.name in locus.organismLocusDict:
            raise StandardError, 'duplicate organism/locus: organism = %s, locus = %s, seqId = %s' % (organism.name, locus.name, seqRecord.id)
        organism.organismLocusDict[locus.name] = self
        locus.organismLocusDict[organism.name] = self
        
    def addSamAlignment(self, samAlignment):
        self.samAlignmentList.append(samAlignment)
            
    def mapqSum(self):
        if len(self.samAlignmentList) == 0:
            return None
        return sum([a.mapq for a in self.samAlignmentList])
    
    def qnameSet(self):
        # FIXME: may have to trim away "/1", "/2"?
        return set([a.qname for a in self.samAlignmentList])


class Organism(object):
    """Represent an organism (in the GenBank / NCBI sense of the term).
    
@ivar name: this organism's name
@type name: C{str}
@ivar organismLocusDict: dictionary of loci in this organism
@type organismLocusDict: C{dict} of C{OrganismLocus} instances with locus names as keys
"""
    
    def __init__(self, name):
        self.name = name
        self.organismLocusDict = {}


class Locus(object):
    """Represent a locus.

@ivar name: the name of this locus
@type name: C{str}
@ivar organismLocusDict: dictionary of organisms with this locus
@type organismLocusDict: C{dict} of C{OrganismLocus} instances with organism names as keys
"""
    
    def __init__(self, name):
        self.name = name
        self.organismLocusDict = {}

    def qnameSet(self):
        s = set()
        for organismLocus in self.organismLocusDict.values():
            s = s | organismLocus.qnameSet()
        return s
    

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

    organismLocusRe = re.compile('([^-]+)-([^-]+)')
    
    def __init__(self, targetsSourcePath, forwardFastq, reverseFastq=None, workdirTgz=None, workDirname='pafpipertmp'):
        super(HybpiperAnalyser, self).__init__(targetsSourcePath, forwardFastq, reverseFastq, workdirTgz, workDirname)
        self.bwaNumThreads = 1
        self.bwaMinSeedLength = 19
        self.bwaScoreThreshold = 30
        self.bwaReseedTrigger = 1.5
        self.spadesCovCutoff = 8
        self.spadesKvalList = None
        self.exoneratePercentIdentityThreshold = 65.0
        self.initOrganismLocusDicts()
        
    def extractOrganismAndLocusNames(self, s):
        m = self.organismLocusRe.match(s)
        if m is not None:
            organismName = m.group(1)
            locusName = m.group(2)
        else:
            organismName = 'unknown'
            locusName = s
        return organismName, locusName
        
    def initOrganismLocusDicts(self):
        if self.targetsSourcePath is None:
            raise StandardError, 'illegal state: cannot init organism and locus dicts with targetsSourcePath = None'
        self.locusDict = {}
        self.organismDict = {}
        for sr in Bio.SeqIO.parse(self.targetsSourcePath, 'fasta'):
            organismName, locusName = self.extractOrganismAndLocusNames(sr.id)
            if not isSane(organismName):
                raise StandardError, 'bad organism name: %s' % organismName
            if not isSane(locusName):
                raise StandardError, 'bad locus name: %s' % locusName
            if organismName not in self.organismDict:
                self.organismDict[organismName] = Organism(organismName)
            if locusName not in self.locusDict:
                self.locusDict[locusName] = Locus(locusName)
            organismLocus = OrganismLocus(self.organismDict[organismName], self.locusDict[locusName], sr)
        logger.info('%s organisms, %s loci' % (len(self.organismDict), len(self.locusDict)))
        self.representativeOrganismLocusDict = None

    def setup(self):
        logger.debug('setting up')
        if self.targetsSourcePath is None:
            raise StandardError, 'illegal state: cannot set up with targetsSourcePath = None'
        self.setupTmpdir()
        shutil.copy(self.targetsSourcePath, self.makeTargetsFname(True))
            
    def cleanup(self):
        self.cleanupTmpdir()
    
    def bwaIndexReference(self):
        bwaIndexArgv = ['bwa', 'index', self.makeTargetsFname(True)]
        logger.debug('%s', ' '.join(bwaIndexArgv))
        subprocess.check_call(bwaIndexArgv)
        
    def mapReadsBwa(self):
        """Map reads to locus sequences (from multiple organisms possibly).
"""
        logger.debug('mapping reads to locus sequences')
        self.bwaIndexReference()
        fastqArgs = [os.path.join(os.getcwd(), self.forwardFastq)]
        if self.reverseFastq is not None:
            fastqArgs.append(os.path.join(os.getcwd(), self.reverseFastq))
        # bwa parameters for tweaking considerations: -k, -r, -T
        bwaArgv = ['bwa', 'mem', '-M', '-k', '%d' % self.bwaMinSeedLength, '-T', '%d' % self.bwaScoreThreshold, '-r', '%f' % self.bwaReseedTrigger, '-t', '%d' % self.bwaNumThreads, self.makeTargetsFname()] + fastqArgs
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
                organismName, locusName = self.extractOrganismAndLocusNames(samAlignment.rname)
                if organismName not in self.organismDict:
                    raise StandardError, 'unknown organism: %s' % organismName
                if locusName not in self.locusDict:
                    raise standardError, 'unknown locus: %s' % locusName
                if locusName not in self.organismDict[organismName].organismLocusDict:
                    raise StandardError, 'no entry for locus %s in organism %s' % (locusName, organismName)
                organismLocus = self.organismDict[organismName].organismLocusDict[locusName]
                organismLocus.addSamAlignment(samAlignment)
            samLine = samtoolsProcess.stdout.readline()
        bwaProcess.stdout.close()
        samtoolsProcess.stdout.close()
        bwaReturncode = bwaProcess.wait()
        samtoolsReturncode = samtoolsProcess.wait()
        if bwaReturncode != 0:
            raise StandardError, 'process "%s" returned %d' % (' '.join(bwaArgv), bwaReturncode)
        if samtoolsReturncode != 0:
            raise StandardError, 'process "%s" returned %d' % (' '.join(samtoolsArgv), samtoolsReturncode)
    
    def setRepresentativeLoci(self):
        """Roughly equivalent to "distribute targets" in HybPiper."""
        self.representativeOrganismLocusDict = {}
        for locusName in self.locusDict:
            representativeOrganismLocus = None
            maxMapqSum = None
            for organismName in self.locusDict[locusName].organismLocusDict:
                mapqSum = self.locusDict[locusName].organismLocusDict[organismName].mapqSum()
                if representativeOrganismLocus is None or (mapqSum is not None and mapqSum > maxMapqSum):
                    representativeOrganismLocus = self.locusDict[locusName].organismLocusDict[organismName]
                    maxMapqSum = mapqSum
            self.representativeOrganismLocusDict[locusName] = representativeOrganismLocus
            if representativeOrganismLocus is None:
                logger.debug('represenative for %s: none', locusName)
            else:
                logger.debug('representative for %s: %s', representativeOrganismLocus.locus.name, representativeOrganismLocus.organism.name)
    
    def distributeSingle(self):
        fForward = open(self.forwardFastq, 'r')
        fqiForward = Bio.SeqIO.QualityIO.FastqGeneralIterator(fForward)
        for fwdReadTitle, fwdReadSeq, fwdReadQual in fqiForward:
            readName = fwdReadTitle.split()[0]
            for locus in self.locusDict.values():
                if readName in locus.qnameSet():
                    f = open(self.makeLocusFname(locus.name, True), 'a')
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
                raise StandardError, 'paired read files %s / %s out of sync at read %s / %s' % (self.forwardFastq, self.reverseFastq, fwdReadTitle, revReadTitle)
            for locus in self.locusDict.values():
                if readName in locus.qnameSet():
                    f = open(self.makeLocusFname(locus.name, True), 'a')
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
            for locusName in self.locusNameSet:
                parallelSpadesProcess.stdin.write('%s\n' % locusName)
            parallelSpadesProcess.stdin.close()
            os._exit(0)
        parallelSpadesProcess.stdin.close()
        wPid, wExit = os.waitpid(pid, 0)
        if pid != wPid:
            raise StandardError, 'wait returned pid %s (expected %d)' % (wPid, pid)
        if wExit != 0:
            raise StandardError, 'wait on forked process returned %d' % wExit
        parallelSpadesReturncode = parallelSpadesProcess.wait()
        if parallelSpadesReturncode != 0:
            raise StandardError, 'parallel spades process exited with %d' % parallelSpadesReturncode
        
    def makeLocusDirname(self, locusName):
        return 'spades-%s' % locusName
    
    def makeLocusDirPath(self, locusName):
        return os.path.join(self.makeWorkDirname(), self.makeLocusDirname(locusName))
            
    def assembleLocusSpades(self, locusName):
        # FIXME: should return file with contigs / scaffolds upon success, None otherwise
        # consider --fg to ensure wait for all parallel processes?
        # is --eta really of any use here?
        # FIXME: hard-coded fasta pattern '{}_interleaved.fasta' for parallel
        locusFname = self.makeLocusFname(locusName)
        if self.isPaired():
            spadesInputArgs = ['--12', locusFname]
        else:
            spadesInputArgs = ['-s', locusFname]
        if not os.path.exists(os.path.join(self.makeWorkDirname(), locusFname)):
            logger.debug('locus fasta file %s does not exist (no reads?)', locusFname)
            return None
        spadesArgv = ['spades.py', '--only-assembler', '--threads', '1', '--cov-cutoff', '%d' % self.spadesCovCutoff]
        if self.spadesKvalList is not None:
            spadesArgv.extend(['-k', ','.join(['%d' % k for k in self.spadesKvalList])])
        spadesArgv.extend(spadesInputArgs)
        spadesArgv.extend(['-o', self.makeLocusDirname(locusName)])
        logger.debug('%s', ' '.join(spadesArgv))
        spadesProcess = subprocess.Popen(spadesArgv, cwd = self.makeWorkDirname())
        spadesReturncode = spadesProcess.wait()
        if spadesReturncode != 0:
            # raise StandardError, 'spades process "%s" exited with %d' % (' '.join(spadesArgv), spadesReturncode)
            logger.warning('spades process "%s" exited with %d', ' '.join(spadesArgv), spadesReturncode)
        spadesContigFname = os.path.join(self.makeLocusDirPath(locusName), 'contigs.fasta')
        # logger.debug('spadesContigFname: %s', spadesContigFname)
        if os.path.exists(spadesContigFname):
            spadesContigList = list(Bio.SeqIO.parse(spadesContigFname, 'fasta'))
            # logger.debug('spadesContigFname: %s, %d contigs', spadesContigFname, len(spadesContigList))
        else:
            spadesContigList = None
            # logger.debug('spadesContigFname: %s, no contigs', spadesContigFname)
        return spadesContigList
    
    def translateLocus(self, locusDna):
        # FIXME: add support for locus specific translation table setting
        l = len(locusDna) - (len(locusDna) % 3)
        if l < len(locusDna):
            logger.warning('locus %s: length %d is not an integer multiple of 3 -- not a CDS?', locusDna.id, len(locusDna))
        locusProtein = Bio.SeqRecord.SeqRecord(locusDna.seq[:l].translate(), id='%s-pep' % locusDna.id, description='%s, translated' % locusDna.description)
        return locusProtein
    
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
            # FIXME: resolving tie by arbitrarily preferring target start position
            if exonerateResult.targetAlignmentStart < other.targetAlignmentStart:
                return False
            # FIXME: resolving tie using contig id, consider using more meaningful criteria but be mindful of biases...???
            if exonerateResult.targetId is None:
                raise StandardError, 'cannot break tie when exonerateResult.targetId is None'
            if other.targetId is None:
                raise StandardError, 'cannot break tie when other.targetId is None'
            if exonerateResult.targetId < other.targetId:
                return False
            elif other.targetId < exonerateResult.targetId:
                return True
            raise StandardError, 'cannot break tie: exonerateResult = %s, other = %s' % (str(exonerateResult), str(other))
            
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
        nonOverlappingExonerateResultList = []
        for exonerateResult in exonerateResultList:
            for other in exonerateResultList:
                if exonerateResult is not other:
                    if exonerateResult.overlapsQueryAlignmentRange(other):
                        logger.debug('overlap found: %s, %s', str(exonerateResult), str(other))
            nonOverlappingExonerateResultList.append(exonerateResult)
        return nonOverlappingExonerateResultList
    
    def filterExonerateResultList(self, locusName, exonerateResultList):
        logger.debug('locus %s: %d exonerate results', locusName, len(exonerateResultList))
        exonerateResultList = self.filterByPercentIdentity(exonerateResultList)
        logger.debug('locus %s: %d sufficiently close exonerate results', locusName, len(exonerateResultList))
        exonerateResultList = self.filterByContainment(exonerateResultList)
        logger.debug('locus %s: %d non-contained exonerate results', locusName, len(exonerateResultList))
        return exonerateResultList
    
    def reconstructCds(self, locusName):
        logger.debug('reconstructing CDS for locus %s', locusName)
        if self.representativeOrganismLocusDict is None:
            raise StandardError, 'illegal state: no represesentative loci'
        if self.representativeOrganismLocusDict[locusName] is None:
            raise StandardError, 'no representative for locus %s' % locusName
        os.mkdir(self.makeLocusDirPath(locusName))
        contigList = self.assembleLocusSpades(locusName)
        if contigList is None:
            logger.warning('locus %s: no spades contigs', locusName)
            return None
        if len(contigList) == 0:
            logger.warning('locus %s: empty contig list', locusName)
            return None
        logger.debug('locus %s: %d spades contigs', locusName, len(contigList))
        locusProtein = self.translateLocus(self.representativeOrganismLocusDict[locusName].seqRecord)
        setDiff = set(str(locusProtein.seq).lower()) - set(Bio.Alphabet.IUPAC.protein.letters.lower())
        if len(setDiff) > 0:
            logger.warning('locus %s: invalid amino acids %s' % (locusName, ', '.join(setDiff)))
            return None
        contigFname = os.path.join(self.makeLocusDirPath(locusName), '%s-contigs.fasta' % locusName)
        Bio.SeqIO.write(contigList, contigFname, 'fasta')
        exonerateRunner = paftol.tools.ExonerateRunner()
        exonerateResultList = exonerateRunner.parse(locusProtein, contigFname, 'protein2genome', len(contigList))
        logger.debug('%d contigs, %d exonerate results', len(contigList), len(exonerateResultList))
        if len(exonerateResultList) == 0:
            logger.warning('locus %s: no exonerate results from %d contigs', locusName, len(contigList))
        exonerateResultList.sort(cmpExonerateResultByQueryAlignmentStart)
        for exonerateResult in exonerateResultList:
            if exonerateResult.targetStrand == '-':
                exonerateResult.reverseComplementTarget()
        filteredExonerateResultList = self.filterExonerateResultList(locusName, exonerateResultList)
        if len(filteredExonerateResultList) == 0:
            logger.warning('locus %s: no exonerate results left after filtering', locusName)
            return None
        supercontig = Bio.SeqRecord.SeqRecord(Bio.Seq.Seq(''.join([str(e.targetCdsSeq.seq) for e in filteredExonerateResultList])), id='%s_supercontig' % locusName)
        if len(supercontig) == 0:
            logger.warning('locus %s: empty supercontig', locusName)
            return None
        supercontigFname = os.path.join(self.makeLocusDirPath(locusName), '%s-supercontig.fasta' % locusName)
        Bio.SeqIO.write([supercontig], supercontigFname, 'fasta')
        supercontigErList = exonerateRunner.parse(locusProtein, supercontigFname, 'protein2genome', len(contigList))
        if len(supercontigErList) == 0:
            logger.warning('locus %s: no exonerate results from supercontig', locusName)
            return None
        # not filtering for percent identity to locus again, as that is already done
        if self.reverseFastq is not None:
            readsSpec = '%s, %s' % (self.forwardFastq, self.reverseFastq)
        else :
            readsSpec = self.forwardFastq
        splicedSupercontig = Bio.SeqRecord.SeqRecord(Bio.Seq.Seq(''.join([str(e.targetCdsSeq.seq) for e in supercontigErList])), id=locusName, description='reconstructed CDS computed by paftol.HybpiperAnalyser, targets: %s, reads: %s' % (self.targetsSourcePath, readsSpec))
        return splicedSupercontig

    def analyse(self):
        self.checkTargets()
        try:
            self.setup()
            self.mapReadsBwa()
            self.distribute()
            self.setRepresentativeLoci()
            reconstructedCdsDict = {}
            for locusName in self.locusDict:
                reconstructedCdsDict[locusName] = self.reconstructCds(locusName)
            return reconstructedCdsDict
        finally:
            self.makeTgz()
            self.cleanup()
