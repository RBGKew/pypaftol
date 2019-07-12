import sys
import re
import os
import os.path
import copy
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
import paftol.version


__version__ = paftol.version.__version__


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

    def __init__(self, workdirTgz=None, workDirname='paftoolstmp'):
        """Initialiser

@param workdirTgz: name of the tgz archive of the working directory
@type workdirTgz: C{String}
@param workDirname: name of the working directory
@type workDirname: C{String}
"""
        # FIXME: temporary workdir management should be moved out of analyser -- perhaps into result?
        self.workdirTgz = workdirTgz
        self.workDirname = workDirname
        self.tmpDirname = None
        # parameters for ensuring file names don't clash, e.g. because paftolGene / organism name is same as targets basename etc.
        self.forwardFasta = 'fwd.fasta'
        self.reverseFasta = 'rev.fasta'
        self.targetsFname = 'targets.fasta'
        self.geneAssemblyDirnamePattern = 'targetasm-%s'
        self.geneReadFnamePattern = 'gene-%s.fasta'
        self.geneContigsFnamePattern = '%s-contigs.fasta'
        self.geneRepresentativeFnamePattern = 'generep-%s.fasta'
        self.exoneratePercentIdentityThreshold = 65.0

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

    def cleanup(self):
        self.cleanupTmpdir()

    def makeGeneDirname(self, geneName):
        return self.geneAssemblyDirnamePattern % geneName

    def makeGeneDirPath(self, geneName):
        return self.makeWorkdirPath(self.makeGeneDirname(geneName))

    def makeWorkDirname(self):
        if self.tmpDirname is None:
            raise StandardError('illegal state: no temporary directory and hence no working directory')
        # logger.debug('tmpDirname = %s, workDirname = %s', self.tmpDirname, self.workDirname)
        return os.path.join(self.tmpDirname, self.workDirname)

    def makeWorkdirPath(self, filename):
        return os.path.join(self.makeWorkDirname(), filename)
    
    def makeGeneContigsFname(self, geneName):
        return self.makeWorkdirPath(self.geneContigsFnamePattern % geneName)

    def makeTargetsFname(self, absolutePath=False):
        if absolutePath:
            return self.makeWorkdirPath(self.targetsFname)
        else:
            return self.targetsFname

    def makeGeneReadFname(self, geneName, absolutePath=False):
        # logger.warning('obsolescent -- responsibility for read fnames moving to TargetAssembler')
        geneReadFname = self.geneReadFnamePattern % geneName
        if absolutePath:
            return self.makeWorkdirPath(geneReadFname)
        else:
            return geneReadFname

    def makeGeneRepresentativeFname(self, geneName, absolutePath=False):
        geneRepresentativeFname = self.geneRepresentativeFnamePattern % geneName
        if absolutePath:
            return self.makeWorkdirPath(geneRepresentativeFname)
        else:
            return geneRepresentativeFname

    def translateGene(self, geneDna):
        # FIXME: add support for gene specific translation table setting
        # FIXME: not really an instance method -- static?
        l = len(geneDna) - (len(geneDna) % 3)
        if l < len(geneDna):
            logger.warning('gene %s: length %d is not an integer multiple of 3 -- not a CDS?', geneDna.id, len(geneDna))
        geneProtein = Bio.SeqRecord.SeqRecord(geneDna.seq[:l].translate(), id='%s-pep' % geneDna.id, description='%s, translated' % geneDna.description)
        return geneProtein

    def setRepresentativeGenes(self, result):
        """Roughly equivalent to "distribute targets" in HybPiper."""
        result.representativePaftolTargetDict = {}
        for geneName in result.paftolTargetSet.paftolGeneDict:
            representativePaftolTarget = None
            maxMappingScoreSum = None
            for organismName in result.paftolTargetSet.paftolGeneDict[geneName].paftolTargetDict:
                mappingScoreSum = result.paftolTargetSet.paftolGeneDict[geneName].paftolTargetDict[organismName].mappingScoreSum()
                if representativePaftolTarget is None or (mappingScoreSum is not None and mappingScoreSum > maxMappingScoreSum):
                    representativePaftolTarget = result.paftolTargetSet.paftolGeneDict[geneName].paftolTargetDict[organismName]
                    maxMappingScoreSum = mappingScoreSum
            result.representativePaftolTargetDict[geneName] = representativePaftolTarget
            if representativePaftolTarget is None:
                logger.debug('represenative for %s: none', geneName)
            else:
                logger.debug('representative for %s: %s', representativePaftolTarget.paftolGene.name, representativePaftolTarget.organism.name)

    def writeRepresentativeGenes(self, result):
        for geneName in result.representativePaftolTargetDict:
            paftolTarget = result.representativePaftolTargetDict[geneName]
            paftolTarget.writeFasta(self.makeGeneRepresentativeFname(geneName, True))

    def readMappedReadsSingle(self, result):
        readNameMappedReadDict = result.paftolTargetSet.makeReadNameMappedReadDict()
        with open(result.forwardFastq, 'r') as forwardFile:
            forwardParser = Bio.SeqIO.parse(forwardFile, 'fastq')
            for forwardRead in fowardParser:
                readName = forwardRead.id
                if readName in readNameMappedReadDict:
                    for mappedRead in readNameMappedReadDict[readName]:
                        if mappedRead.fowardRead is not None:
                            raise StandardError, 'duplicate forward read for %s' % readName
                        mappedRead.forwardRead = forwardRead

    def readMappedReadsPaired(self, result):
        readNameMappedReadDict = result.paftolTargetSet.makeReadNameMappedReadDict()
        with open(result.forwardFastq, 'r') as forwardFile:
            forwardParser = Bio.SeqIO.parse(forwardFile, 'fastq')
            with open(result.reverseFastq, 'r') as reverseFile:
                reverseParser = Bio.SeqIO.parse(reverseFile, 'fastq')
                for forwardRead in forwardParser:
                    reverseRead = reverseParser.next()
                    forwardRead.id = MappedRead.readBasename(forwardRead.id)
                    reverseRead.id = MappedRead.readBasename(reverseRead.id)
                    if reverseRead.id != forwardRead.id:
                        raise StandardError('paired read files %s / %s out of sync at read %s / %s' % (result.forwardFastq, result.reverseFastq, forwardRead.id, reverseRead.id))
                    readName = forwardRead.id
                    if readName in readNameMappedReadDict:
                        for mappedRead in readNameMappedReadDict[readName]:
                            if mappedRead.forwardRead is not None:
                                raise StandardError, 'duplicate forward read for %s' % readName
                            mappedRead.forwardRead = forwardRead
                            mappedRead.reverseRead = reverseRead
                # FIXME: check for dangling stuff in reverse: reverse.next() should trigger exception (StopIteration?)

    def writeMappedReadsFasta(self, result, maxNumReadsPerGene):
        for paftolGene in result.paftolTargetSet.paftolGeneDict.values():
            with open(self.makeGeneReadFname(paftolGene.name, True), 'w') as fastaFile:
                paftolGene.writeMappedReadsFasta(fastaFile, True, result.reverseFastq is not None, maxNumReadsPerGene)

    def distributeSingle(self, result, maxNumReadsPerGene):
        self.readMappedReadsSingle(result)
        self.writeMappedReadsFasta(result, maxNumReadsPerGene)

    def distributePaired(self, result, maxNumReadsPerGene):
        self.readMappedReadsPaired(result)
        self.writeMappedReadsFasta(result, maxNumReadsPerGene)

    def distributeSingleOld(self, result):
        fForward = open(result.forwardFastq, 'r')
        fqiForward = Bio.SeqIO.QualityIO.FastqGeneralIterator(fForward)
        readNameGeneDict = result.paftolTargetSet.makeReadNameGeneDict()
        for fwdReadTitle, fwdReadSeq, fwdReadQual in fqiForward:
            readName = fwdReadTitle.split()[0]
            if readName in readNameGeneDict:
                for paftolGene in readNameGeneDict[readName]:
                    with open(self.makeGeneReadFname(paftolGene.name, True), 'a') as f:
                        f.write('>%s\n%s\n' % (fwdReadTitle, fwdReadSeq))
        fForward.close()

    def distributePairedOld(self, result):
        # FIXME: consider try...finally to ensure files are closed
        fForward = open(result.forwardFastq, 'r')
        fqiForward = Bio.SeqIO.QualityIO.FastqGeneralIterator(fForward)
        fReverse = open(result.reverseFastq, 'r')
        fqiReverse = Bio.SeqIO.QualityIO.FastqGeneralIterator(fReverse)
        readNameGeneDict = result.paftolTargetSet.makeReadNameGeneDict()
        for fwdReadTitle, fwdReadSeq, fwdReadQual in fqiForward:
            readName = fwdReadTitle.split()[0]
            # FIXME: premature end of reverse fastq will trigger
            # StopIteration and premature end of forward will leave
            # rest of reverse ignored
            revReadTitle, revReadSeq, revReadQual = fqiReverse.next()
            if readName != revReadTitle.split()[0]:
                raise StandardError('paired read files %s / %s out of sync at read %s / %s' % (result.forwardFastq, result.reverseFastq, fwdReadTitle, revReadTitle))
            if readName in readNameGeneDict:
                for paftolGene in readNameGeneDict[readName]:
                    with open(self.makeGeneReadFname(paftolGene.name, True), 'a') as f:
                        f.write('>%s\n%s\n' % (fwdReadTitle, fwdReadSeq))
                        f.write('>%s\n%s\n' % (revReadTitle, revReadSeq))
        # FIXME: check for dangling stuff in reverse: should trigger
        # an exception:
        # revReadTitle, revReadSeq, revReadQual = fqiReverse.next()
        fForward.close()
        fReverse.close()

    def distribute(self, result, maxNumReadsPerGene):
        if result.isPaired():
            self.distributePaired(result, maxNumReadsPerGene)
        else:
            self.distributeSingle(result, maxNumReadsPerGene)

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
    def filterByOverlap(self, exonerateResultList, strictOverlapFiltering):

        def isPreferred(exonerateResult, other):
            if exonerateResult.queryAlignmentLength != other.queryAlignmentLength:
                return exonerateResult.queryAlignmentLength > other.queryAlignmentLength
            if exonerateResult.rawScore != other.rawScore:
                return exonerateResult.rawScore > other.rawScore
            if exonerateResult.queryAlignmentLength != other.queryAlignmentLength:
                return exonerateResult.queryAlignmentLength > other.queryAlignmentLength
            # FIXME: arbitrary tie breaking
            if exonerateResult.queryAlignmentStart != other.queryAlignmentStart:
                return exonerateResult.queryAlignmentStart < other.queryAlignmentStart
            if exonerateResult.queryAlignmentEnd != other.queryAlignmentEnd:
                return exonerateResult.queryAlignmentEnd < other.queryAlignmentEnd
            if exonerateResult.targetId != other.targetId:
                return exonerateResult.targetId < other.targetId
            raise StanddardError('cannot break tie of overlapping contigs: exonerateResult = %s, other = %s' % (str(exonerateResult), str(other)))

        logger.warning('scanning for overlaps but not resolving them, pending development of concept')
        nonOverlappingExonerateResultList = []
        for exonerateResult in exonerateResultList:
            isOverlapping = False
            for other in exonerateResultList:
                if exonerateResult is not other:
                    if exonerateResult.overlapsQueryAlignmentRange(other):
                        if strictOverlapFiltering:
                            isOverlapping = isOverlapping or (not isPreferred(exonerateResult, other))
                        else:
                            logger.warning('overlap found, but not resolved: %s, %s', str(exonerateResult), str(other))
            if not isOverlapping:
                nonOverlappingExonerateResultList.append(exonerateResult)
        return nonOverlappingExonerateResultList

    def filterExonerateResultList(self, geneName, exonerateResultList, strictOverlapFiltering):
        logger.debug('gene %s: %d exonerate results', geneName, len(exonerateResultList))
        exonerateResultList = self.filterByPercentIdentity(exonerateResultList)
        logger.debug('gene %s: %d sufficiently close exonerate results', geneName, len(exonerateResultList))
        exonerateResultList = self.filterByContainment(exonerateResultList)
        logger.debug('gene %s: %d non-contained exonerate results', geneName, len(exonerateResultList))
        # logger.debug('gene %s: non-contained contig list: %s', geneName, ', '.join([e.targetId for e in exonerateResultList]))
        exonerateResultList = self.filterByOverlap(exonerateResultList, strictOverlapFiltering)
        # logger.debug('gene %s: %d non-overlapping exonerate results', geneName, len(exonerateResultList))
        logger.debug('gene %s: %d non-overlapping contig list [strict=%s]: %s', geneName, len(exonerateResultList), str(strictOverlapFiltering), ', '.join(['%s; tcdsLen=%d' % (e.targetId, len(e.targetCdsSeq.seq)) for e in exonerateResultList]))
        return exonerateResultList

    def makeTgz(self):
        if self.workdirTgz is not None:
            if self.tmpDirname is None:
                raise StandardError('illegal state: no temporary directory generated')
            tmpTgz = os.path.join(self.tmpDirname, '%s.tgz' % self.workDirname)
            tgzArgv = ['tar', '-zcf', tmpTgz, self.workDirname]
            tgzProcess = subprocess.Popen(tgzArgv, cwd=self.tmpDirname)
            tgzReturncode = tgzProcess.wait()
            if tgzReturncode != 0:
                raise StandardError('process "%s" returned %d' % (' '.join(tgzArgv), tgzReturncode))
            # FIXME: clumsy to first create tgz in temp dir and then
            # moving it to final destination, compute absolute path to
            # final destination and use that directly?
            shutil.move(os.path.join(self.tmpDirname, tmpTgz), self.workdirTgz)

            
class TargetMapper(object):

    def __init__(self):
        self.workdir = None
        self.targetsFname = 'targets.fasta'
        
    def makePath(self, baseFname):
        return os.path.join(self.workdir, baseFname)
        
    def makeTargetsPath(self):
        return self.makePath(self.targetsFname)
    
    def writeTargetsFile(self, paftolTargetSet):
        paftolTargetSet.writeFasta(self.makeTargetsPath())
        
    def makeWorkdir(self):
        os.mkdir(self.workdir)
        
    def isSetup(self):
        return self.workdir is not None
        
    def setup(self, workdir):
        self.workdir = workdir
        self.makeWorkdir()
        logger.debug('set up workdir: %s', self.workdir)

    def cleanup(self):
        shutil.rmtree(self.workdir)
        self.workdir = None

    def mapReads(self, paftolTargetSet, forwardReadsFname, reverseReadsFname):
        raise StandardError, 'abstract method called'


class TargetMapperTblastn(TargetMapper):
    
    def __init__(self, tblastnRunner=None):
        super(TargetMapperTblastn, self).__init__()
        if tblastnRunner is None:
            tblastnRunner = paftol.tools.TblastnRunner()
        self.tblastnRunner = tblastnRunner
        self.forwardFasta = 'fwd.fasta'
        self.reverseFasta = 'rev.fasta'

    def mapReads(self, paftolTargetSet, forwardReadsFname, reverseReadsFname):
        logger.debug('mapping gene sequences to reads')
        if not self.isSetup():
            raise StandardError, 'illegal state: TargetMapperTblastn instance not set up'
        self.writeTargetsFile(paftolTargetSet)
        referenceFname = self.makeTargetsPath()
        targetProteinList = [paftol.tools.translateSeqRecord(geneSr) for geneSr in paftolTargetSet.getSeqRecordList()]
        # FIXME: check these parameters, consider numAlignments?
        self.tblastnRunner.maxTargetSeqs = 10000000
        self.tblastnRunner.maxHsps = 1
        # jtk: analyser logic needs developing to mange numOfftargetsReads generically
        # result.paftolTargetSet.numOfftargetReads = None
        forwardFastaPath = self.makePath(self.forwardFasta)
        paftol.tools.fastqToFasta(forwardReadsFname, forwardFastaPath)
        self.tblastnRunner.indexDatabase(forwardFastaPath)
        self.tblastnRunner.processTblastn(paftolTargetSet, forwardFastaPath, targetProteinList)
        if reverseReadsFname is not None:
            reverseFastaPath = self.makePath(self.reverseFasta)
            paftol.tools.fastqToFasta(reverseReadsFname, reverseFastaPath)
            self.tblastnRunner.indexDatabase(reverseFastaPath)
            self.tblastnRunner.processTblastn(paftolTargetSet, reverseFastaPath, targetProteinList)


class TargetMapperBwa(TargetMapper):
    
    def __init__(self, bwaRunner=None):
        super(TargetMapperBwa, self).__init__()
        if bwaRunner is None:
            bwaRunner = paftol.tools.BwaRunner()
        self.bwaRunner = bwaRunner
 
    def mapReads(self, paftolTargetSet, forwardReadsFname, reverseReadsFname):
        """Map reads to genes.
"""
        logger.debug('mapping reads to gene sequences')
        if self.workdir is None:
            raise StandardError, 'illegal state: no workdir'
        referenceFname = self.makeTargetsPath()
        self.bwaRunner.indexReference(referenceFname)
        forwardReadsFname = os.path.join(os.getcwd(), result.forwardFastq)
        if result.reverseFastq is None:
            reverseReadsFname = None
        else:
            reverseReadsFname = os.path.join(os.getcwd(), result.reverseFastq)
        paftolTargetSet.numOfftargetReads = 0
        self.bwaRunner.processBwa(paftolTargetSet, referenceFname, forwardReadsFname, reverseReadsFname)

        
class TargetAssembler(object):
    
    def __init__(self):
        self.workdir = None
        self.geneReadFnamePattern = 'gene-%s.fasta'

    def makeWorkdir(self):
        os.mkdir(self.workdir)
        
    def isSetup(self):
        return self.workdir is not None
        
    def setup(self, workdir):
        self.workdir = workdir
        self.makeWorkdir()

    def cleanup(self):
        shutil.rmtree(self.workdir)
        self.workdir = None

    def makeWorkdirPath(self, filename):
        return os.path.join(self.workdir, filename)

    def makeGeneReadFname(self, geneName):
        geneReadFname = self.geneReadFnamePattern % geneName
        return self.makeWorkdirPath(geneReadFname)


class TargetAssemblerOverlapSerial(TargetAssembler):

    def __init__(self, semiglobalAlignmentRunner=None):
        super(TargetAssemblerOverlapSerial, self).__init__()
        self.overlapCsvFnamePattern = 'overlap-%s.csv'
        self.positionedReadFnamePattern = 'posread-%s.fasta'
        self.semiglobalAlignmentRunner = semiglobalAlignmentRunner
        self.windowSizeReference = None
        self.relIdentityThresholdReference = None
        self.windowSizeReadOverlap = None
        self.relIdentityThresholdReadOverlap = None
        
    def makeWorkdirPath(self, filename):
        return os.path.join(self.workdir, filename)

    def checkControlParameters(self):
        controlParameterList = ['windowSizeReference', 'relIdentityThresholdReference', 'windowSizeReadOverlap', 'relIdentityThresholdReadOverlap', 'semiglobalAlignmentRunner']
        missingControlParameterList = []
        for controlParameter in controlParameterList:
            if self.__dict__[controlParameter] is None:
                missingControlParameterList.append(controlParameter)
        if len(missingControlParameterList) > 0:
            raise StandardError, 'missing control parameters: %s' % ', '.join(missingControlParameterList)

    def assembleGene(self, result, geneName):
        # logger.debug('tracking: starting with gene %s' % geneName)
        if not self.isSetup():
            raise StandardError, 'cannot assemble gene %s: TargetAssemblerOverlapSerial instance not set up' % geneName
        self.checkControlParameters()
        overlapCsvFname = self.makeWorkdirPath('overlap-%s.csv' % geneName)
        positionedReadDirname = self.makeWorkdirPath('posread-%s' % geneName)
        positionedReadFname = self.makeWorkdirPath('posread-%s.fasta' % geneName)
        os.mkdir(positionedReadDirname)
        readSrFwdList = copy.deepcopy(result.paftolTargetSet.paftolGeneDict[geneName].makeMappedReadsUniqueList(includeForward=True, includeReverse=False))
        readSrRevList = copy.deepcopy(result.paftolTargetSet.paftolGeneDict[geneName].makeMappedReadsUniqueList(includeForward=False, includeReverse=True))
        readSrList = []
        for readSr in readSrFwdList:
            readSr.id = '%s-fwd' % readSr.id
            readSrList.append(readSr)
        for readSr in readSrRevList:
            readSr.id = '%s-rev' % readSr.id
            readSrList.append(readSr)
        readSrList.extend(paftol.tools.reverseComplementSeqRecordList(readSrList))
        repGene = result.representativePaftolTargetDict[geneName].seqRecord
        # repGeneProtein = self.translateGene(repGene)
        geneReadFname = self.makeGeneReadFname(geneName)
        # for sr in readSrList:
        #     sys.stderr.write('%s\n' % sr.id)
        readSrDict = Bio.SeqIO.to_dict(readSrList)
        alignmentList = paftol.tools.semiglobalOneVsAll(repGene, readSrList)
        numReads = len(readSrList)
        if len(alignmentList) != numReads:
            raise StandardError, 'readSrList / alignment mismatch'
        positionedReadList = []
        for i in xrange(numReads):
            # sys.stderr.write('%s / %s: %f\n' % (alignmentList[i][0].id, alignmentList[i][1].id, findMaxRelativeIdentity(alignmentList[i], self.windowSize)))
            maxRelativeIdentity = paftol.tools.findMaxRelativeIdentity(alignmentList[i], self.windowSizeReference)
            if maxRelativeIdentity >= self.relIdentityThresholdReference:
                position = paftol.tools.findReadPosition(alignmentList[i])
                # coreAlignment = findOverlapAlignment(alignmentList[i])
                # coreLength = coreAlignment.get_alignment_length()
                # coreMatch = findRelativeIdentity(coreAlignment)
                coreLength = None
                coreMatch = None
                positionedReadList.append(paftol.tools.PositionedRead(readSrList[i], position, maxRelativeIdentity, coreLength, coreMatch))
            else:
                pass
                # sys.stderr.write('skipping read %s with maxRelativeIdentity %f\n' % (readSrList[i].id, maxRelativeIdentity))
        positionedReadList.sort()
        # logger.debug('tracking: positioned reads')
        positionedSrList = [positionedRead.readSr for positionedRead in positionedReadList]
        if positionedReadDirname is not None:
            for i in xrange(len(positionedSrList)):
                Bio.SeqIO.write([positionedSrList[i]], '%s/p%03d.fasta' % (positionedReadDirname, i), 'fasta')
        if positionedReadFname is not None:
            Bio.SeqIO.write(positionedSrList, positionedReadFname, 'fasta')
        if overlapCsvFname is not None:
            overlapDataFrame = paftol.tools.DataFrame(['read0', 'read1', 'read1pos', 'maxRelId', 'coreLength', 'coreMatch', 'overlapLength', 'overlapMatch'])
        else:
            overlapDataFrame = None
        contigList = []
        currentContig = paftol.tools.Contig(self.windowSizeReadOverlap, self.relIdentityThresholdReadOverlap, self.semiglobalAlignmentRunner)
        for i in xrange(len(positionedReadList)):
            if overlapDataFrame is not None:
                if i == 0:
                    overlapRow = {'read0': None, 'read1': positionedReadList[i].readSr.id, 'read1pos': positionedReadList[i].position, 'maxRelId': positionedReadList[i].maxRelativeIdentity, 'coreLength': positionedReadList[i].coreLength, 'coreMatch': positionedReadList[i].coreMatch, 'overlapLength': None, 'overlapMatch': None}
                else:
                    alignmentList = self.semiglobalAlignmentRunner.align(positionedReadList[i - 1].readSr, [positionedReadList[i].readSr])
                    alignment = alignmentList[0]
                    overlapAlignment = paftol.tools.findOverlapAlignment(alignment)
                    if overlapAlignment.get_alignment_length() == 0:
                        overlapMatch = None
                    else:
                        overlapMatch = paftol.tools.findRelativeIdentity(overlapAlignment)
                    overlapRow = {'read0': positionedReadList[i - 1].readSr.id, 'read1': positionedReadList[i].readSr.id, 'read1pos': positionedReadList[i].position, 'maxRelId': positionedReadList[i].maxRelativeIdentity, 'coreLength': positionedReadList[i].coreLength, 'coreMatch': positionedReadList[i].coreMatch, 'overlapLength': overlapAlignment.get_alignment_length(), 'overlapMatch': overlapMatch}
                overlapDataFrame.addRow(overlapRow)
            if currentContig.addRead(positionedReadList[i].readSr):
                logger.debug('added read %s to current contig', positionedReadList[i].readSr.id)
            else:
                logger.debug('started new contig with read %s', positionedReadList[i].readSr.id)
                currentContig.removeTerminalGaps()
                contigList.append(currentContig)
                currentContig = paftol.tools.Contig(self.windowSizeReadOverlap, self.relIdentityThresholdReadOverlap, self.semiglobalAlignmentRunner)
                currentContig.addRead(positionedReadList[i].readSr)
        currentContig.removeTerminalGaps()
        contigList.append(currentContig)
        if overlapCsvFname is not None:
            with open(overlapCsvFname, 'w') as f:
                overlapDataFrame.writeCsv(f)
        consensusList = []
        contigNumber = 0
        for contig in contigList:
            consensus = contig.getConsensus()
            if consensus is not None:
                consensus.id = '%s--contig%05d' % (geneName, contigNumber)
                contigNumber = contigNumber + 1
                consensusList.append(consensus)
        return consensusList

    
class TargetRecoverer(HybseqAnalyser):

    def __init__(self, workdirTgz, workDirname, trimmomaticRunner=None, targetMapper=None, targetAssembler=None):
        super(TargetRecoverer, self).__init__(workdirTgz, workDirname)
        self.trimmomaticRunner = trimmomaticRunner
        self.targetMapper = targetMapper
        self.targetAssembler = targetAssembler
        self.targetMapperWorkdir = 'targetmapper'
        self.targetAssemblerWorkdir = 'targetassembler'
        self.trimmedPairedFwd = 'trimmed_paired_fwd.fastq'
        self.trimmedPairedRev = 'trimmed_paired_rev.fastq'
        self.trimmedUnpairedFwd = 'trimmed_unpaired_fwd.fastq'
        self.trimmedUnpairedRev = 'trimmed_unpaired_rev.fastq'
        self.trimlogFname = 'trimlog.txt'

    def setup(self, result):
        logger.debug('setting up')
        self.setupTmpdir()
        self.targetMapper.setup(self.makeWorkdirPath(self.targetMapperWorkdir))
        self.targetAssembler.setup(self.makeWorkdirPath(self.targetAssemblerWorkdir))
        # FIXME: is writing the targets fasta file really part of setup?
	# result.paftolTargetSet.writeFasta(self.makeTargetsFname(True))
        
    # ideas for hybrid / consensus sequence for (multiple) re-mapping
    # reference CDS:     atgtac------catacagaagagacgtga
    # reconstructed CDS:    cactcatttcat---gga
    # "consensus"        atgCACTCAATTCAT   GGAgagacgtga
    # principe: Where reconstructed symbol is available, use that in preference.
    #   * gap in reference: use symbols from reconstructed (must be non-gap if pairwise alignment)
    #   * gap in reconstructed: skip symbols from reference
    #   * ends / portions with no alignment to reconstructed: fill in from reference
    # Problem: avoid non-homologous alignment portions (e.g. around borders of reconstructed)?

    def recoverContigs(self, result, geneName):
        logger.debug('reconstructing CDS for gene %s', geneName)
        if result.representativePaftolTargetDict is None:
            raise StandardError('illegal state: no represesentative genes')
        if result.representativePaftolTargetDict[geneName] is None:
            raise StandardError('no representative for gene %s' % geneName)
        os.mkdir(self.makeGeneDirPath(geneName))
        contigList = self.targetAssembler.assembleGene(result, geneName)
        if contigList is None:
            logger.warning('gene %s: no contigs', geneName)
            return None
        if len(contigList) == 0:
            logger.warning('gene %s: empty contig list', geneName)
            return None
        logger.debug('gene %s: %d contigs', geneName, len(contigList))
        for contig in contigList:
            contig.description = 'representativeGene=%s, targetsFasta=%s, %s' % (result.representativePaftolTargetDict[geneName].getName(), result.paftolTargetSet.fastaHandleStr, contig.description)
        return contigList
    
    def reconstructCds(self, result, geneName, strictOverlapFiltering):
        if geneName not in result.contigDict:
            raise StandardError, 'no contig recovery result for gene %s' % geneName
        contigList = result.contigDict[geneName]
        if contigList is None:
            logger.warning('gene %s: no cds reconstruction possible because no contigs were recovered' % geneName)
            return None
        geneProtein = self.translateGene(result.representativePaftolTargetDict[geneName].seqRecord)
        exonerateRunner = paftol.tools.ExonerateRunner()
        Bio.SeqIO.write([geneProtein], self.makeWorkdirPath('%s-protein.fasta' % geneName), 'fasta')
        aminoAcidSet = set(Bio.Alphabet.IUPAC.protein.letters.lower())
        # allow stop translation
        aminoAcidSet.add('*')
        setDiff = set(str(geneProtein.seq).lower()) - aminoAcidSet
        if len(setDiff) > 0:
            logger.warning('gene %s: invalid amino acids %s' % (geneName, ', '.join(setDiff)))
            return None
        contigFname = self.makeGeneContigsFname(geneName)
        Bio.SeqIO.write(contigList, contigFname, 'fasta')
        exonerateResultList = exonerateRunner.parse(geneProtein, contigFname, 'protein2genome', bestn=len(contigList))
        logger.debug('gene %s: %d contigs, %d exonerate results', geneName, len(contigList), len(exonerateResultList))
        if len(exonerateResultList) == 0:
            logger.warning('gene %s: no exonerate results from %d contigs', geneName, len(contigList))
        exonerateResultList.sort(paftol.cmpExonerateResultByQueryAlignmentStart)
        # reverse complementing extraneous as that is done by exonerate itself
        # for exonerateResult in exonerateResultList:
        #     logger.debug('gene %s, contig %s: targetStrand = %s', geneName, exonerateResult.targetId, exonerateResult.targetStrand)
        #     logger.debug('gene %s, contig %s, raw: %d -> %d, %s', geneName, exonerateResult.targetId, exonerateResult.targetCdsStart, exonerateResult.targetCdsEnd, str(exonerateResult.targetAlignmentSeq.seq))
        #     if exonerateResult.targetStrand == '-':
        #         exonerateResult.reverseComplementTarget()
        #     logger.debug('gene %s, contig %s, can: %d -> %d, %s', geneName, exonerateResult.targetId, exonerateResult.targetCdsStart, exonerateResult.targetCdsEnd, str(exonerateResult.targetAlignmentSeq.seq))
        # logger.warning('provisional filtering and supercontig construction, handling of overlapping contigs not finalised')
        filteredExonerateResultList = self.filterExonerateResultList(geneName, exonerateResultList, strictOverlapFiltering)
        logger.debug('gene %s: %d exonerate results after filtering', geneName, len(filteredExonerateResultList))
        Bio.SeqIO.write([e.targetCdsSeq for e in filteredExonerateResultList], os.path.join(self.makeGeneDirPath(geneName), '%s-fecds.fasta' % geneName), 'fasta')
        if len(filteredExonerateResultList) == 0:
            logger.warning('gene %s: no exonerate results left after filtering', geneName)
            return None
        supercontig = Bio.SeqRecord.SeqRecord(Bio.Seq.Seq(''.join([str(e.targetCdsSeq.seq) for e in filteredExonerateResultList])), id='%s_supercontig' % geneName)
        logger.debug('gene %s: supercontig length %d', geneName, len(supercontig))
        if len(supercontig) == 0:
            logger.warning('gene %s: empty supercontig', geneName)
            return None
        supercontigFname = os.path.join(self.makeGeneDirPath(geneName), '%s-supercontig.fasta' % geneName)
        Bio.SeqIO.write([supercontig], supercontigFname, 'fasta')
        Bio.SeqIO.write([geneProtein], os.path.join(self.makeGeneDirPath(geneName), '%s-supercontigref.fasta' % geneName), 'fasta')
        # FIXME: use exonerate to align "supercontig" to reference and
        # retrieve coding sequence of exonerate result with highest
        # score. In case of tied highest score, select result with
        # shortest CDS, as this is indicative of highest
        # "concentration" of matches and fewest gaps.
        supercontigErList = exonerateRunner.parse(geneProtein, supercontigFname, 'protein2genome', bestn=1)
        logger.debug('gene %s: %d supercontig exonerate results', geneName, len(supercontigErList))
        splicedSupercontigEr = None
        if len(supercontigErList) == 0:
            logger.warning('gene %s: no exonerate results from supercontig', geneName)
            return None
        if len(supercontigErList) > 1:
            splicedSupercontigEr = supercontigErList[0]
            minLength = len(splicedSupercontigEr.targetCdsSeq)
            for supercontigEr in supercontigErList:
                if len(supercontigEr.targetCdsSeq) < minLength:
                    splicedSupercontigEr = supercontigEr
                    minLength = len(splicedSupercontigEr.targetCdsSeq)
            contigStats = ', '.join(['raw=%d, cdsLen=%d' % (e.rawScore, len(e.targetCdsSeq)) for e in supercontigErList])
            logger.warning('gene %s: received %d supercontig exonerate results despite bestn=1 (%s), selected raw=%d, cdsLen=%d', geneName, len(supercontigErList), contigStats, splicedSupercontigEr.rawScore, len(splicedSupercontigEr.targetCdsSeq))
        else:
            splicedSupercontigEr = supercontigErList[0]
        # not filtering for percent identity to gene again, as that is already done
        if result.reverseFastq is not None:
            readsSpec = '%s, %s' % (result.forwardFastq, result.reverseFastq)
        else:
            readsSpec = result.forwardFastq
        splicedSupercontig = Bio.SeqRecord.SeqRecord(Bio.Seq.Seq(str(splicedSupercontigEr.targetCdsSeq.seq)), id=geneName, description='reconstructed CDS computed by paftol.overlapAnalyser, targets: %s, reads: %s' % (result.paftolTargetSet.fastaHandleStr, readsSpec))
        logger.debug('gene %s: splicedSupercontig length %d', geneName, len(splicedSupercontig))
        splicedSupercontigFname = os.path.join(self.makeGeneDirPath(geneName), '%s-splicedsupercontig.fasta' % geneName)
        Bio.SeqIO.write([splicedSupercontig], splicedSupercontigFname, 'fasta')
        return splicedSupercontig
    
    def analyse(self, targetsSourcePath, forwardFastq, reverseFastq, allowInvalidBases, strictOverlapFiltering, maxNumReadsPerGene):
        raise StandardError, 'obsolete -- use recoverTargets'

    def recoverTargets(self, targetsSourcePath, forwardFastq, reverseFastq, allowInvalidBases, strictOverlapFiltering, maxNumReadsPerGene):
        logger.debug('starting')
	paftolTargetSet = PaftolTargetSet()
	paftolTargetSet.readFasta(targetsSourcePath)
        # FIXME: put allowInvalidBases in result for subsequent reference?
	paftolTargetSet.sanityCheck(allowInvalidBases)
        result = HybpiperResult(paftolTargetSet, forwardFastq, reverseFastq)
	try:
            self.setup(result)
            logger.debug('setup done')
            if self.trimmomaticRunner is None:
                logger.debug('running without trimming')
                trimmedForwardPairedFastqPath = forwardFastq
                trimmedReversePairedFastqPath = reverseFastq
            else:
                logger.debug('running with trimmomatic trimming')
                trimmedForwardPairedFastqPath = self.makeWorkdirPath(self.trimmedPairedFwd)
                trimmedReversePairedFastqPath = self.makeWorkdirPath(self.trimmedPairedRev)
                trimmedForwardUnpairedFastqPath = self.makeWorkdirPath(self.trimmedUnpairedFwd)
                trimmedReverseUnpairedFastqPath = self.makeWorkdirPath(self.trimmedUnpairedRev)
                trimlogPath = self.makeWorkdirPath(self.trimlogFname)
                self.trimmomaticRunner.runTrimmomaticPaired(forwardFastq, reverseFastq, trimmedForwardPairedFastqPath, trimmedReversePairedFastqPath, trimmedForwardUnpairedFastqPath, trimmedReverseUnpairedFastqPath, trimlogFname = trimlogPath)
                result.forwardFastqTrimmedPaired = trimmedForwardPairedFastqPath
                result.reverseFastqTrimmedPaired = trimmedReversePairedFastqPath
                result.forwardFastqTrimmedUnpaired = trimmedForwardUnpairedFastqPath
                result.reverseFastqTrimmedUnpaired = trimmedReverseUnpairedFastqPath
                result.generateFastqcStats()
            self.targetMapper.mapReads(paftolTargetSet, trimmedForwardPairedFastqPath, trimmedReversePairedFastqPath)
            logger.debug('mapping done')
            self.distribute(result, maxNumReadsPerGene)
            logger.debug('read distribution done')
            self.setRepresentativeGenes(result)
            self.writeRepresentativeGenes(result)
            logger.debug('representative genes selected')
            result.contigDict = {}
            result.reconstructedCdsDict = {}
            for geneName in result.paftolTargetSet.paftolGeneDict:
                result.contigDict[geneName] = self.recoverContigs(result, geneName)
                result.reconstructedCdsDict[geneName] = self.reconstructCds(result, geneName, strictOverlapFiltering)
	    logger.debug('CDS reconstruction done')
            logger.debug('finished')
            return result
        finally:
            self.makeTgz()
            logger.debug('tgz file made')
            self.targetMapper.cleanup()
            self.targetAssembler.cleanup()
            self.cleanup()
            logger.debug('cleanup done')


class MappedRead(object):
    """Represent a mapping of an NGS read to a PaftolTarget.
"""

    slashNumberReadIdRe = re.compile('([^/]+)/([12])')

    def __init__(self, paftolTarget):
        self.paftolTarget = paftolTarget
        self.forwardRead = None
        self.reverseRead = None

    def getReadName(self):
        raise StandardError, 'abstract method not overridden'

    def getMappingScore(self):
        raise StandardError, 'abstract method not overridden'

    @staticmethod
    def readBasename(rawReadName):
        m = MappedRead.slashNumberReadIdRe.match(rawReadName)
        if m is None:
            return rawReadName
        else:
            return m.group(1)


class SamMappedRead(MappedRead):

    def __init__(self, paftolTarget, samAlignment):
        super(SamMappedRead, self).__init__(paftolTarget)
	self.samAlignment = samAlignment

    def getReadName(self):
        return self.readBasename(self.samAlignment.qname)

    def getMappingScore(self):
        return self.samAlignment.mapq


class BlastMappedRead(MappedRead):
    """Represent a mapping of a read to a PaftolTarget based on BLAST.

The mapping score is based on the best HSP only, other HSPs are not
taken into account.
"""

    def __init__(self, paftolTarget, blastAlignment):
        super(BlastMappedRead, self).__init__(paftolTarget)
	self.blastAlignment = blastAlignment

    def getReadName(self):
        return self.readBasename(self.blastAlignment.hit_id)

    def getMappingScore(self):
        # FIXME: returning score of HSP #0, ignoring all subsequent HSPs
        return self.blastAlignment.hsps[0].score


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
@ivar mappedReadList: the list of reads mapped to this target
@type seqRecord: C{list} of C{MappedRead}
"""

    csvFieldNames = ['organism', 'gene', 'seqLength', 'numMappedReads']

    def __init__(self, organism, paftolGene, seqRecord):
        self.organism = organism
        self.paftolGene = paftolGene
        self.seqRecord = seqRecord
        self.mappedReadList = []
	# self.readAssociationList = []
        if paftolGene.name in organism.paftolTargetDict or organism.name in paftolGene.paftolTargetDict:
            raise StandardError('duplicate organism/gene: organism = %s, gene = %s, seqId = %s' % (organism.name, paftolGene.name, seqRecord.id))
        organism.paftolTargetDict[paftolGene.name] = self
        paftolGene.paftolTargetDict[organism.name] = self

    def getName(self):
        return '%s-%s' % (self.organism.name, self.paftolGene.name)

    def addMappedRead(self, mappedRead):
        self.mappedReadList.append(mappedRead)

    def mappingScoreSum(self):
        if len(self.mappedReadList) == 0:
            return None
        return sum([mr.getMappingScore() for mr in self.mappedReadList])

    def getReadNameSet(self):
        # FIXME: may have to trim away "/1", "/2"?
        return set([mr.getReadName() for mr in self.mappedReadList])

    def numMappedReads(self):
        # FIXME: need to check for duplicates
        # return len(self.mappedReadList)
        return len(self.getReadNameSet())

    def csvRowDict(self):
        d = {}
        d['organism'] = self.organism.name
        d['gene'] = self.paftolGene.name
        if self.seqRecord is None:
            d['seqRecord'] = None
        else:
            d['seqLength'] = len(self.seqRecord)
        d['numMappedReads'] = self.numMappedReads()
        return d

    def writeFasta(self, fastaHandle):
        Bio.SeqIO.write([self.seqRecord], fastaHandle, 'fasta')

class Organism(object):

    """Represent an organism (in the GenBank / NCBI sense of the term).

@ivar name: this organism's name
@type name: C{str}
@ivar paftolTargetDict: dictionary of genes in this organism
@type paftolTargetDict: C{dict} of C{PaftolTarget} instances with PAFTOL gene names as keys
"""

    csvFieldNames = ['organism', 'numGenes', 'numMappedReads']

    def __init__(self, name):
        self.name = name
        self.paftolTargetDict = {}

    def numMappedReads(self):
        return sum([t.numMappedReads() for t in self.paftolTargetDict.values()])

    def csvRowDict(self):
        d = {}
        d['organism'] = self.name
        d['numGenes'] = len(self.paftolTargetDict)
        d['numMappedReads'] = self.numMappedReads()
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

    csvFieldNames = ['gene', 'numOrganisms', 'meanSeqLength', 'numMappedReads']

    def __init__(self, name):
        self.name = name
        self.paftolTargetDict = {}
        self.annotationDict = {}

    def getReadNameSet(self):
        s = set()
        for paftolTarget in self.paftolTargetDict.values():
            s = s | paftolTarget.getReadNameSet()
        return s

    def meanSequenceLength(self):
        if len(self.paftolTargetDict) == 0:
            return None
        else:
            return float(sum([len(t.seqRecord) for t in self.paftolTargetDict.values()])) / float(len(self.paftolTargetDict))

    def numMappedReads(self):
        return sum([t.numMappedReads() for t in self.paftolTargetDict.values()])

    def csvRowDict(self):
        d = {}
        d['gene'] = self.name
        d['numOrganisms'] = len(self.paftolTargetDict)
        d['meanSeqLength'] = self.meanSequenceLength()
        d['numMappedReads'] = self.numMappedReads()
        return d

    def makeMappedReadsUniqueList(self, includeForward=True, includeReverse=True):
        readNameSet = set()
        srList = []
        numReads = 0
        for paftolTarget in self.paftolTargetDict.values():
            for mappedRead in paftolTarget.mappedReadList:
                readName = mappedRead.getReadName()
                if readName not in readNameSet:
                    numReads = numReads + 1
                    readNameSet.add(readName)
                    if includeForward:
                        if mappedRead.forwardRead is None:
                            raise StandardError, 'mapped read %s: no forward read SeqRecord' % mappedRead.getReadName()
                        srList.append(mappedRead.forwardRead)
                    if includeReverse:
                        if mappedRead.reverseRead is None:
                            raise StandardError, 'mapped read %s: no reverse read SeqRecord' % mappedRead.getReadName()
                        srList.append(mappedRead.reverseRead)
        return srList

    def writeMappedReadsFasta(self, fastaHandle, writeForward=True, writeReverse=True, maxNumReads=None):
        readsList = self.makeMappedReadsUniqueList(writeForward, writeReverse)
        if maxNumReads is not None and len(readsList) > maxNumReads:
            selectedReadsList = tools.selectLongestReads(readsList, maxNumReads)
        else:
            selectedReadsList = readsList
        logger.debug('maxNumReads: %s, numReads: %d, selected: %d', str(maxNumReads), len(readsList), len(selectedReadsList))
        Bio.SeqIO.write(selectedReadsList, fastaHandle, 'fasta')


def extractOrganismAndGeneNames(s):
    # FIXME: should tighten this up to fail on dangling garbage (?)
    paftolTargetRe = re.compile('([^-]+)-([^-]+)')
    m = paftolTargetRe.match(s)
    if m is not None:
        organismName = m.group(1)
        geneName = m.group(2)
    else:
        organismName = 'unknown'
        geneName = s
    return organismName, geneName


class PaftolTargetSet(object):
    """Represent a set of PAFTOL targets.

This class supports mapping using C{bwa} and C{tblastn} by
implementing the processSamAlignment and processBlastAlignment
methods, respectively.
"""


    def __init__(self):
        self.paftolGeneDict = {}
        self.organismDict = {}
        self.numOfftargetReads = None
        self.fastaHandleStr = None

    # FIXME: static?
    def makeFastaId(self, organismName, geneName):
        return '%s-%s' % (organismName, geneName)

    # deprecated, use getSeqRecordSelection instead
    def getGeneSeqRecordList(self, geneNameList):
        srList = []
        for geneName in geneNameList:
            if geneName not in self.paftolGeneDict:
                raise StandardError, 'gene %s not found in this target set' % geneName
            for paftolTarget in self.paftolGeneDict[geneName].paftolTargetDict.values():
                srList.append(paftolTarget.seqRecord)
        return srList

    def getSeqRecordSelection(self, organismNameList=None, geneNameList=None):
        if organismNameList is None:
            organismList = self.organismDict.values()
        else:
            organismList = []
            for organismName in organismNameList:
                if organismName not in self.organismDict:
                    raise StandardError, 'organism %s not found in this target set' % organismName
                organismList.append(self.organismDict[organismName])
        if geneNameList is not None:
            for geneName in geneNameList:
                if geneName not in self.paftolGeneDict:
                    raise StandardError, 'gene %s not found in this target set' % geneName
        srList = []
        for organism in organismList:
            for geneName in organism.paftolTargetDict:
                if geneNameList is None or geneName in geneNameList:
                    srList.append(organism.paftolTargetDict[geneName].seqRecord)
        return srList

    def readFasta(self, fastaHandle):
        # FIXME: add provision to control tolerance for invalid bases -- checkTargets type functionality?
        self.paftolGeneDict = {}
        self.organismDict = {}
        self.fastaHandleStr = str(fastaHandle)
        for sr in Bio.SeqIO.parse(fastaHandle, 'fasta', alphabet=Bio.Alphabet.IUPAC.ambiguous_dna):
            organismName, geneName = extractOrganismAndGeneNames(sr.id)
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

    def checkOrganismAndGene(self, organismName, geneName):
        if organismName not in self.organismDict:
            raise StandardError('unknown organism: %s' % organismName)
        if geneName not in self.paftolGeneDict:
            raise StandardError('unknown gene: %s' % geneName)
        if geneName not in self.organismDict[organismName].paftolTargetDict:
            raise StandardError('no entry for gene %s in organism %s' % (geneName, organismName))

    def processSamAlignment(self, samAlignment):
        if samAlignment.isMapped():
            organismName, geneName = extractOrganismAndGeneNames(samAlignment.rname)
            self.checkOrganismAndGene(organismName, geneName)
            paftolTarget = self.organismDict[organismName].paftolTargetDict[geneName]
	    mappedRead = SamMappedRead(paftolTarget, samAlignment)
	    paftolTarget.addMappedRead(mappedRead)
        else:
            self.numOfftargetReads = self.numOfftargetReads + 1

    def processBlastAlignment(self, query, blastAlignment):
        organismName, geneName = extractOrganismAndGeneNames(query)
        self.checkOrganismAndGene(organismName, geneName)
        paftolTarget = self.organismDict[organismName].paftolTargetDict[geneName]
        mappedRead = BlastMappedRead(paftolTarget, blastAlignment)
        paftolTarget.addMappedRead(mappedRead)

    def makeReadNameGeneDict(self):
        readNameGeneDict = {}
        for paftolGene in self.paftolGeneDict.values():
            for readName in paftolGene.getReadNameSet():
                if readName not in readNameGeneDict:
                    readNameGeneDict[readName] = []
                readNameGeneDict[readName].append(paftolGene)
        return readNameGeneDict

    def makeReadNameMappedReadDict(self):
        # FIXME: not exactly exemplary for following law of Demeter -- inner parts of loop probably want to be PaftolTarget or MappedRead methods
        readNameMappedReadDict = {}
        for paftolGene in self.paftolGeneDict.values():
            for paftolTarget in paftolGene.paftolTargetDict.values():
                for mappedRead in paftolTarget.mappedReadList:
                    readName = mappedRead.getReadName()
                    if readName not in readNameMappedReadDict:
                        readNameMappedReadDict[readName] = []
                    readNameMappedReadDict[readName].append(mappedRead)
        return readNameMappedReadDict

    def targetStats(self):
        dataFrame = paftol.tools.DataFrame(PaftolTarget.csvFieldNames)
        for organism in self.organismDict.values():
            for paftolTarget in organism.paftolTargetDict.values():
                dataFrame.addRow(paftolTarget.csvRowDict())
        return dataFrame

    def geneStats(self):
        dataFrame = paftol.tools.DataFrame(PaftolGene.csvFieldNames)
        for paftolGene in self.paftolGeneDict.values():
            dataFrame.addRow(paftolGene.csvRowDict())
        return dataFrame

    def organismStats(self):
        dataFrame = paftol.tools.DataFrame(Organism.csvFieldNames)
        for organism in self.organismDict.values():
            dataFrame.addRow(organism.csvRowDict())
        return dataFrame

    def getMappedReadNameSet(self):
        mappedReadNameSet = set()
        for paftolGene in self.paftolGeneDict.values():
            mappedReadNameSet = mappedReadNameSet | paftolGene.getReadNameSet()
        return mappedReadNameSet

    def numMappedReads(self):
        return len(self.getMappedReadNameSet())
        # FIXME: no check for multiple counting -- should change to len(self.getMappedReadNameSet())
        # n = 0
        # for organism in self.organismDict.values():
        #     for paftolTarget in organism.paftolTargetDict.values():
        #         n = n + paftolTarget.numMappedReads()
        # return n

    # FIXME: apparently obsolete -- delete?
    def writeMappedReadsFasta(self, fastaFname):
        with open(fastaFname, 'w') as fastaFile:
            for organism in self.organismDict.values():
                for paftolTarget in organism.paftolTargetDict.values():
                    paftolTarget.writeMappedReadsFasta(fastaFile)

    # FIXME: untested after transplant from HybpiperAnalyser
    def sanityCheck(self, allowInvalidBases=False):
        for organism in self.organismDict.values():
            for paftolTarget in organism.paftolTargetDict.values():
		if not allowInvalidBases:
		    setDiff = set(str(paftolTarget.seqRecord.seq).lower()) - set('acgt')
		    if len(setDiff) != 0:
			raise StandardError('target %s: illegal base(s) %s' % (paftolTarget.seqRecord.id, ', '.join(setDiff)))


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

    def makeGenomicSeqRecord(self, seqId=None):
        s = self.geneFeature.extract(self.seqRecord.seq)
        if seqId is None:
            seqId = self.geneId
        return Bio.SeqRecord.SeqRecord(s, id=seqId, description='sequence extracted from gene feature')

    def makeMrnaSeqRecord(self, seqId=None):
        if len(self.mrnaFeatureList) == 0:
            return None
        if len(self.mrnaFeatureList) > 1:
            logger.warning('reference gene %s: %d mRNA features, using #0 to extract sequence', self.geneId, len(self.mrnaFeatureList))
        mrnaFeature = self.mrnaFeatureList[0]
        s = mrnaFeature.extract(self.seqRecord.seq)
        if seqId is None:
            seqId = self.geneId
        return Bio.SeqRecord.SeqRecord(s, id=seqId, description='sequence extracted from mRNA feature')

    def makeCdsSeqRecord(self, seqId=None):
        if len(self.cdsFeatureList) == 0:
            return None
        if len(self.cdsFeatureList) > 1:
            logger.warning('reference gene %s: %d CDS features, using #0 to extract sequence', self.geneId, len(self.cdsFeatureList))
        cdsFeature = self.cdsFeatureList[0]
        s = cdsFeature.extract(self.seqRecord.seq)
        if seqId is None:
            seqId = self.geneId
        return Bio.SeqRecord.SeqRecord(s, id=seqId, description='sequence extracted from CDS feature')


class ReferenceGenomeMappingProcessor(object):

    def __init__(self, referenceGenome):
        if referenceGenome.genomeLength is None:
            raise StandardError, 'reference genome length is None'
        self.referenceGenome = referenceGenome
        self.intergenicId = 'intergenic'
        self.unmappedId = 'unmapped'
        self.geneHitDict = {}
        self.intergenicLength = referenceGenome.genomeLength
        for gene in referenceGenome.geneList:
            geneLength = gene.getLength()
            self.geneHitDict[gene.geneId] = {'geneId': gene.geneId, 'geneLength': geneLength, 'numHits': 0}
            self.intergenicLength = self.intergenicLength - geneLength
        self.geneHitDict[self.intergenicId] = {'geneId': self.intergenicId, 'geneLength': self.intergenicLength, 'numHits': 0}
        self.geneHitDict[self.unmappedId] = {'geneId': self.unmappedId, 'geneLength': None, 'numHits': 0}
        self.rawmapTable = paftol.tools.DataFrame(['qname', 'rname', 'pos'])

    def getStatsTable(self):
        statsTable = paftol.tools.DataFrame(['geneId', 'geneLength', 'numHits'])
        for gene in self.referenceGenome.geneList:
            statsTable.addRow(self.geneHitDict[gene.geneId])
        statsTable.addRow(self.geneHitDict[self.intergenicId])
        statsTable.addRow(self.geneHitDict[self.unmappedId])
        return statsTable

    def processSamAlignment(self, samAlignment):
        if samAlignment.isMapped():
            self.rawmapTable.addRow({'qname': samAlignment.qname, 'rname': samAlignment.rname, 'pos': samAlignment.pos})
            geneId = self.referenceGenome.findGeneIdForSamAlignment(samAlignment)
            if geneId is None:
                geneId = self.intergenicId
        else:
            geneId = self.unmappedId
        self.geneHitDict[geneId]['numHits'] = self.geneHitDict[geneId]['numHits'] + 1


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

    def makeCdsListGeneric(self):

        def extractCdsId(qualifiers, qualifierId):
            qualifier = qualifiers[qualifierId]
            if len(qualifier) == 0:
                raise StandardError, 'qualifier %s has length 0' % qualifierId
            elif len(qualifier) > 1:
                logger.warning('qualifier %s has %d values, using [0] (%s)', qualifierId, len(qualifier), ', '.join(qualifier))
            return qualifier[0]

        cdsIdNumberDict = {}
        cdsList = []
        with open(self.genbankFname, 'r') as f:
            for seqRecord in Bio.SeqIO.parse(f, 'genbank'):
                for seqFeature in seqRecord.features:
                    if seqFeature.type == 'CDS':
                        if 'gene' in seqFeature.qualifiers:
                            cdsId = extractCdsId(seqFeature.qualifiers, 'gene')
                        elif 'locus_tag' in seqFeature.qualifiers:
                            cdsId = extractCdsId(seqFeature.qualifiers, 'locus_tag')
                        else:
                            cdsId = '%s%s' % (seqRecord.id, str(seqFeature.location))
                        # FIXME: ad hoc sanitising of cdsId, replacing spaces with underscores to ensure entire cdsId ends up in FASTA ID portion
                        # necessitated by makeblastdb, which may otherwise see entries with identical cdsId
                        cdsId = cdsId.upper().replace(' ', '_')
                        if cdsId not in cdsIdNumberDict:
                            cdsIdNumberDict[cdsId] = 0
                        cdsIdNumberDict[cdsId] = cdsIdNumberDict[cdsId] + 1
                        cdsIdNumbered = '%s_%d' % (cdsId, cdsIdNumberDict[cdsId])
                        cdsSeq = seqFeature.extract(seqRecord.seq)
                        cdsList.append(Bio.SeqRecord.SeqRecord(cdsSeq, id=cdsIdNumbered, description=str(seqFeature.location)))
        return cdsList

    def makeCdsList(self):
        return self.makeCdsListGeneric()

    def scanGenesAth(self):
        if self.genbankFname is None:
            raise StandardError('no GenBank file name, cannot scan genes (ath method)')
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
        blastnProcess = subprocess.Popen(blastnArgv, stdin=subprocess.PIPE, stdout=subprocess.PIPE)
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
        targetGeneTable = paftol.tools.DataFrame(['targetId', 'geneId', 'geneName', 'geneNote', 'mrnaProduct', 'cdsProduct'])
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
        cdsList = []
        for targetId in targetIdToGeneDict:
            for gene in targetIdToGeneDict[targetId]:
                cds = gene.makeCdsSeqRecord()
                if cds is None:
                    logger.debug('gene %s (target %s): no CDS', gene.geneId, targetId)
                else:
                    cds.description = '%s, matched to target %s' % (cds.description, targetId)
                    cdsList.append(cds)
        return targetGeneTable, cdsList

    def findGeneIdForSamAlignment(self, samAlignment):
        # FIXME: clumsy linear search
        for gene in self.geneList:
            if gene.containsSamAlignment(samAlignment):
                return gene.geneId

    def mapReadsStatsBwaMem(self, bwaRunner, forwardReadsFname, reverseReadsFname=None):
        referenceGenomeMappingProcessor = ReferenceGenomeMappingProcessor(self)
        bwaRunner.processBwa(referenceGenomeMappingProcessor, forwardReadsFname, reverseReadsFname)
        return referenceGenomeMappingProcessor.getStatsTable(), referenceGenomeMappingProcessor.rawmapTable


class PaftolTargetSeqRetriever(object):

    def __init__(self):
        self.blastAlignmentDict = None

    def processBlastAlignment(self, query, blastAlignment):
        organismName, geneName = extractOrganismAndGeneNames(query)
        if geneName in self.blastAlignmentDict:
            if blastAlignment.hsps[0].expect < self.blastAlignmentDict[geneName].hsps[0].expect:
                self.blastAlignmentDict[geneName] = blastAlignment
        else:
            self.blastAlignmentDict[geneName] = blastAlignment

    def retrievePaftolTargetList(self, genomeName, fastaFname, paftolTargetSet, blastnRunner=None):
        if blastnRunner is None:
            blastnRunner = tools.BlastnRunner()
        self.blastAlignmentDict = {}
        blastnRunner.processBlast(self, fastaFname, paftolTargetSet.getSeqRecordList())
        seqIdGeneDict = {}
        for geneName in self.blastAlignmentDict:
            seqId = self.blastAlignmentDict[geneName].hit_id
            if seqId in seqIdGeneDict:
                raise StandardError, 'multiple PAFTOL genes for %s: %s, %s' % (seqId, seqIdGeneDict[seqId], geneName)
            seqIdGeneDict[seqId] = geneName
        paftolTargetList = []
        for seqRecord in Bio.SeqIO.parse(fastaFname, 'fasta'):
            if seqRecord.id in seqIdGeneDict:
                seqId = seqRecord.id
                geneName = seqIdGeneDict[seqId]
                evalue = self.blastAlignmentDict[geneName].hsps[0].expect
                seqRecord.description = '%s, original ID: %s, evalue: %1.12g' % (seqRecord.description, seqId, evalue)
                # seqRecord.id = '%s-%s' % (genomeName, geneName, )
                seqRecord.id = '%s' % geneName
                paftolTargetList.append(seqRecord)
        return paftolTargetList


class HybpiperAnalyser(HybseqAnalyser):
    """L{HybseqAnalyser} subclass for HybPiper style analysis (SPAdes based).

@ivar spadesRunner: SPAdes runner, providing instance variables for configuring SPAdes
@type spadesRunner: C{paftol.tools.SpadesRunner}
"""

    def __init__(self, workdirTgz, workDirname, spadesRunner):
        """Initialiser.

@param workdirTgz: name of the tgz archive of the working directory
@type workdirTgz: C{String}
@param workDirname: name of the working directory
@type workDirname: C{String}
@param spadesRunner: SPAdes runner for this analyser to use, C{None} to create a default instance
@type spadesRunner: C{paftol.tools.SpadesRunner} instance
"""
        super(HybpiperAnalyser, self).__init__(workdirTgz, workDirname)
        if spadesRunner is None:
            self.spadesRunner = paftol.tools.SpadesRunner()
        else:
            self.spadesRunner = spadesRunner

    def assembleGeneSpades(self, result, geneName):
        geneReadFname = self.makeGeneReadFname(geneName)
        if not os.path.exists(self.makeWorkdirPath(geneReadFname)):
            logger.debug('gene fasta file %s does not exist (no reads?)', geneReadFname)
            return None
        if result.isPaired():
            # FIXME: tight implicit coupling with distributeSingle / distributePaired
            libraryType = paftol.tools.SpadesRunner.INTERLACED
        else:
            libraryType = paftol.tools.SpadesRunner.SINGLE
        spadesOutputDirname = self.makeGeneDirPath(geneName)
        spadesContigList = self.spadesRunner.assemble(geneReadFname, libraryType, spadesOutputDirname, self.makeWorkDirname())
        return spadesContigList

    def reconstructCds(self, result, geneName, strictOverlapFiltering):
        logger.debug('reconstructing CDS for gene %s', geneName)
        if result.representativePaftolTargetDict is None:
            raise StandardError('illegal state: no represesentative genes')
        if result.representativePaftolTargetDict[geneName] is None:
            raise StandardError('no representative for gene %s' % geneName)
        os.mkdir(self.makeGeneDirPath(geneName))
        contigList = self.assembleGeneSpades(result, geneName)
        if contigList is None:
            logger.warning('gene %s: no spades contigs', geneName)
            return None
        if len(contigList) == 0:
            logger.warning('gene %s: empty contig list', geneName)
            return None
        logger.debug('gene %s: %d spades contigs', geneName, len(contigList))
        geneProtein = self.translateGene(result.representativePaftolTargetDict[geneName].seqRecord)
        Bio.SeqIO.write([geneProtein], self.makeWorkdirPath('%s-protein.fasta' % geneName), 'fasta')
        aminoAcidSet = set(Bio.Alphabet.IUPAC.protein.letters.lower())
        # allow stop translation
        aminoAcidSet.add('*')
        setDiff = set(str(geneProtein.seq).lower()) - aminoAcidSet
        if len(setDiff) > 0:
            logger.warning('gene %s: invalid amino acids %s' % (geneName, ', '.join(setDiff)))
            return None
        contigFname = self.makeGeneContigsFname(geneName)
        Bio.SeqIO.write(contigList, contigFname, 'fasta')
        exonerateRunner = paftol.tools.ExonerateRunner()
        exonerateResultList = exonerateRunner.parse(geneProtein, contigFname, 'protein2genome', bestn=len(contigList))
        logger.debug('gene %s: %d contigs, %d exonerate results', geneName, len(contigList), len(exonerateResultList))
        if len(exonerateResultList) == 0:
            logger.warning('gene %s: no exonerate results from %d contigs', geneName, len(contigList))
        exonerateResultList.sort(cmpExonerateResultByQueryAlignmentStart)
        # reverse complementing extraneous as that is done by exonerate itself
        # for exonerateResult in exonerateResultList:
        #     logger.debug('gene %s, contig %s: targetStrand = %s', geneName, exonerateResult.targetId, exonerateResult.targetStrand)
        #     logger.debug('gene %s, contig %s, raw: %d -> %d, %s', geneName, exonerateResult.targetId, exonerateResult.targetCdsStart, exonerateResult.targetCdsEnd, str(exonerateResult.targetAlignmentSeq.seq))
        #     if exonerateResult.targetStrand == '-':
        #         exonerateResult.reverseComplementTarget()
        #     logger.debug('gene %s, contig %s, can: %d -> %d, %s', geneName, exonerateResult.targetId, exonerateResult.targetCdsStart, exonerateResult.targetCdsEnd, str(exonerateResult.targetAlignmentSeq.seq))
        # logger.warning('provisional filtering and supercontig construction, handling of overlapping contigs not finalised')
        filteredExonerateResultList = self.filterExonerateResultList(geneName, exonerateResultList, strictOverlapFiltering)
        logger.debug('gene %s: %d exonerate results after filtering', geneName, len(filteredExonerateResultList))
        Bio.SeqIO.write([e.targetCdsSeq for e in filteredExonerateResultList], os.path.join(self.makeGeneDirPath(geneName), '%s-fecds.fasta' % geneName), 'fasta')
        if len(filteredExonerateResultList) == 0:
            logger.warning('gene %s: no exonerate results left after filtering', geneName)
            return None
        supercontig = Bio.SeqRecord.SeqRecord(Bio.Seq.Seq(''.join([str(e.targetCdsSeq.seq) for e in filteredExonerateResultList])), id='%s_supercontig' % geneName)
        logger.debug('gene %s: supercontig length %d', geneName, len(supercontig))
        if len(supercontig) == 0:
            logger.warning('gene %s: empty supercontig', geneName)
            return None
        supercontigFname = os.path.join(self.makeGeneDirPath(geneName), '%s-supercontig.fasta' % geneName)
        Bio.SeqIO.write([supercontig], supercontigFname, 'fasta')
        Bio.SeqIO.write([geneProtein], os.path.join(self.makeGeneDirPath(geneName), '%s-supercontigref.fasta' % geneName), 'fasta')
        # FIXME: use exonerate to align "supercontig" to reference and
        # retrieve coding sequence of exonerate result with highest
        # score. In case of tied highest score, select result with
        # shortest CDS, as this is indicative of highest
        # "concentration" of matches and fewest gaps.
        supercontigErList = exonerateRunner.parse(geneProtein, supercontigFname, 'protein2genome', bestn=1)
        logger.debug('gene %s: %d supercontig exonerate results', geneName, len(supercontigErList))
        splicedSupercontigEr = None
        if len(supercontigErList) == 0:
            logger.warning('gene %s: no exonerate results from supercontig', geneName)
            return None
        if len(supercontigErList) > 1:
            splicedSupercontigEr = supercontigErList[0]
            minLength = len(splicedSupercontigEr.targetCdsSeq)
            for supercontigEr in supercontigErList:
                if len(supercontigEr.targetCdsSeq) < minLength:
                    splicedSupercontigEr = supercontigEr
                    minLength = len(splicedSupercontigEr.targetCdsSeq)
            contigStats = ', '.join(['raw=%d, cdsLen=%d' % (e.rawScore, len(e.targetCdsSeq)) for e in supercontigErList])
            logger.warning('gene %s: received %d supercontig exonerate results despite bestn=1 (%s), selected raw=%d, cdsLen=%d', geneName, len(supercontigErList), contigStats, splicedSupercontigEr.rawScore, len(splicedSupercontigEr.targetCdsSeq))
        else:
            splicedSupercontigEr = supercontigErList[0]
        # not filtering for percent identity to gene again, as that is already done
        if result.reverseFastq is not None:
            readsSpec = '%s, %s' % (result.forwardFastq, result.reverseFastq)
        else:
            readsSpec = result.forwardFastq
        splicedSupercontig = Bio.SeqRecord.SeqRecord(Bio.Seq.Seq(str(splicedSupercontigEr.targetCdsSeq.seq)), id=geneName, description='reconstructed CDS computed by paftol.HybpiperAnalyser, targets: %s, reads: %s' % (result.paftolTargetSet.fastaHandleStr, readsSpec))
        logger.debug('gene %s: splicedSupercontig length %d', geneName, len(splicedSupercontig))
        splicedSupercontigFname = os.path.join(self.makeGeneDirPath(geneName), '%s-splicedsupercontig.fasta' % geneName)
        Bio.SeqIO.write([splicedSupercontig], splicedSupercontigFname, 'fasta')
        return splicedSupercontig


class HybpiperBwaAnalyser(HybpiperAnalyser):

    """L{HybpiperAnalyser} subclass that implements an analysis process
close to the HybPiper pipeline.

Some parameters to SPAdes can be controlled via instance variables as
documented below. Defaults of these parameters correspond to the
defaults provided by SPAdes, respectively (at the time of developing
this).

@ivar bwaRunner: SPAdes runner, providing instance variables for configuring BWA
@type bwaRunner: C{paftol.tools.BwaRunner}

    """

    def __init__(self, workdirTgz=None, workDirname='pafpipertmp', bwaRunner=None, spadesRunner=None):
        super(HybpiperBwaAnalyser, self).__init__(workdirTgz, workDirname, spadesRunner)
        if bwaRunner is None:
            self.bwaRunner = paftol.tools.BwaRunner()
        else:
            self.bwaRunner = bwaRunner

    def setup(self, result):
        logger.debug('setting up')
        self.setupTmpdir()
        # FIXME: is writing the targets fasta file really part of setup?
	result.paftolTargetSet.writeFasta(self.makeTargetsFname(True))

    def mapReadsBwa(self, result):
        """Map reads to gene sequences (from multiple organisms possibly).
"""
        logger.debug('mapping reads to gene sequences')
        referenceFname = self.makeTargetsFname(True)
        self.bwaRunner.indexReference(referenceFname)
        forwardReadsFname = os.path.join(os.getcwd(), result.forwardFastq)
        if result.reverseFastq is None:
            reverseReadsFname = None
        else:
            reverseReadsFname = os.path.join(os.getcwd(), result.reverseFastq)
        result.paftolTargetSet.numOfftargetReads = 0
        self.bwaRunner.processBwa(result.paftolTargetSet, referenceFname, forwardReadsFname, reverseReadsFname)

    # ideas for hybrid / consensus sequence for (multiple) re-mapping
    # reference CDS:     atgtac------catacagaagagacgtga
    # reconstructed CDS:    cactcatttcat---gga
    # "consensus"        atgCACTCAATTCAT   GGAgagacgtga
    # principe: Where reconstructed symbol is available, use that in preference.
    #   * gap in reference: use symbols from reconstructed (must be non-gap if pairwise alignment)
    #   * gap in reconstructed: skip symbols from reference
    #   * ends / portions with no alignment to reconstructed: fill in from reference
    # Problem: avoid non-homologous alignment portions (e.g. around borders of reconstructed)?

    def analyse(self, targetsSourcePath, forwardFastq, reverseFastq, allowInvalidBases, strictOverlapFiltering, maxNumReadsPerGene):
        logger.debug('starting')
	paftolTargetSet = PaftolTargetSet()
	paftolTargetSet.readFasta(targetsSourcePath)
        # FIXME: put allowInvalidBases in result for subsequent reference?
	paftolTargetSet.sanityCheck(allowInvalidBases)
        result = HybpiperResult(paftolTargetSet, forwardFastq, reverseFastq)
	try:
            self.setup(result)
            logger.debug('setup done')
            self.mapReadsBwa(result)
            logger.debug('BWA mapping done')
            self.distribute(result, maxNumReadsPerGene)
            logger.debug('read distribution done')
            self.setRepresentativeGenes(result)
            self.writeRepresentativeGenes(result)
            logger.debug('representative genes selected')
            result.reconstructedCdsDict = {}
            for geneName in result.paftolTargetSet.paftolGeneDict:
                result.reconstructedCdsDict[geneName] = self.reconstructCds(result, geneName, strictOverlapFiltering)
	    logger.debug('CDS reconstruction done')
            logger.debug('finished')
            return result
        finally:
            self.makeTgz()
            logger.debug('tgz file made')
            self.cleanup()
            logger.debug('cleanup done')


class HybpiperTblastnAnalyser(HybpiperAnalyser):

    """L{HybseqAnalyser} subclass that implements an analysis process
close to the HybPiper pipeline.

The C{tblastnRunner} and C{spadesRunner} should be considered "owned" by
the analyser, i.e. they may be modified by the analyser (e.g. to set parameters
as required), so they should not be modified or otherwise be used by clients.

Some parameters to SPAdes can be controlled via instance variables as
documented below. Defaults of these parameters correspond to the
defaults provided by SPAdes, respectively (at the time of developing
this).

@ivar tblastnRunner: tblastn runner, providing instance variables for configuring tblastn
@type tblastnRunner: C{paftol.tools.TblastnRunner}
@ivar spadesRunner: SPAdes runner, providing instance variables for configuring SPAdes
@type spadesRunner: C{paftol.tools.SpadesRunner}

    """

    def __init__(self, workdirTgz=None, workDirname='pafpipertmp', tblastnRunner=None, spadesRunner=None):
        super(HybpiperTblastnAnalyser, self).__init__(workdirTgz, workDirname, spadesRunner)
        if tblastnRunner is None:
            self.tblastnRunner = paftol.tools.TblastnRunner()
        else:
            self.tblastnRunner = tblastnRunner

    def setup(self, result):
        logger.debug('setting up')
        self.setupTmpdir()
	result.paftolTargetSet.writeFasta(self.makeTargetsFname(True))
        forwardFastaPath = self.makeWorkdirPath(self.forwardFasta)
        paftol.tools.fastqToFasta(result.forwardFastq, forwardFastaPath)
        self.tblastnRunner.indexDatabase(forwardFastaPath)
        if result.reverseFastq is not None:
            reverseFastaPath = self.makeWorkdirPath(self.reverseFasta)
            paftol.tools.fastqToFasta(result.reverseFastq, reverseFastaPath)
            self.tblastnRunner.indexDatabase(reverseFastaPath)

    def mapReadsTblastn(self, result):
        """Map gene sequences to reads (from multiple organisms possibly).
"""
        logger.debug('mapping gene sequences to reads')
        referenceFname = self.makeTargetsFname(True) ## check this holds
        targetProteinList = [self.translateGene(geneSr) for geneSr in result.paftolTargetSet.getSeqRecordList()]
        # FIXME: check these parameters, consider numAlignments?
        self.tblastnRunner.maxTargetSeqs = 10000000
        self.tblastnRunner.maxHsps = 1
        result.paftolTargetSet.numOfftargetReads = None
        self.tblastnRunner.processTblastn(result.paftolTargetSet, self.makeWorkdirPath(self.forwardFasta), targetProteinList)
        # FIXME: should be not None (!!!)
        if result.reverseFastq is not None:
            self.tblastnRunner.processTblastn(result.paftolTargetSet, self.makeWorkdirPath(self.reverseFasta), targetProteinList)


    # ideas for hybrid / consensus sequence for (multiple) re-mapping
    # reference CDS:     atgtac------catacagaagagacgtga
    # reconstructed CDS:    cactcatttcat---gga
    # "consensus"        atgCACTCAATTCAT   GGAgagacgtga
    # principe: Where reconstructed symbol is available, use that in preference.
    #   * gap in reference: use symbols from reconstructed (must be non-gap if pairwise alignment)
    #   * gap in reconstructed: skip symbols from reference
    #   * ends / portions with no alignment to reconstructed: fill in from reference
    # Problem: avoid non-homologous alignment portions (e.g. around borders of reconstructed)?

    def analyse(self, targetsSourcePath, forwardFastq, reverseFastq, allowInvalidBases, strictOverlapFiltering, maxNumReadsPerGene):
        logger.debug('starting')
	paftolTargetSet = PaftolTargetSet()
	paftolTargetSet.readFasta(targetsSourcePath)
        # FIXME: put allowInvalidBases in result for subsequent reference?
	paftolTargetSet.sanityCheck(allowInvalidBases)
        result = HybpiperResult(paftolTargetSet, forwardFastq, reverseFastq)
	try:
            self.setup(result)
            logger.debug('setup done')
            self.mapReadsTblastn(result)
            logger.debug('tblastn mapping done')
            self.distribute(result, maxNumReadsPerGene)
            logger.debug('read distribution done')
            self.setRepresentativeGenes(result)
            self.writeRepresentativeGenes(result)
            logger.debug('representative genes selected')
            result.reconstructedCdsDict = {}
            for geneName in result.paftolTargetSet.paftolGeneDict:
                result.reconstructedCdsDict[geneName] = self.reconstructCds(result, geneName, strictOverlapFiltering)
	    logger.debug('CDS reconstruction done')
            logger.debug('finished')
            return result
        finally:
            self.makeTgz()
            logger.debug('tgz file made')
            self.cleanup()
            logger.debug('cleanup done')


class OverlapAnalyser(HybseqAnalyser):

    def __init__(self, workdirTgz=None, workDirname='pafpipertmp', tblastnRunner=None, spadesRunner=None):
        super(OverlapAnalyser, self).__init__(workdirTgz, workDirname)
        if tblastnRunner is None:
            self.tblastnRunner = paftol.tools.TblastnRunner()
        else:
            self.tblastnRunner = tblastnRunner
        self.windowSizeReference = None
        self.relIdentityThresholdReference = None
        self.windowSizeReadOverlap = None
        self.relIdentityThresholdReadOverlap = None
        # hard-coded alignment runner while API is incomplete...
        self.alignmentRunner = tools.SemiglobalAlignmentRunner()

    def setup(self, result):
        logger.debug('setting up')
        self.setupTmpdir()
	result.paftolTargetSet.writeFasta(self.makeTargetsFname(True))
        forwardFastaPath = self.makeWorkdirPath(self.forwardFasta)
        paftol.tools.fastqToFasta(result.forwardFastq, forwardFastaPath)
        self.tblastnRunner.indexDatabase(forwardFastaPath)
        if result.reverseFastq is not None:
            reverseFastaPath = self.makeWorkdirPath(self.reverseFasta)
            paftol.tools.fastqToFasta(result.reverseFastq, reverseFastaPath)
            self.tblastnRunner.indexDatabase(reverseFastaPath)

    def mapReadsTblastn(self, result):
        """Map gene sequences to reads (from multiple organisms possibly).
"""
        logger.debug('mapping gene sequences to reads')
        referenceFname = self.makeTargetsFname(True) ## check this holds
        targetProteinList = [self.translateGene(geneSr) for geneSr in result.paftolTargetSet.getSeqRecordList()]
        # FIXME: check these parameters, consider numAlignments?
        self.tblastnRunner.maxTargetSeqs = 10000000
        self.tblastnRunner.maxHsps = 1
        result.paftolTargetSet.numOfftargetReads = None
        self.tblastnRunner.processTblastn(result.paftolTargetSet, self.makeWorkdirPath(self.forwardFasta), targetProteinList)
        # FIXME: should be not None (!!!)
        if result.reverseFastq is not None:
            self.tblastnRunner.processTblastn(result.paftolTargetSet, self.makeWorkdirPath(self.reverseFasta), targetProteinList)

    def assembleGeneSerialOverlap(self, result, geneName):
        # logger.debug('tracking: starting with gene %s' % geneName)
        overlapCsvFname = self.makeWorkdirPath('overlap-%s.csv' % geneName)
        positionedReadDirname = self.makeWorkdirPath('posread-%s' % geneName)
        positionedReadFname = self.makeWorkdirPath('posread-%s.fasta' % geneName)
        os.mkdir(positionedReadDirname)
        readSrFwdList = copy.deepcopy(result.paftolTargetSet.paftolGeneDict[geneName].makeMappedReadsUniqueList(includeForward=True, includeReverse=False))
        readSrRevList = copy.deepcopy(result.paftolTargetSet.paftolGeneDict[geneName].makeMappedReadsUniqueList(includeForward=False, includeReverse=True))
        readSrList = []
        for readSr in readSrFwdList:
            readSr.id = '%s-fwd' % readSr.id
            readSrList.append(readSr)
        for readSr in readSrRevList:
            readSr.id = '%s-rev' % readSr.id
            readSrList.append(readSr)
        readSrList.extend(paftol.tools.reverseComplementSeqRecordList(readSrList))
        repGene = result.representativePaftolTargetDict[geneName].seqRecord
        # repGeneProtein = self.translateGene(repGene)
        geneReadFname = self.makeGeneReadFname(geneName)
        # for sr in readSrList:
        #     sys.stderr.write('%s\n' % sr.id)
        readSrDict = Bio.SeqIO.to_dict(readSrList)
        alignmentList = paftol.tools.semiglobalOneVsAll(repGene, readSrList)
        numReads = len(readSrList)
        if len(alignmentList) != numReads:
            raise StandardError, 'readSrList / alignment mismatch'
        positionedReadList = []
        for i in xrange(numReads):
            # sys.stderr.write('%s / %s: %f\n' % (alignmentList[i][0].id, alignmentList[i][1].id, findMaxRelativeIdentity(alignmentList[i], self.windowSize)))
            maxRelativeIdentity = paftol.tools.findMaxRelativeIdentity(alignmentList[i], self.windowSizeReference)
            if maxRelativeIdentity >= self.relIdentityThresholdReference:
                position = paftol.tools.findReadPosition(alignmentList[i])
                # coreAlignment = findOverlapAlignment(alignmentList[i])
                # coreLength = coreAlignment.get_alignment_length()
                # coreMatch = findRelativeIdentity(coreAlignment)
                coreLength = None
                coreMatch = None
                positionedReadList.append(paftol.tools.PositionedRead(readSrList[i], position, maxRelativeIdentity, coreLength, coreMatch))
            else:
                pass
                # sys.stderr.write('skipping read %s with maxRelativeIdentity %f\n' % (readSrList[i].id, maxRelativeIdentity))
        positionedReadList.sort()
        # logger.debug('tracking: positioned reads')
        positionedSrList = [positionedRead.readSr for positionedRead in positionedReadList]
        if positionedReadDirname is not None:
            for i in xrange(len(positionedSrList)):
                Bio.SeqIO.write([positionedSrList[i]], '%s/p%03d.fasta' % (positionedReadDirname, i), 'fasta')
        if positionedReadFname is not None:
            Bio.SeqIO.write(positionedSrList, positionedReadFname, 'fasta')
        if overlapCsvFname is not None:
            overlapDataFrame = paftol.tools.DataFrame(['read0', 'read1', 'read1pos', 'maxRelId', 'coreLength', 'coreMatch', 'overlapLength', 'overlapMatch'])
        else:
            overlapDataFrame = None
        contigList = []
        currentContig = paftol.tools.Contig(self.windowSizeReadOverlap, self.relIdentityThresholdReadOverlap, self.alignmentRunner)
        for i in xrange(len(positionedReadList)):
            if overlapDataFrame is not None:
                if i == 0:
                    overlapRow = {'read0': None, 'read1': positionedReadList[i].readSr.id, 'read1pos': positionedReadList[i].position, 'maxRelId': positionedReadList[i].maxRelativeIdentity, 'coreLength': positionedReadList[i].coreLength, 'coreMatch': positionedReadList[i].coreMatch, 'overlapLength': None, 'overlapMatch': None}
                else:
                    alignmentList = self.alignmentRunner.align(positionedReadList[i - 1].readSr, [positionedReadList[i].readSr])
                    alignment = alignmentList[0]
                    overlapAlignment = paftol.tools.findOverlapAlignment(alignment)
                    if overlapAlignment.get_alignment_length() == 0:
                        overlapMatch = None
                    else:
                        overlapMatch = paftol.tools.findRelativeIdentity(overlapAlignment)
                    overlapRow = {'read0': positionedReadList[i - 1].readSr.id, 'read1': positionedReadList[i].readSr.id, 'read1pos': positionedReadList[i].position, 'maxRelId': positionedReadList[i].maxRelativeIdentity, 'coreLength': positionedReadList[i].coreLength, 'coreMatch': positionedReadList[i].coreMatch, 'overlapLength': overlapAlignment.get_alignment_length(), 'overlapMatch': overlapMatch}
                overlapDataFrame.addRow(overlapRow)
            if currentContig.addRead(positionedReadList[i].readSr):
                logger.debug('added read %s to current contig', positionedReadList[i].readSr.id)
            else:
                logger.debug('started new contig with read %s', positionedReadList[i].readSr.id)
                currentContig.removeTerminalGaps()
                contigList.append(currentContig)
                currentContig = paftol.tools.Contig(self.windowSizeReadOverlap, self.relIdentityThresholdReadOverlap, self.alignmentRunner)
                currentContig.addRead(positionedReadList[i].readSr)
        currentContig.removeTerminalGaps()
        contigList.append(currentContig)
        if overlapCsvFname is not None:
            with open(overlapCsvFname, 'w') as f:
                overlapDataFrame.writeCsv(f)
        consensusList = []
        contigNumber = 0
        for contig in contigList:
            consensus = contig.getConsensus()
            if consensus is not None:
                consensus.id = '%s--contig%05d' % (geneName, contigNumber)
                contigNumber = contigNumber + 1
                consensusList.append(consensus)
        return consensusList

    def reconstructCds(self, result, geneName, strictOverlapFiltering):
        logger.debug('reconstructing CDS for gene %s', geneName)
        if result.representativePaftolTargetDict is None:
            raise StandardError('illegal state: no represesentative genes')
        if result.representativePaftolTargetDict[geneName] is None:
            raise StandardError('no representative for gene %s' % geneName)
        os.mkdir(self.makeGeneDirPath(geneName))
        contigList = self.assembleGeneSerialOverlap(result, geneName)
        if contigList is None:
            logger.warning('gene %s: no overlap contigs', geneName)
            return None
        if len(contigList) == 0:
            logger.warning('gene %s: empty overlap contig list', geneName)
            return None
        logger.debug('gene %s: %d spades contigs', geneName, len(contigList))
        geneProtein = self.translateGene(result.representativePaftolTargetDict[geneName].seqRecord)
        exonerateRunner = paftol.tools.ExonerateRunner()

        Bio.SeqIO.write([geneProtein], self.makeWorkdirPath('%s-protein.fasta' % geneName), 'fasta')
        aminoAcidSet = set(Bio.Alphabet.IUPAC.protein.letters.lower())
        # allow stop translation
        aminoAcidSet.add('*')
        setDiff = set(str(geneProtein.seq).lower()) - aminoAcidSet
        if len(setDiff) > 0:
            logger.warning('gene %s: invalid amino acids %s' % (geneName, ', '.join(setDiff)))
            return None
        contigFname = self.makeGeneContigsFname(geneName)
        Bio.SeqIO.write(contigList, contigFname, 'fasta')
        exonerateResultList = exonerateRunner.parse(geneProtein, contigFname, 'protein2genome', bestn=len(contigList))
        logger.debug('gene %s: %d contigs, %d exonerate results', geneName, len(contigList), len(exonerateResultList))
        if len(exonerateResultList) == 0:
            logger.warning('gene %s: no exonerate results from %d contigs', geneName, len(contigList))
        exonerateResultList.sort(paftol.cmpExonerateResultByQueryAlignmentStart)
        # reverse complementing extraneous as that is done by exonerate itself
        # for exonerateResult in exonerateResultList:
        #     logger.debug('gene %s, contig %s: targetStrand = %s', geneName, exonerateResult.targetId, exonerateResult.targetStrand)
        #     logger.debug('gene %s, contig %s, raw: %d -> %d, %s', geneName, exonerateResult.targetId, exonerateResult.targetCdsStart, exonerateResult.targetCdsEnd, str(exonerateResult.targetAlignmentSeq.seq))
        #     if exonerateResult.targetStrand == '-':
        #         exonerateResult.reverseComplementTarget()
        #     logger.debug('gene %s, contig %s, can: %d -> %d, %s', geneName, exonerateResult.targetId, exonerateResult.targetCdsStart, exonerateResult.targetCdsEnd, str(exonerateResult.targetAlignmentSeq.seq))
        # logger.warning('provisional filtering and supercontig construction, handling of overlapping contigs not finalised')
        filteredExonerateResultList = self.filterExonerateResultList(geneName, exonerateResultList, strictOverlapFiltering)
        logger.debug('gene %s: %d exonerate results after filtering', geneName, len(filteredExonerateResultList))
        Bio.SeqIO.write([e.targetCdsSeq for e in filteredExonerateResultList], os.path.join(self.makeGeneDirPath(geneName), '%s-fecds.fasta' % geneName), 'fasta')
        if len(filteredExonerateResultList) == 0:
            logger.warning('gene %s: no exonerate results left after filtering', geneName)
            return None
        supercontig = Bio.SeqRecord.SeqRecord(Bio.Seq.Seq(''.join([str(e.targetCdsSeq.seq) for e in filteredExonerateResultList])), id='%s_supercontig' % geneName)
        logger.debug('gene %s: supercontig length %d', geneName, len(supercontig))
        if len(supercontig) == 0:
            logger.warning('gene %s: empty supercontig', geneName)
            return None
        supercontigFname = os.path.join(self.makeGeneDirPath(geneName), '%s-supercontig.fasta' % geneName)
        Bio.SeqIO.write([supercontig], supercontigFname, 'fasta')
        Bio.SeqIO.write([geneProtein], os.path.join(self.makeGeneDirPath(geneName), '%s-supercontigref.fasta' % geneName), 'fasta')
        # FIXME: use exonerate to align "supercontig" to reference and
        # retrieve coding sequence of exonerate result with highest
        # score. In case of tied highest score, select result with
        # shortest CDS, as this is indicative of highest
        # "concentration" of matches and fewest gaps.
        supercontigErList = exonerateRunner.parse(geneProtein, supercontigFname, 'protein2genome', bestn=1)
        logger.debug('gene %s: %d supercontig exonerate results', geneName, len(supercontigErList))
        splicedSupercontigEr = None
        if len(supercontigErList) == 0:
            logger.warning('gene %s: no exonerate results from supercontig', geneName)
            return None
        if len(supercontigErList) > 1:
            splicedSupercontigEr = supercontigErList[0]
            minLength = len(splicedSupercontigEr.targetCdsSeq)
            for supercontigEr in supercontigErList:
                if len(supercontigEr.targetCdsSeq) < minLength:
                    splicedSupercontigEr = supercontigEr
                    minLength = len(splicedSupercontigEr.targetCdsSeq)
            contigStats = ', '.join(['raw=%d, cdsLen=%d' % (e.rawScore, len(e.targetCdsSeq)) for e in supercontigErList])
            logger.warning('gene %s: received %d supercontig exonerate results despite bestn=1 (%s), selected raw=%d, cdsLen=%d', geneName, len(supercontigErList), contigStats, splicedSupercontigEr.rawScore, len(splicedSupercontigEr.targetCdsSeq))
        else:
            splicedSupercontigEr = supercontigErList[0]
        # not filtering for percent identity to gene again, as that is already done
        if result.reverseFastq is not None:
            readsSpec = '%s, %s' % (result.forwardFastq, result.reverseFastq)
        else:
            readsSpec = result.forwardFastq
        splicedSupercontig = Bio.SeqRecord.SeqRecord(Bio.Seq.Seq(str(splicedSupercontigEr.targetCdsSeq.seq)), id=geneName, description='reconstructed CDS computed by paftol.overlapAnalyser, targets: %s, reads: %s' % (result.paftolTargetSet.fastaHandleStr, readsSpec))
        logger.debug('gene %s: splicedSupercontig length %d', geneName, len(splicedSupercontig))
        splicedSupercontigFname = os.path.join(self.makeGeneDirPath(geneName), '%s-splicedsupercontig.fasta' % geneName)
        Bio.SeqIO.write([splicedSupercontig], splicedSupercontigFname, 'fasta')
        return splicedSupercontig

    def analyse(self, targetsSourcePath, forwardFastq, reverseFastq, allowInvalidBases, strictOverlapFiltering, maxNumReadsPerGene):
        if self.windowSizeReference is None:
            raise StandardError, 'illegal state: windowSizeReference not set, not ready to analyse'
        if self.relIdentityThresholdReference is None:
            raise StandardError, 'illegal state: relIdentityThresholdReference not set, not ready to analyse'
        if self.windowSizeReadOverlap is None:
            raise StandardError, 'illegal state: windowSizeReadOverlap not set, not ready to analyse'
        if self.relIdentityThresholdReadOverlap is None:
            raise StandardError, 'illegal state: relIdentityThresholdReadOverlap not set, not ready to analyse'
        logger.debug('starting')
	paftolTargetSet = paftol.PaftolTargetSet()
	paftolTargetSet.readFasta(targetsSourcePath)
        # FIXME: put allowInvalidBases in result for subsequent reference?
	paftolTargetSet.sanityCheck(allowInvalidBases)
        result = paftol.HybpiperResult(paftolTargetSet, forwardFastq, reverseFastq)
	try:
            self.setup(result)
            logger.debug('setup done')
            self.mapReadsTblastn(result)
            logger.debug('tblastn mapping done')
            self.distribute(result, maxNumReadsPerGene)
            logger.debug('read distribution done')
            self.setRepresentativeGenes(result)
            self.writeRepresentativeGenes(result)
            logger.debug('representative genes selected')
            result.reconstructedCdsDict = {}
            for geneName in result.paftolTargetSet.paftolGeneDict:
                result.reconstructedCdsDict[geneName] = self.reconstructCds(result, geneName, strictOverlapFiltering)
	    logger.debug('CDS reconstruction done')
            logger.debug('finished')
            return result
        finally:
            self.makeTgz()
            logger.debug('tgz file made')
            self.cleanup()
            logger.debug('cleanup done')


class HybseqResult(object):

    def __init__(self):
        self.contigDict = None
        self.reconstructedCdsDict = None
        self.contigFastaFname = None
        self.cmdLine = None

    def summaryStats(self):
        raise StandardError, 'not implemented by this abstract class'

    def writeContigFastaFile(self, contigFastaFname):
        allContigList = []
        for contigList in self.contigDict.values():
            if contigList is not None:
                allContigList.extend(contigList)
        Bio.SeqIO.write(allContigList, contigFastaFname, 'fasta')
        self.contigFastaFname = contigFastaFname

        
class HybpiperResult(HybseqResult):

    def __init__(self, paftolTargetSet, forwardFastq, reverseFastq):
        super(HybpiperResult, self).__init__()
        self.paftolTargetSet = paftolTargetSet
        self.forwardFastq = forwardFastq
        self.reverseFastq = reverseFastq
        self.forwardFastqcStats = None
        self.reverseFastqcStats = None
        self.forwardFastqTrimmedPaired = None
        self.reverseFastqTrimmedPaired = None
        self.forwardFastqTrimmedUnpaired = None
        self.reverseFastqTrimmedUnpaired = None
        self.forwardTrimmedPairedFastqcStats = None
        self.reverseTrimmedPairedFastqcStats = None
        self.forwardTrimmedUnpairedFastqcStats = None
        self.reverseTrimmedUnpairedFastqcStats = None
        self.representativePaftolTargetDict = None
        
    def generateFastqcStats(self):
        if self.forwardFastq is not None:
            self.forwardFastqcStats = paftol.tools.generateFastqcStats(self.forwardFastq)
        if self.reverseFastq is not None:
            self.reverseFastqcStats = paftol.tools.generateFastqcStats(self.reverseFastq)
        if self.forwardFastqTrimmedPaired is not None:
            self.forwardTrimmedPairedFastqcStats = paftol.tools.generateFastqcStats(self.forwardFastqTrimmedPaired)
        if self.reverseFastqTrimmedPaired is not None:
            self.reverseTrimmedPairedFastqcStats = paftol.tools.generateFastqcStats(self.reverseFastqTrimmedPaired)
        if self.forwardFastqTrimmedUnpaired is not None:
            self.forwardTrimmedUnpairedFastqcStats = paftol.tools.generateFastqcStats(self.forwardFastqTrimmedUnpaired)
        if self.reverseFastqTrimmedUnpaired is not None:
            self.reverseTrimmedUnpairedFastqcStats = paftol.tools.generateFastqcStats(self.reverseFastqTrimmedUnpaired)

    def isPaired(self):
        return self.reverseFastq is not None

    def summaryStats(self):
        if self.reconstructedCdsDict is None:
	    raise StandardError, 'Illegal state, reconstructedCdsDict not populated'
	summaryColumnList = ['sampleName', 'targetsFile', 'paftolGene', 'paftolOrganism', 'paftolTargetLength', 'numReadsFwd', 'numReadsRev', 'qual28Fwd', 'qual28Rev', 'meanA', 'stddevA', 'meanC', 'stddevC', 'meanG', 'stddevG', 'meanT', 'stddevT', 'meanN', 'stddevN', 'numMappedReads', 'numMappedReadsPerGene', 'totNumMappedReads', 'totNumUnmappedReads', 'hybpiperCdsLength', 'representativeTarget']
        fqDataFrameFwd = paftol.tools.generateFastqcStats(self.forwardFastq)
        perBaseSequenceContentFwd = fqDataFrameFwd.calculateMeanStd(fqDataFrameFwd.perBaseSequenceContent)
        perBaseNContentFwd = fqDataFrameFwd.calculateMeanStd(fqDataFrameFwd.perBaseNContent)
        if self.reverseFastq is not None:
            fqDataFrameRev = paftol.tools.generateFastqcStats(self.reverseFastq)
            perBaseSequenceContentRev = fqDataFrameRev.calculateMeanStd(fqDataFrameRev.perBaseSequenceContent)
            perBaseNContentRev = fqDataFrameRev.calculateMeanStd(fqDataFrameRev.perBaseNContent)
        summaryDataFrame = paftol.tools.DataFrame(summaryColumnList)
        rowDict = {}
        for columnName in summaryColumnList:
            rowDict[columnName] = None
        rowDict['numReadsFwd'] = paftol.tools.countSeqRecords(self.forwardFastq, 'fastq')
        paftolSampleId = extractPaftolSampleId(self.forwardFastq)
        rowDict['sampleName'] = paftolSampleId
        rowDict['targetsFile'] = self.paftolTargetSet.fastaHandleStr
        rowDict['qual28Fwd'] = paftol.tools.getQual28(fqDataFrameFwd.perBaseSequenceQuality)
        if self.reverseFastq is not None:
            rowDict['numReadsRev'] = paftol.tools.countSeqRecords(self.reverseFastq, 'fastq')
            rowDict['qual28Rev'] = paftol.tools.getQual28(fqDataFrameRev.perBaseSequenceQuality)
            rowDict['meanA'] = (perBaseSequenceContentFwd['a'].mean + perBaseSequenceContentRev['a'].mean) / 2.0
            rowDict['stddevA'] = (perBaseSequenceContentFwd['a'].std + perBaseSequenceContentRev['a'].std) / 2.0
            rowDict['meanC'] = (perBaseSequenceContentFwd['c'].mean + perBaseSequenceContentRev['c'].mean) / 2.0
            rowDict['stddevC'] = (perBaseSequenceContentFwd['c'].std + perBaseSequenceContentRev['c'].std) / 2.0
            rowDict['meanG'] = (perBaseSequenceContentFwd['g'].mean + perBaseSequenceContentRev['g'].mean) / 2.0
            rowDict['stddevG'] = (perBaseSequenceContentFwd['g'].std + perBaseSequenceContentRev['g'].std) / 2.0
            rowDict['meanT'] = (perBaseSequenceContentFwd['t'].mean + perBaseSequenceContentRev['t'].mean) / 2.0
            rowDict['stddevT'] = (perBaseSequenceContentFwd['t'].std + perBaseSequenceContentRev['t'].std) / 2.0
            rowDict['meanN'] = (perBaseNContentFwd['nCount'].mean + perBaseNContentRev['nCount'].mean) / 2.0
            rowDict['stddevN'] = (perBaseNContentFwd['nCount'].std + perBaseNContentRev['nCount'].std) / 2.0
        else:
            rowDict['meanA'] = perBaseSequenceContentFwd['a'].mean
            rowDict['stddevA'] = perBaseSequenceContentFwd['a'].std
            rowDict['meanC'] = perBaseSequenceContentFwd['c'].mean
            rowDict['stddevC'] = perBaseSequenceContentFwd['c'].std
            rowDict['meanG'] = perBaseSequenceContentFwd['g'].mean
            rowDict['stddevG'] = perBaseSequenceContentFwd['g'].std
            rowDict['meanT'] = perBaseSequenceContentFwd['t'].mean
            rowDict['stddevT'] = perBaseSequenceContentFwd['t'].std
            rowDict['meanN'] = perBaseNContentFwd['nCount'].mean
            rowDict['stddevN'] = perBaseNContentFwd['nCount'].std
        rowDict['totNumMappedReads'] = len(self.paftolTargetSet.getMappedReadNameSet())
        rowDict['totNumUnmappedReads'] = self.paftolTargetSet.numOfftargetReads
        targetSeqRecordList = self.paftolTargetSet.getSeqRecordList()
        for paftolOrganism in self.paftolTargetSet.organismDict.values():
            rowDict['paftolOrganism'] = paftolOrganism.name
            for paftolTarget in paftolOrganism.paftolTargetDict.values():
                geneName = paftolTarget.paftolGene.name
                rowDict['paftolGene'] = paftolTarget.paftolGene.name
                rowDict['paftolTargetLength'] = len(paftolTarget.seqRecord)
                rowDict['numMappedReadsPerGene'] = len(paftolTarget.paftolGene.getReadNameSet())
                rowDict['numMappedReads'] = paftolTarget.numMappedReads()
                hybpiperCdsLength = None
                if paftolTarget.paftolGene.name in self.reconstructedCdsDict and self.reconstructedCdsDict[paftolTarget.paftolGene.name] is not None:
                    hybpiperCdsLength = len(self.reconstructedCdsDict[paftolTarget.paftolGene.name])
                rowDict['hybpiperCdsLength'] = hybpiperCdsLength
                representativeTarget = None
                if geneName in self.representativePaftolTargetDict:
                    representativeTarget = self.representativePaftolTargetDict[geneName]
                rowDict['representativeTarget'] = representativeTarget.seqRecord.id
                summaryDataFrame.addRow(rowDict)
        return summaryDataFrame


# FIXME: also extracts more generic sample IDs now, so name is not quite right
def extractPaftolSampleId(fastqName):
    paftolFastqRe = re.compile('([A-Z0-9_-]+)_L001_R([12])')
    illuminaPairedFastqRe = re.compile('([^.])_R1[^.]*\\.fastq')
    genericFastqRe = re.compile('([^.]+)\\.fastq')
    for r in [paftolFastqRe, illuminaPairedFastqRe, genericFastqRe]:
        m = r.match(fastqName)
        if m is not None:
            return m.group(1)
    raise StandardError, 'invalid PAFTOL fastq sample file name: %s' % fastqName


def paftolSummary(paftolTargetFname, fastqPairList, bwaRunner):
    summaryColumnList = ['sampleName', 'targetsFile', 'paftolGene', 'paftolOrganism', 'paftolTargetLength', 'numReadsFwd', 'numReadsRev', 'qual28Fwd', 'qual28Rev', 'meanA', 'stddevA', 'meanC', 'stddevC', 'meanG', 'stddevG', 'meanT', 'stddevT', 'meanN', 'stddevN', 'numMappedReads', 'totNumMappedReads', 'totNumUnmappedReads', 'hybpiperCdsLength']
    summaryDataFrame = paftol.tools.DataFrame(summaryColumnList)
    for fastqFwd, fastqRev in fastqPairList:
        logger.debug('fastqPair: %s, %s' % (fastqFwd, fastqRev))
        rowDict = {}
        for columnName in summaryColumnList:
            rowDict[columnName] = None
        rowDict['numReadsFwd'] = paftol.tools.countSeqRecords(fastqFwd, 'fastq')
        rowDict['numReadsRev'] = paftol.tools.countSeqRecords(fastqRev, 'fastq')
        paftolSampleId = extractPaftolSampleId(fastqFwd)
        rowDict['sampleName'] = paftolSampleId
        rowDict['targetsFile'] = paftolTargetFname
        fqDataFrameFwd = paftol.tools.generateFastqcStats(fastqFwd)
        fqDataFrameRev = paftol.tools.generateFastqcStats(fastqRev)
        rowDict['qual28Fwd'] = getQual28(fqDataFrameFwd.perBaseSequenceQuality)
        rowDict['qual28Rev'] = getQual28(fqDataFrameRev.perBaseSequenceQuality)
        perBaseSequenceContentFwd = fqDataFrameFwd.calculateMeanStd(fqDataFrameFwd.perBaseSequenceContent)
        perBaseSequenceContentRev = fqDataFrameRev.calculateMeanStd(fqDataFrameRev.perBaseSequenceContent)
        rowDict['meanA'] = (perBaseSequenceContentFwd['a'].mean + perBaseSequenceContentRev['a'].mean) / 2.0
        rowDict['stddevA'] = (perBaseSequenceContentFwd['a'].std + perBaseSequenceContentRev['a'].std) / 2.0
        rowDict['meanC'] = (perBaseSequenceContentFwd['c'].mean + perBaseSequenceContentRev['c'].mean) / 2.0
        rowDict['stddevC'] = (perBaseSequenceContentFwd['c'].std + perBaseSequenceContentRev['c'].std) / 2.0
        rowDict['meanG'] = (perBaseSequenceContentFwd['g'].mean + perBaseSequenceContentRev['g'].mean) / 2.0
        rowDict['stddevG'] = (perBaseSequenceContentFwd['g'].std + perBaseSequenceContentRev['g'].std) / 2.0
        rowDict['meanT'] = (perBaseSequenceContentFwd['t'].mean + perBaseSequenceContentRev['t'].mean) / 2.0
        rowDict['stddevT'] = (perBaseSequenceContentFwd['t'].std + perBaseSequenceContentRev['t'].std) / 2.0
        perBaseNContentFwd = fqDataFrameFwd.calculateMeanStd(fqDataFrameFwd.perBaseNContent)
        perBaseNContentRev = fqDataFrameRev.calculateMeanStd(fqDataFrameRev.perBaseNContent)
        rowDict['meanN'] = (perBaseNContentFwd['nCount'].mean + perBaseNContentRev['nCount'].mean) / 2.0
        rowDict['stddevN'] = (perBaseNContentFwd['nCount'].std + perBaseNContentRev['nCount'].std) / 2.0
        hybpiperAnalyser = HybpiperAnalyser(paftolTargetFname, fastqFwd, fastqRev, bwaRunner=bwaRunner)
        hybpiperResult = hybpiperAnalyser.analyse()
        rowDict['totNumMappedReads'] = len(hybpiperAnalyser.paftolTargetSet.getMappedReadNameSet())
        rowDict['totNumUnmappedReads'] = hybpiperAnalyser.paftolTargetSet.numOfftargetReads
        targetSeqRecordList = hybpiperAnalyser.paftolTargetSet.getSeqRecordList()
        for paftolOrganism in hybpiperAnalyser.paftolTargetSet.organismDict.values():
            rowDict['paftolOrganism'] = paftolOrganism.name
            for paftolTarget in paftolOrganism.paftolTargetDict.values():
                rowDict['paftolGene'] = paftolTarget.paftolGene.name
                rowDict['paftolTargetLength'] = len(paftolTarget.seqRecord)
                rowDict['numMappedReads'] = paftolTarget.numMappedReads()
                hybpiperCdsLength = None
                if paftolTarget.paftolGene.name in hybpiperResult.reconstructedCdsDict and hybpiperResult.reconstructedCdsDict[paftolTarget.paftolGene.name] is not None:
                    hybpiperCdsLength = len(hybpiperResult.reconstructedCdsDict[paftolTarget.paftolGene.name])
                rowDict['hybpiperCdsLength'] = hybpiperCdsLength
                summaryDataFrame.addRow(rowDict)
        tmpCsvFname = 'tmp_%s.csv' % paftolSampleId
        with open(tmpCsvFname, 'w') as tmpCsv:
            summaryDataFrame.writeCsv(tmpCsv)
    return summaryDataFrame


    summaryColumnList = ['sampleName', 'targetsFile', 'paftolGene', 'paftolOrganism', 'paftolTargetLength', 'numReadsFwd', 'numReadsRev', 'qual28Fwd', 'qual28Rev', 'meanA', 'stddevA', 'meanC', 'stddevC', 'meanG', 'stddevG', 'meanT', 'stddevT', 'meanN', 'stddevN', 'numMappedReads', 'numMappedReadsPerGene', 'totNumUnmappedReads', 'hybpiperCdsLength', 'representativeTarget']
    geneSetStatsDataFrame = paftol.tools.DataFrame(summaryColumnList)
    srDict = Bio.SeqIO.to_dict(Bio.SeqIO.parse(f, 'fasta'))
    for paftolGene in paftolTargetSet.paftolGeneDict.values():
        for paftolTarget in paftolGene.paftolTargetDict.values():
            rowDict = {}
            for summaryColumn in summaryColumnList:
                rowDict[summaryColumn] = None
            rowDict['sampleName'] = sampleName
            rowDict['paftolGene'] = paftolGene.name
            rowDict['paftolOrganism'] = paftolTarget.organism.name
            rowDict['targetsFile'] = paftolTargetSet.fastaHandleStr
            rowDict['paftolTargetLength'] = len(paftolTarget.seqRecord)
            if paftolGene.name in srDict:
                rowDict['hybpiperCdsLength'] = len(srDict[paftolGene.name])
            geneSetStatsDataFrame.addRow(rowDict)
    return geneSetStatsDataFrame
