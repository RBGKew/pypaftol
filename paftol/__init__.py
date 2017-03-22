import sys
import re
import os
import tempfile
import subprocess
import shutil
import multiprocessing

import Bio
import Bio.SeqIO
import Bio.SeqIO.QualityIO
import Bio.Alphabet.IUPAC

import paftol.tools


verbose = 0
keepTmp = False


def isSane(filename):
    """Check whether a file name is sane, in the sense that it does not contain any "funny" characters"""
    if filename == '':
        return False
    funnyCharRe = re.compile('/ ;,$#')
    m = funnyCharRe.search(filename)
    if m is not None:
        return False
    if filename[0] == '-':
        return False
    return True

            
def cmpExonerateResultByQueryAlignmentStart(e0, e1):
    if e0.queryAlignmentStart < e1.queryAlignmentStart:
        return -1
    elif e0.queryAlignmentStart > e1.queryAlignmentStart:
        return 1
    return 0

    
class HybseqAnalyser(object):
    
    def __init__(self, targetsSourcePath, outFname, forwardFastq, reverseFastq=None, workdirTgz=None, workDirname='paftoolstmp'):
        self.targetsSourcePath = targetsSourcePath
        self.outFname = outFname
        self.forwardFastq = forwardFastq
        self.reverseFastq = reverseFastq
        self.workdirTgz = workdirTgz
        self.workDirname = workDirname
        self.tmpDirname = None
        # parameters for ensuring file names don't clash, e.g. because locus / organism name is same as targets basename etc.
        self.targetsFname = 'targets.fasta'
        self.locusFnamePattern = 'locus-%s.fasta'
    
    def __str__(self):
        return 'HybseqAnalyser(targetsSourcePath=%s, outFname=%s, forwardFastq=%s, reverseFastq=%s)' % (repr(self.targetsSourcePath), repr(self.outFname), repr(self.forwardFastq), repr(self.reverseFastq))
        
    def checkTargets(self):
        # FIXME: merge with __init__()?
        for targetSr in Bio.SeqIO.parse(self.targetsSourcePath, 'fasta', alphabet = Bio.Alphabet.IUPAC.ambiguous_dna):
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
                sys.stderr.write('not removing temporary directory %s\n' % self.tmpDirname)
            else:
                shutil.rmtree(self.tmpDirname)
            self.tmpDirname = None

    def makeWorkDirname(self):
        if self.tmpDirname is None:
            raise StandardError, 'illegal state: no temporary directory and hence no working directory'
        # sys.stderr.write('tmpDirname = %s, workDirname = %s\n' % (self.tmpDirname, self.workDirname))
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
to provide fields required for HybSeq analysis only."""
    
    def __init__(self, samLine):
        if samLine[-1] == '\n':
            samLine = samLine[:-1]
        w = samLine.split('\t')
        self.qname = w[0]
        self.rname = w[2]
        self.mapq = int(w[4])
        self.seq = w[9]


class OrganismLocus(object):
    
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
    
    def __init__(self, name):
        self.name = name
        self.organismLocusDict = {}


class Locus(object):
    
    def __init__(self, name):
        self.name = name
        self.organismLocusDict = {}

    def qnameSet(self):
        s = set()
        for organismLocus in self.organismLocusDict.values():
            s = s | organismLocus.qnameSet()
        return s
    

class HybpiperAnalyser(HybseqAnalyser):
    
    organismLocusRe = re.compile('([^-]+)-([^-]+)')
    
    def __init__(self, targetsSourcePath, outFname, forwardFastq, reverseFastq=None, workdirTgz=None, workDirname='pafpipertmp'):
        super(HybpiperAnalyser, self).__init__(targetsSourcePath, outFname, forwardFastq, reverseFastq, workdirTgz, workDirname)
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
        self.representativeOrganismLocusDict = None

    def setup(self):
        if self.targetsSourcePath is None:
            raise StandardError, 'illegal state: cannot set up with targetsSourcePath = None'
        self.setupTmpdir()
        shutil.copy(self.targetsSourcePath, self.makeTargetsFname(True))
            
    def cleanup(self):
        self.cleanupTmpdir()
    
    def bwaIndexReference(self):
        bwaIndexArgv = ['bwa', 'index', self.makeTargetsFname(True)]
        if verbose:
            sys.stderr.write('%s\n' % ' '.join(bwaIndexArgv))
        subprocess.check_call(bwaIndexArgv)
        
    def mapReadsBwa(self):
        self.bwaIndexReference()
        fastqArgs = [os.path.join(os.getcwd(), self.forwardFastq)]
        if self.reverseFastq is not None:
            fastqArgs.append(os.path.join(os.getcwd(), self.reverseFastq))
        # bwa parameters for tweaking considerations: -k, -r, -T
        bwaArgv = ['bwa', 'mem', '-M', '-k', '%d' % self.bwaMinSeedLength, '-T', '%d' % self.bwaScoreThreshold, '-r', '%f' % self.bwaReseedTrigger, '-t', '%d' % self.bwaNumThreads, self.makeTargetsFname()] + fastqArgs
        samtoolsArgv = ['samtools', 'view', '-h', '-S', '-F', '4']
        sys.stderr.write('%s\n' % ' '.join(bwaArgv))
        bwaProcess = subprocess.Popen(bwaArgv, stdout=subprocess.PIPE, cwd = self.makeWorkDirname())
        sys.stderr.write('%s\n' % ' '.join(samtoolsArgv))
        samtoolsProcess = subprocess.Popen(samtoolsArgv, stdin=bwaProcess.stdout.fileno(), stdout=subprocess.PIPE, cwd = self.makeWorkDirname())
        samLine = samtoolsProcess.stdout.readline()
        while samLine != '':
            # sys.stderr.write(samLine)
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
        """Roughly equivalent to "distribute targets" in HybPiper"""
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
                sys.stderr.write('represenative for %s: none\n' % locusName)
            else:
                sys.stderr.write('representative for %s: %s\n' % (representativeOrganismLocus.locus.name, representativeOrganismLocus.organism.name))
    
    def distributeSingle(self):
        fForward = open(self.forwardFastq, 'r')
        fqiForward = Bio.SeqIO.QualityIO.FastqGeneralIterator(fForward)
        for fwdReadTitle, fwdReadSeq, fwdReadQual in fqiForward:
            readName = fwdReadTitle.split()[0]
            for locus in self.locusDict.values():
                if readName in locus.qnameSet():
                    f = open(self.makeLocusFname(locus.name, True), 'a')
                    sys.stderr.write('appending to %s\n' % f.name)
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
        sys.stderr.write('%s\n' % ' '.join(parallelSpadesArgv))
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
            sys.stderr.write('locus fasta file %s does not exist (no reads?)\n' % locusFname)
            return False
        spadesArgv = ['spades.py', '--only-assembler', '--threads', '1', '--cov-cutoff', '%d' % self.spadesCovCutoff]
        if self.spadesKvalList is not None:
            spadesArgv.extend(['-k', ','.join(['%d' % k for k in self.spadesKvalList])])
        spadesArgv.extend(spadesInputArgs)
        spadesArgv.extend(['-o', self.makeLocusDirname(locusName)])
        sys.stderr.write('%s\n' % ' '.join(spadesArgv))
        spadesProcess = subprocess.Popen(spadesArgv, cwd = self.makeWorkDirname())
        spadesReturncode = spadesProcess.wait()
        if spadesReturncode != 0:
            # raise StandardError, 'spades process "%s" exited with %d' % (' '.join(spadesArgv), spadesReturncode)
            sys.stderr.write('spades process "%s" exited with %d\n' % (' '.join(spadesArgv), spadesReturncode))
        spadesContigFname = os.path.join(self.makeLocusDirPath(locusName), 'contigs.fasta')
        # sys.stderr.write('spadesContigFname: %s\n' % spadesContigFname)
        if os.path.exists(spadesContigFname):
            spadesContigList = list(Bio.SeqIO.parse(spadesContigFname, 'fasta'))
            # sys.stderr.write('spadesContigFname: %s, %d contigs\n' % (spadesContigFname, len(spadesContigList)))
        else:
            spadesContigList = None
            # sys.stderr.write('spadesContigFname: %s, no contigs\n' % spadesContigFname)
        return spadesContigList
    
    def translateLocus(self, locusDna):
        # FIXME: add support for locus specific translation table setting
        l = len(locusDna) - (len(locusDna) % 3)
        if l < len(locusDna):
            sys.stderr.write('locus %s: length %d is not an integer multiple of 3 -- not a CDS?\n' % (locusDna.id, len(locusDna)))
        locusProtein = Bio.SeqRecord.SeqRecord(locusDna.seq[:l].translate(), id='%s-pep' % locusDna.id, description='%s, translated' % locusDna.description)
        return locusProtein

    def exonerateContigs(self, locusName, contigList):
        if self.representativeOrganismLocusDict is None:
            raise StandardError, 'illegal state: no represesentative loci'
        if self.representativeOrganismLocusDict[locusName] is None:
            raise StandardError, 'no representative for locus %s' % locusName
        locusProtein = self.translateLocus(self.representativeOrganismLocusDict[locusName].seqRecord)
        contigFname = os.path.join(self.makeLocusDirPath(locusName), '%s-contigs.fasta' % locusName)
        Bio.SeqIO.write(contigList, contigFname, 'fasta')
        exonerateRunner = paftol.tools.ExonerateRunner()
        # FIXME: need to translate query
        exonerateResultList = exonerateRunner.parse(locusProtein, contigFname, 'protein2genome', len(contigList))
        sys.stderr.write('%d contigs, %d exonerate results\n' % (len(contigList), len(exonerateResultList)))
        exonerateResultList.sort(cmpExonerateResultByQueryAlignmentStart)
        for exonerateResult in exonerateResultList:
            sys.stderr.write('  %s\n' % str(exonerateResult))
        return exonerateResultList
    
    def filterByPercentIdentity(self, exonerateResultList):
        return [e for e in exonerateResultList if e.percentIdentity >= self.exoneratePercentIdentityThreshold]
    
    def filterByContainment(self, exonerateResultList):
        
        def isContainedWithTiebreak(exonerateResult, other):
            if not other.containsQueryAlignmentRange(exonerateResult):
                return False
            if not exonerateResult.containsQueryAlignmentRange(other):
                return True
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
    
    def filterByOverlap(self, exonerateResultList):
        nonOverlappingExonerateResultList = []
        for exonerateResult in exonerateResultList:
            for other in exonerateResultList:
                if exonerateResult is not other:
                    if exonerateResult.overlapsQueryAlignmentRange(other):
                        sys.stderr.write('  overlap found: %s, %s\n' % (str(exonerateResult), str(other)))
            nonOverlappingExonerateResultList.append(exonerateResult)
        return nonOverlappingExonerateResultList
        
    def analyse(self):
        self.checkTargets()
        try:
            self.setup()
            self.mapReadsBwa()
            self.distribute()
            self.setRepresentativeLoci()
            reconstructedLocusDict = {}
            for locusName in self.locusDict:
                os.mkdir(self.makeLocusDirPath(locusName))
                spadesContigList = self.assembleLocusSpades(locusName)
                if spadesContigList is None:
                    sys.stderr.write('locus %s: no spades contigs\n' % locusName)
                    reconstructedLocusDict[locusName] = None
                else:
                    sys.stderr.write('locus %s: %d spades contigs\n' % (locusName, len(spadesContigList)))
                    exonerateResultList = self.exonerateContigs(locusName, spadesContigList)
                    sys.stderr.write('locus %s: %d exonerate results\n' % (locusName, len(exonerateResultList)))
                    exonerateResultList = self.filterByPercentIdentity(exonerateResultList)
                    sys.stderr.write('locus %s: %d sufficiently close exonerate results\n' % (locusName, len(exonerateResultList)))
                    exonerateResultList = self.filterByContainment(exonerateResultList)
                    sys.stderr.write('locus %s: %d non-contained exonerate results\n' % (locusName, len(exonerateResultList)))
                    exonerateResultList = self.filterByOverlap(exonerateResultList)
                    sys.stderr.write('locus %s: %d non-overlapping exonerate results\n' % (locusName, len(exonerateResultList)))
                    # TODO: concatenate and 'exonerate' to remove introns (?)
            self.makeTgz()
        finally:
            self.cleanup()
