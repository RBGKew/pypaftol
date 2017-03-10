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


verbose = 0


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

    
class HybseqAnalyser(object):
    
    def __init__(self, targetsSourcePath, outFname, forwardFastq, reverseFastq=None, workdirTgz=None, workDirname='paftoolstmp'):
        self.targetsSourcePath = targetsSourcePath
        self.outFname = outFname
        self.forwardFastq = forwardFastq
        self.reverseFastq = reverseFastq
        self.workdirTgz = workdirTgz
        self.workDirname = workDirname
        self.tmpDirname = None
        self.keepTmpDir = False
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
            if self.keepTmpDir:
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
        self.seqRecord = seqRecord
        self.samAlignmentList = []
        if locus.name in organism.organismLocusDict or organism.name in locus.organismLocusDict:
            raise StandardError, 'duplicate organism/locus: organism = %s, locus = %s, seqId = %s' % (organism.name, locus.name, seqRecord.id)
        organism.organismLocusDict[locus.name] = self
        locus.organismLocusDict[organism.name] = self
        
    def addSamAlignment(self, samAlignment):
        self.samAlignmentList.append(samAlignment)
            
    def mapqSum(self):
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
            
    def assembleSpadesParallel(self, spadesCovCutoff=8, spadesKvals=None):
        # consider --fg to ensure wait for all parallel processes?
        # is --eta really of any use here?
        # FIXME: hard-coded fasta pattern '{}_interleaved.fasta' for parallel
        if self.isPaired():
            spadesInputArgs = ['--12', '{}_interleaved.fasta']
        else:
            spadesInputArgs = ['-s', '{}.fasta']
        parallelSpadesArgv = ['parallel', 'spades.py', '--only-assembler', '--threads', '1', '--cov-cutoff', '%d' % spadesCovCutoff]
        if spadesKvals is not None:
            parallelSpadesArgv.extend(['-k', ','.join(spadesKvals)])
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
        
    def makeSpadesLocusDirname(self, locusName):
        return 'spades-%s' % locusName
            
    def assembleLocusSpades(self, locusName, spadesCovCutoff, spadesKvals):
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
        spadesArgv = ['spades.py', '--only-assembler', '--threads', '1', '--cov-cutoff', '%d' % spadesCovCutoff]
        if spadesKvals is not None:
            spadesArgv.extend(['-k', ','.join(spadesKvals)])
        spadesArgv.extend(spadesInputArgs)
        spadesArgv.extend(['-o', self.makeSpadesLocusDirname(locusName)])
        sys.stderr.write('%s\n' % ' '.join(spadesArgv))
        spadesProcess = subprocess.Popen(spadesArgv, cwd = self.makeWorkDirname())
        spadesReturncode = spadesProcess.wait()
        if spadesReturncode != 0:
            # raise StandardError, 'spades process "%s" exited with %d' % (' '.join(spadesArgv), spadesReturncode)
            sys.stderr.write('spades process "%s" exited with %d\n' % (' '.join(spadesArgv), spadesReturncode))

    def assembleSpades(self, spadesCovCutoff=8, spadesKvals=None):
        for locusName in self.locusDict:
            self.assembleLocusSpades(locusName, spadesCovCutoff, spadesKvals)
        
    def exonerate(self):
        pass
        
    def analyse(self):
        self.checkTargets()
        try:
            self.setup()
            self.mapReadsBwa()
            self.distribute()
            self.assembleSpades()
            self.exonerate()
            self.makeTgz()
        finally:
            self.cleanup()
