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
    
    
class HybseqAnalyser(object):
    
    def __init__(self, targetsFname, outFname, forwardFastq, reverseFastq=None, workdirTgz=None, workDirname='paftoolstmp'):
        self.targetsFname = targetsFname
        self.outFname = outFname
        self.forwardFastq = forwardFastq
        self.reverseFastq = reverseFastq
        self.workdirTgz = workdirTgz
        self.workDirname = workDirname
        self.tmpDirname = None
        self.numThreads = 1
    
    def __str__(self):
        return 'HybseqAnalyser(targetsFname=%s, outFname=%s, forwardFastq=%s, reverseFastq=%s)' % (repr(self.targetsFname), repr(self.outFname), repr(self.forwardFastq), repr(self.reverseFastq))
        
    def checkTargets(self):
        for targetSr in Bio.SeqIO.parse(self.targetsFname, 'fasta', alphabet = Bio.Alphabet.IUPAC.ambiguous_dna):
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
            # shutil.rmtree(self.tmpDirname)
            sys.stderr.write('not removing temporary directory %s\n' % self.tmpDirname)
            self.tmpDirname = None

    def makeWorkDirname(self):
        if self.tmpDirname is None:
            raise StandardError, 'illegal state: no temporary directory and hence no working directory'
        # sys.stderr.write('tmpDirname = %s, workDirname = %s\n' % (self.tmpDirname, self.workDirname))
        return os.path.join(self.tmpDirname, self.workDirname)
            
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
            shutil.move(os.path.join(self.tmpDirname, tmpTgz), self.workdirTgz)


class HybpiperAnalyser(HybseqAnalyser):
    
    taxonLocusRe = re.compile('([^-]+)-([^-]+)')
    
    def __init__(self, targetsFname, outFname, forwardFastq, reverseFastq=None, workdirTgz=None, workDirname='pafpipertmp'):
        super(HybpiperAnalyser, self).__init__(targetsFname, outFname, forwardFastq, reverseFastq, workdirTgz, workDirname)
        self.initLocusNameSet()
        
    def extractLocusName(self, seqId, allowNew=False):
        m = self.taxonLocusRe.match(seqId)
        if m is not None:
            locusName = m.group(2)
        else:
            locusName = seqId
        if not allowNew and locusName not in self.locusNameSet:
            raise StandardError, 'unknown locus "%s" (found in seqId "%s")' % (locusName, seqId)
        return locusName
        
    def initLocusNameSet(self):
        if self.targetsFname is None:
            raise StandardError, 'illegal state: cannot init locus name with targetsFname = None'
        self.locusNameSet = set()
        for sr in Bio.SeqIO.parse(self.targetsFname, 'fasta'):
            self.locusNameSet.add(self.extractLocusName(sr.id, allowNew=True))

    def setup(self):
        self.setupTmpdir()
        shutil.copy(self.targetsFname, self.makeWorkDirname())
            
    def cleanup(self):
        self.cleanupTmpdir()
    
    def bwaIndexReference(self):
        bwaIndexArgv = ['bwa', 'index', os.path.join(self.makeWorkDirname(), self.targetsFname)]
        sys.stderr.write('%s\n' % ' '.join(bwaIndexArgv))
        subprocess.check_call(bwaIndexArgv)
        
    def bwaReadLocusDict(self):
        self.bwaIndexReference()
        fastqArgs = [os.path.join(os.getcwd(), self.forwardFastq)]
        if self.reverseFastq is not None:
            fastqArgs.append(os.path.join(os.getcwd(), self.reverseFastq))
        # bwa parameters for tweaking considerations: -k, -r, -T
        bwaArgv = ['bwa', 'mem', '-M', '-k', '11', '-T', '10', '-t', '%d' % self.numThreads, self.targetsFname] + fastqArgs
        samtoolsArgv = ['samtools', 'view', '-h', '-S', '-F', '4']
        sys.stderr.write('%s\n' % ' '.join(bwaArgv))
        bwaProcess = subprocess.Popen(bwaArgv, stdout=subprocess.PIPE, cwd = self.makeWorkDirname())
        sys.stderr.write('%s\n' % ' '.join(samtoolsArgv))
        samtoolsProcess = subprocess.Popen(samtoolsArgv, stdin=bwaProcess.stdout.fileno(), stdout=subprocess.PIPE, cwd = self.makeWorkDirname())
        readLocusDict = {}
        samLine = samtoolsProcess.stdout.readline()
        while samLine != '':
            # sys.stderr.write(samLine)
            if samLine[0] != '@':
                w = samLine[:-1].split('\t')
                qname = w[0]
                rname = w[2]
                locusName = self.extractLocusName(rname)
                if qname not in readLocusDict:
                    readLocusDict[qname] = set()
                readLocusDict[qname].add(locusName)
            samLine = samtoolsProcess.stdout.readline()
        bwaProcess.stdout.close()
        samtoolsProcess.stdout.close()
        bwaReturncode = bwaProcess.wait()
        samtoolsReturncode = samtoolsProcess.wait()
        if bwaReturncode != 0:
            raise StandardError, 'process "%s" returned %d' % (' '.join(bwaArgv), bwaReturncode)
        if samtoolsReturncode != 0:
            raise StandardError, 'process "%s" returned %d' % (' '.join(samtoolsArgv), samtoolsReturncode)
        return readLocusDict
    
    def locusFastaPath(self, locusName, suffix = ''):
        return os.path.join(self.makeWorkDirname(), '%s%s.fasta' % (locusName, suffix))
    
    def locusInterleavedFastaPath(self, locusName):
        return os.path.join(self.makeWorkDirname(), '%s_interleaved.fasta' % locusName)
    
    def distributeSingle(self, readLocusDict):
        fForward = open(self.forwardFastq, 'r')
        fqiForward = Bio.SeqIO.QualityIO.FastqGeneralIterator(fForward)
        for fwdReadTitle, fwdReadSeq, fwdReadQual in fqiForward:
            readName = fwdReadTitle.split()[0]
            if readName in readLocusDict:
                for locusName in readLocusDict[readName]:
                    f = open(self.locusFastaPath(locusName), 'a')
                    f.write('>%s\n%s\n' % (fwdReadTitle, fwdReadSeq))
                    f.close()
        fForward.close()
    
    def distributePaired(self, readLocusDict):
        # FIXME: consider try...finally to ensure files are closed
        for readName in readLocusDict.keys():
            sys.stderr.write('%s\n' % readName)
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
            if readName in readLocusDict:
                for locusName in readLocusDict[readName]:
                    f = open(self.locusInterleavedFastaPath(locusName), 'a')
                    f.write('>%s\n%s\n' % (fwdReadTitle, fwdReadSeq))
                    f.write('>%s\n%s\n' % (revReadTitle, revReadSeq))
                    f.close()
                    # f = open(self.locusFastaPath(locusName, '_R1'), 'a')
                    # f.write('>%s\n%s\n' % (fwdReadTitle, fwdReadSeq))
                    # f.close()
                    # f = open(self.locusFastaPath(locusName, '_R2'), 'a')
                    # f.write('>%s\n%s\n' % (revReadTitle, revReadSeq))
                    # f.close()
        # should trigger an exception: revReadTitle, revReadSeq, revReadQual = fqiReverse.next()
        fForward.close()
        fReverse.close()
    
    def distributeBwa(self):
        readLocusDict = self.bwaReadLocusDict()
        if self.isPaired():
            self.distributePaired(readLocusDict)
        else:
            self.distributeSingle(readLocusDict)
            
            
    def assembleSpades(self, spadesNumThreads=1, spadesCovCutoff=8, spadesKvals=None):
        # consider --fg to ensure wait for all parallel processes?
        # is --eta really of any use here?
        # FIXME: hard-coded fasta pattern '{}_interleaved.fasta' for parallel
        if self.isPaired():
            spadesInputArgs = ['--12', '{}_interleaved.fasta']
        else:
            spadesInputArgs = ['-s', '{}.fasta']
            
        parallelSpadesArgv = ['parallel', 'spades.py', '--only-assembler', '--threads', '%d' % spadesNumThreads, '--cov-cutoff', '%d' % spadesCovCutoff]
        if spadesKvals is not None:
            parallelSpadesArgv.extend(['-k', ','.join(spadesKvals)])
        parallelSpadesArgv.extend(spadesInputArgs)
        parallelSpadesArgv.extend(['-o', '{}_spades'])
        # time parallel --eta spades.py --only-assembler --threads 1 --cov-cutoff 8 --12 {}/{}_interleaved.fasta -o {}/{}_spades :::: spades_genelist.txt > spades.log
        parallelSpadesProcess = subprocess.Popen(parallelSpadesArgv, stdin=subprocess.PIPE, cwd = self.makeWorkDirname())
        sys.stderr.write('%s\n' % ' '.join(parallelSpadesArgv))
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
    
    def analyse(self):
        self.checkTargets()
        try:
            self.setup()
            self.distributeBwa()
            self.assembleSpades()
            self.makeTgz()
        finally:
            self.cleanup()


def runHybseq(argNamespace):
    hybseqAnalyser = HybseqAnalyser(argNamespace.targetsfile, argNamespace.outfile, argNamespace.forwardreads, argNamespace.reversereads, argNamespace.tgz)
    sys.stderr.write('%s\n' % str(hybseqAnalyser))


def runHybpiper(argNamespace):
    hybpiperAnalyser = HybpiperAnalyser(argNamespace.targetsfile, argNamespace.outfile, argNamespace.forwardreads, argNamespace.reversereads, argNamespace.tgz)
    hybpiperAnalyser.analyse()
