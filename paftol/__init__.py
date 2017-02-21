import sys
import os
import tempfile
import subprocess
import shutil
import multiprocessing

import Bio

verbose = 0
    
    
class HybseqAnalyser(object):
    
    def __init__(self, targetsFname, outFname, forwardFastq, reverseFastq=None, tempDirname=None):
        self.targetsFname = targetsFname
        self.outFname = outFname
        self.forwardFastq = forwardFastq
        self.reverseFastq = reverseFastq
        self.tempDirname = tempDirname
        self.generateTempdir = tempDirname is None
    
    def __str__(self):
        return 'HybseqAnalyser(targetsFname=%s, outFname=%s, forwardFastq=%s, reverseFastq=%s)' % (repr(self.targetsFname), repr(self.outFname), repr(self.forwardFastq), repr(self.reverseFastq))
        
    def isPaired(self):
        return self.reverseFastq is not None
    
    def analyse(self):
        raise StandardError, 'not implemented in this "abstract" base class'
    

class HybpiperAnalyser(HybseqAnalyser):
    
    def __init__(self, targetsFname, outFname, forwardFastq, reverseFastq=None, tempDirname=None):
        super(HybpiperAnalyser, self).__init__(targetsFname, outFname, forwardFastq, reverseFastq, tempDirname)

    def setup(self):
        if self.generateTempdir:
            self.tempDirname = tempfile.mkdtemp(prefix='paftools')
        shutil.copy(self.targetsFname, self.tempDirname)
            
    def cleanup(self):
        if self.generateTempdir:
            # shutil.rmtree(self.tempDirname)
            sys.stderr.write('not removing temporary directory %s\n' % self.tempDirname)
    
    def bwaIndexReference(self):
        bwaIndexArgs = ['bwa', 'index', os.path.join(self.tempDirname, self.targetsFname)]
        subprocess.check_call(bwaIndexArgs)
        
    def distributeBwa(self):
        self.bwaIndexReference()
        
    def analyse(self):
        self.setup()
        self.distributeBwa()
        self.cleanup()
        
    
def runHybseq(argNamespace):
    hybseqAnalyser = HybseqAnalyser(argNamespace.targetsfile, argNamespace.outfile, argNamespace.forwardreads, argNamespace.reversereads)
    sys.stderr.write('%s\n' % str(hybseqAnalyser))


def runHybpiper(argNamespace):
    hybpiperAnalyser = HybpiperAnalyser(argNamespace.targetsfile, argNamespace.outfile, argNamespace.forwardreads, argNamespace.reversereads)
    hybpiperAnalyser.analyse()
