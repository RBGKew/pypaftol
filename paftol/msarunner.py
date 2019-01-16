import sys
import Bio.SeqIO
import csv
import subprocess
import Bio.SeqRecord
import unittest
import Bio.AlignIO
import logging
import os

logger = logging.getLogger(__name__)

class MultipleSequenceAlignmentRunner(object):

    """Wrapper class for running Multiple Sequence Alignment (MSA) programs
     such as Clustal Omega and MAFFT.

This is a base class for runners that wrap specific MSA programs.

"""

    def __init__(self):
        pass

    def align(self, seqRecordList):
        raise StandardError, 'abstract method'


class MafftRunner(MultipleSequenceAlignmentRunner):

    def __init__(self):
        pass

    def align(self, seqRecordList):
        p = self.makeMafftSubprocess()
        mafftArgv = self.makeMafftArgv()
        logger.debug('%s', ' '.join(mafftArgv))
        pid = os.fork()
        if pid == 0:
            p.stdout.close()
            Bio.SeqIO.write(seqRecordList, p.stdin, 'fasta')
            p.stdin.close()
            os._exit(0)
        p.stdin.close()
        alignment = Bio.AlignIO.read(p.stdout, 'fasta')
        p.stdout.close()
        wPid, wExit = os.waitpid(pid, 0)
        if pid != wPid:
            raise StandardError, 'wait returned pid %s (expected %d)' % (wPid, pid)
        if wExit != 0:
            raise StandardError, 'wait on forked process returned %d' % wExit
        mafftReturncode = p.wait()
        if mafftReturncode != 0:
            raise StandardError, 'process "%s" returned %d' % (' '.join(mafftArgv), mafftReturncode)
        return alignment

    def makeMafftSubprocess(self):
        mafftArgv = self.makeMafftArgv()
        p = subprocess.Popen(mafftArgv, stdin=subprocess.PIPE, stdout=subprocess.PIPE, universal_newlines=True)
        return p

    def makeMafftArgv(self):
        mafftArgv = ['mafft', '--quiet', '-']
        return mafftArgv
    
class ClustaloRunner(MultipleSequenceAlignmentRunner):

    def __init__(self):
        pass

    def align(self, seqRecordList):
        p = self.makeClustaloSubprocess()
        clustaloArgv = self.makeClustaloArgv()
        logger.debug('%s', ' '.join(clustaloArgv))
        pid = os.fork()
        if pid == 0:
            p.stdout.close()
            Bio.SeqIO.write(seqRecordList, p.stdin, 'fasta')
            p.stdin.close()
            os._exit(0)
        p.stdin.close()
        alignment = Bio.AlignIO.read(p.stdout, 'fasta')
        p.stdout.close()
        wPid, wExit = os.waitpid(pid, 0)
        if pid != wPid:
            raise StandardError, 'wait returned pid %s (expected %d)' % (wPid, pid)
        if wExit != 0:
            raise StandardError, 'wait on forked process returned %d' % wExit
        clustaloReturncode = p.wait()
        if clustaloReturncode != 0:
            raise StandardError, 'process "%s" returned %d' % (' '.join(clustaloArgv), clustaloReturncode)
        return alignment

    def makeClustaloSubprocess(self):
    clustaloArgv = self.makeClustaloArgv()
    p = subprocess.Popen(clustaloArgv, stdin=subprocess.PIPE, stdout=subprocess.PIPE, universal_newlines=True)
        return p

    def makeClustaloArgv(self):
        clustaloArgv = ['clustalo', '-i', '-']
        return clustaloArgv
        
