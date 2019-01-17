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
        #raise StandardError, 'abstract method'
        p = self.makeSubprocess()
        argv = self.makeArgv()
        logger.debug('%s', ' '.join(argv))
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
        returncode = p.wait()
        if returncode != 0:
            raise StandardError, 'process "%s" returned %d' % (' '.join(argv), returncode)
        return alignment

    def makeSubprocess(self):
        argv = self.makeArgv()
        p = subprocess.Popen(argv, stdin=subprocess.PIPE, stdout=subprocess.PIPE, universal_newlines=True)
        return p


class MafftRunner(MultipleSequenceAlignmentRunner):

    def __init__(self):
        self.localpair = None
        self.maxiterate = None

    def makeArgv(self):
        argv = ['mafft', '--quiet']
        if self.localpair is not None:
            argv.append('--localpair')
        if self.maxiterate is not None:
            argv.append('--maxiterate')
            if self.maxiterate is not None:
                argv.append('1000')
        argv.append('-')
        return argv


class ClustaloRunner(MultipleSequenceAlignmentRunner):

    def __init__(self):
        pass

    def makeArgv(self):
        argv = ['clustalo', '-i', '-']
        return argv
