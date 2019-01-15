import sys
import Bio.SeqIO
import csv
import subprocess
import Bio.SeqRecord
import unittest
import Bio.AlignIO


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
        pid = os.fork()
        if pid == 0:
            p.stdout.close()
            for sequence in sequenceList:
                p.stdin.write(query.format('fasta'))
            p.stdin.close()
            os._exit(0)
        p.stdin.close()
        wPid, wExit = os.waitpid(pid, 0)
        if pid != wPid:
            raise StandardError, 'wait returned pid %s (expected %d)' % (wPid, pid)
        if wExit != 0:
            raise StandardError, 'wait on forked process returned %d' % wExit
        mafftReturncode = p.wait()
        if mafftReturncode != 0:
            raise StandardError, 'process "%s" returned %d' % (' '.join(mafftArgv), mafftReturncode)
        Bio.SeqIO.write(seqRecordList, p.stdin, 'fasta')
        p.stdin.close()
        alignment = Bio.AlignIO.read(p.stdout, 'fasta')
        return alignment

    def makeMafftSubprocess(self):
        mafftArgv = self.makeMafftArgv()
        logger.debug('%s', ' '.join(mafftArgv))
        p = subprocess.Popen(mafftArgv, stdin=subprocess.PIPE, stdout=subprocess.PIPE, universal_newlines=True)
        return p

    def makeMafftArgv(self):
        mafftArgv = ['mafft', '-']
        return mafftArgv
    
class ClustaloRunner(MultipleSequenceAlignmentRunner):

    def __init__(self):
        pass

    def clustalo():
        mergedSequences = mergeSequencesAndConvertToFasta(sequenceList, fastaFile)
        p = subprocess.Popen(['clustalo' '-i'], stdin=subprocess.PIPE, stdout=subprocess.PIPE, universal_newlines=True)
