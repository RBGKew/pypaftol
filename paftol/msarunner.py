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
        Bio.SeqIO.write(seqRecordList, p.stdin, 'fasta')
        p.stdin.close()
        alignment = Bio.AlignIO.parse(p.stdout, 'fasta')
        return alignment

    def makeMafftSubprocess(self, sequenceList):
        p = subprocess.Popen(['mafft', '-'], stdin=subprocess.PIPE, stdout=subprocess.PIPE, universal_newlines=True)
        return p

    def makeMafftArgv(self, mafftProgram, ):
        mafftArgv = [mafftProgram]

    def processMafft(self, mafftProgram, makeMafftSubprocess):
        mafftArgv = self.makeMafftArgv(mafftProgram)  # databaseFname?)
        logger.debug('%s', ' '.join(mafftArgv))
        mafftProcess = subprocess.Popen(mafftArgv, stdin=subprocess.PIPE, stdout=subprocess.PIPE)
        pid = os.fork()
        if pid == 0:
            mafftProcess.stdout.close()
            for sequence in sequenceList:
                mafftProcess.stdin.write(query.format('fasta'))
            mafftProcess.stdin.close()
            os._exit(0)
        mafftProcess.stdin.close()
        for alignedSequence in alignedSequenceList: #change alignedSequenceList to smthg including mafftProcess.stdout
            for alignment in 
                mafftAlignmentProcessor.processMafftAlignment(mafftRecord.query, alignment)
        mafftProcess.stdout.close()
        wPid, wExit = os.waitpid(pid, 0)
        if pid != wPid:
            raise StandardError, 'wait returned pid %s (expected %d)' % (wPid, pid)
        if wExit != 0:
            raise StandardError, 'wait on forked process returned %d' % wExit
        mafftReturncode = mafftProcess.wait()
        if mafftReturncode != 0:
            raise StandardError, 'process "%s" returned %d' % (' '.join(mafftArgv), mafftReturncode)

class ClustaloRunner(MultipleSequenceAlignmentRunner):

    def __init__(self):
        pass

    def clustalo():
        mergedSequences = mergeSequencesAndConvertToFasta(sequenceList, fastaFile)
        p = subprocess.Popen(['clustalo' '-i'], stdin=subprocess.PIPE, stdout=subprocess.PIPE, universal_newlines=True)
