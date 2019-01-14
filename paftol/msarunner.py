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
        p = subprocess.Popen(['mafft', '--auto', '--reorder', '-'], stdin=subprocess.PIPE, stdout=subprocess.PIPE, universal_newlines=True)
        Bio.SeqIO.write(seqRecordList, p.stdin, 'fasta')
        p.stdin.close()
        alignment = Bio.AlignIO.parse(p.stdout, 'fasta')
        return alignment

    def makeSubprocess(self, sequenceList):
       p = subprocess.Popen(stdin=subprocess.PIPE, stdout=subprocess.PIPE)
       pid = os.fork()
       mafftProcess = 
       if pid == 0:
        mafftProcess.stdout.close()
        for sequence in sequenceList:
            mafftProcess.stdin.write(query.format('fasta'))
        mafftProcess.stdin.close()
        os._exit(0)
        mafftProcess.stdin.close()
        for alignedSequence in alignedSequenceList:


class ClustaloRunner(MultipleSequenceAlignmentRunner):

    def __init__(self):
        pass

    def clustalo():
        mergedSequences = mergeSequencesAndConvertToFasta(sequenceList, fastaFile)
        p = subprocess.Popen(['clustalo' '-i'], stdin=subprocess.PIPE, stdout=subprocess.PIPE, universal_newlines=True)
