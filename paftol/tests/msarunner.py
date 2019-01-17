import sys
import unittest
import subprocess

import Bio
import Bio.Seq
import Bio.SeqRecord
import Bio.Align

import paftol
import paftol.msarunner


class MultipleSequenceAlignmentRunnerTestCase(unittest.TestCase):
    
    def test_instantiation(self):
        multipleSequenceAlignmentRunner = paftol.msarunner.MultipleSequenceAlignmentRunner()
        self.assertIsNotNone(multipleSequenceAlignmentRunner)

        
class MafftRunnerTestCase(unittest.TestCase):
    
    # def test_subprocess(self):
    #     mafftRunner = paftol.msarunner.MafftRunner()
    #     p = mafftRunner.makeSubprocess()
    #     self.assertIsInstance(p, subprocess.Popen, 'p is not a Popen object')
        
    def test_align(self):
        seqRecordList = [
            Bio.SeqRecord.SeqRecord(Bio.Seq.Seq('CGTGATACATTACTTTTTA'), id='seq0', description=''),
            Bio.SeqRecord.SeqRecord(Bio.Seq.Seq('GTGGACTTGACGCGTCATGGAAAGTACAAGATACTTCGACCTGGCAGTGCAAG'), id='seq1', description='')
        ]
        sys.stderr.write('\n')
        Bio.SeqIO.write(seqRecordList, sys.stderr, 'fasta')
        mafftRunner = paftol.msarunner.MafftRunner()
        alignment = mafftRunner.align(seqRecordList)
        self.assertIsInstance(alignment, Bio.Align.MultipleSeqAlignment)


    def test_mafftattributes(self):
        mafftRunner = paftol.msarunner.MafftRunner()
        self.assertIsNone(mafftRunner.localpair)
        self.assertTrue('--localpair' not in mafftRunner.makeMafftArgv())
        mafftRunner.localpair = True
        self.assertTrue('--localpair' in mafftRunner.makeMafftArgv())
        self.assertIsNone(mafftRunner.maxiterate)
        self.assertTrue('--maxiterate' not in mafftRunner.makeMafftArgv())
        mafftRunner.maxiterate = 1000
        self.assertTrue('1000' in mafftRunner.makeMafftArgv())
        self.assertTrue('--maxiterate' in mafftRunner.makeMafftArgv())
        
