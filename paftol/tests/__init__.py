import sys
import unittest
import copy

import Bio
import Bio.Seq
import Bio.SeqRecord

import paftol
import paftol.tools

# atggattacagca
# atggattacagca

def getConsensusFromCopy(contig):
    contigCopy = copy.deepcopy(contig)
    contigCopy.removeTerminalGaps()
    return contigCopy.getConsensus()


class PaftolTestCase(unittest.TestCase):

    def setUp(self):
        self.seq0 = Bio.SeqRecord.SeqRecord(Bio.Seq.Seq('CGTGATACATTACTTTTTA'), id='seq0', description='x1a')
        self.seq1 = Bio.SeqRecord.SeqRecord(Bio.Seq.Seq('GTGGACTTGACGCGTCATGGAAAGTACAAGATACTTCGACCTGGCAGTGCAAG'), id='seq1', description='x1b')

    def test_SemiglobalAlignmentRunner(self):
        semiglobalAlignmentRunner = paftol.tools.SemiglobalAlignmentRunner()
        aList = semiglobalAlignmentRunner.align(self.seq0, [self.seq1])
        a = aList[0]
        self.assertEqual('------------CGTGA-------TACA--TTACTTTTTA-----------------', str(a[0].seq))
        self.assertEqual('GTGGACTTGACGCGTCATGGAAAGTACAAGATACTT----CGACCTGGCAGTGCAAG', str(a[1].seq))

    def test_Contig(self):
        alignmentRunner = paftol.tools.SemiglobalAlignmentRunner()
        contig = paftol.tools.Contig(5, 0.7, alignmentRunner)
        r = Bio.SeqRecord.SeqRecord(Bio.Seq.Seq('gattaca', alphabet=Bio.Alphabet.IUPAC.ambiguous_dna), id='r1')
        contig.addRead(r)
        self.assertEqual('gattaca', str(getConsensusFromCopy(contig).seq), 'unexpected consensus after read %s' % r.id)
        r = Bio.SeqRecord.SeqRecord(Bio.Seq.Seq('atggatta', alphabet=Bio.Alphabet.IUPAC.ambiguous_dna), id='r0')
        contig.addRead(r)
        self.assertEqual('atggattaca', str(getConsensusFromCopy(contig).seq), 'unexpected consensus after read %s' % r.id)
        r = Bio.SeqRecord.SeqRecord(Bio.Seq.Seq('gattaca', alphabet=Bio.Alphabet.IUPAC.ambiguous_dna), id='r1a')
        contig.addRead(r)
        self.assertEqual('atggattaca', str(getConsensusFromCopy(contig).seq), 'unexpected consensus after read %s' % r.id)
        r = Bio.SeqRecord.SeqRecord(Bio.Seq.Seq('ttaca', alphabet=Bio.Alphabet.IUPAC.ambiguous_dna), id='r2')
        contig.addRead(r)
        self.assertEqual('atggattaca', str(getConsensusFromCopy(contig).seq), 'unexpected consensus after read %s' % r.id)
        r = Bio.SeqRecord.SeqRecord(Bio.Seq.Seq('ttacagc', alphabet=Bio.Alphabet.IUPAC.ambiguous_dna), id='r3')
        contig.addRead(r)
        self.assertEqual('atggattacagc', str(getConsensusFromCopy(contig).seq), 'unexpected consensus after read %s' % r.id)
        r = Bio.SeqRecord.SeqRecord(Bio.Seq.Seq('ttaccagc', alphabet=Bio.Alphabet.IUPAC.ambiguous_dna), id='r4')
        contig.addRead(r)
        self.assertEqual('atggattacagc', str(getConsensusFromCopy(contig).seq), 'unexpected consensus after read %s' % r.id)
        r = Bio.SeqRecord.SeqRecord(Bio.Seq.Seq('ttacgc', alphabet=Bio.Alphabet.IUPAC.ambiguous_dna), id='r5')
        contig.addRead(r)
        self.assertEqual('atggattacagc', str(getConsensusFromCopy(contig).seq), 'unexpected consensus after read %s' % r.id)
        r = Bio.SeqRecord.SeqRecord(Bio.Seq.Seq('gattacca', alphabet=Bio.Alphabet.IUPAC.ambiguous_dna), id='r6')
        contig.addRead(r)
        self.assertEqual('atggattacagca', str(getConsensusFromCopy(contig).seq), 'unexpected consensus after read %s' % r.id)
        contig.removeTerminalGaps()
        consensus = contig.getConsensus()
        alignment = contig.getAlignment(terminalGapChar='~')
        contigDepthProfile = contig.getDepthProfile()
        consensusDepthProfile = consensus.letter_annotations['depth']
        self.assertEqual([1, 1, 1, 4, 4, 8, 8, 8, 3, 5, 5, 3, 4, 1], contigDepthProfile, 'unexpected contig depth profile: %s' % str(contigDepthProfile))
        self.assertEqual([1, 1, 1, 4, 4, 8, 8, 8, 5, 5, 3, 4, 1], consensusDepthProfile, 'unexpected consensus depth profile: %s' % str(consensusDepthProfile))
