#!/usr/bin/env python

import types
import sys
import re
import os
import os.path
import subprocess
import csv
import tempfile
import logging
import copy
import math

import Bio
import Bio.Alphabet
import Bio.Alphabet.IUPAC
import Bio.Seq
import Bio.SeqRecord
import Bio.Align
import Bio.SeqIO
import Bio.AlignIO
import Bio.Data.CodonTable
import Bio.Blast
import Bio.Blast.NCBIXML

import paftol
import paftol.clib


logger = logging.getLogger(__name__)


def fastqToFasta(fastqFname, fastaFname):
    with open(fastqFname, 'r') as fastqFile:
        with open(fastaFname, 'w') as fastaFile:
            for record in Bio.SeqIO.parse(fastqFile, 'fastq'):
                Bio.SeqIO.write(record, fastaFile, 'fasta')


def alignCdsByProtein(alignedProteinSeqRecord, cdsSeqRecord):
    u = str(alignedProteinSeqRecord.seq)
    c = str(cdsSeqRecord.seq)
    # if len(c) != 3 * len(u):
    #     raise StandardError, 'incompatible sequences: CDS %s (length %d), ungapped prtotein sequence %s (length %d)' % (cdsSeqRecord.id, len(c), alignedProteinSeqRecord.id, len(u))
    proteinAlphabet = alignedProteinSeqRecord.seq.alphabet
    if not isinstance(proteinAlphabet, Bio.Alphabet.Gapped):
        proteinAlphabet = Bio.Alphabet.Gapped(proteinAlphabet)
    cdsAlphabet = cdsSeqRecord.seq.alphabet
    if not isinstance(cdsAlphabet, Bio.Alphabet.Gapped):
        cdsAlphabet = Bio.Alphabet.Gapped(cdsAlphabet)
    a = str(alignedProteinSeqRecord.seq)
    s = ''
    i = 0
    for aa in a:
        if aa == proteinAlphabet.gap_char:
            s = s + 3 * cdsAlphabet.gap_char
        else:
            s = s + c[i:i + 3]
            i = i + 3
    return Bio.SeqRecord.SeqRecord(Bio.Seq.Seq(s, cdsAlphabet), id=cdsSeqRecord.id, description='%s, gapped along %s' % (cdsSeqRecord.description, alignedProteinSeqRecord.id))


def countSeqRecords(fName, fFormat):
    count = 0
    for seqRecord in Bio.SeqIO.parse(fName, fFormat):
        count = count + 1
    return count

def translateGapped(seq, table='Standard'):
    """Translate a gapped sequence.

All triplets in the sequence must either contain only non-gap symbols
or only gap symbols. Runs of non-gap triplets are translated as per
the specified translation table. Runs of gap symbol triplets are
translated into as many single gap symbols. An exception is raised
if the sequence contains triplets with both gap and non-gap symbols.

@param seq: sequence to be translated
@type seq: C{Bio.Seq.Seq}
@param table: translation table, passed to BioPython C{translate} method
@type table: C{int} or C{str} (any type suitable for C{translate})
@return: translated sequence
@rtype: C{Bio.Seq.Seq}
"""
    if not isinstance(seq.alphabet, Bio.Alphabet.Gapped):
        return seq.translate(table=table)
    ungappedAlphabet = seq.alphabet.alphabet
    s = str(seq)
    gapStart = None
    gapList = []
    for i in xrange(len(s)):
        if s[i] == seq.alphabet.gap_char and gapStart is None:
            gapStart = i
        if s[i] != seq.alphabet.gap_char and gapStart is not None:
            gapList.append((gapStart, i, ))
            gapStart = None
    if gapStart is not None:
        gapList.append((gapStart, len(s), ))
    # logger.debug('gaplist: %s', str(gapList))
    t = ''
    nongapStart = 0
    for gapStart, gapEnd in gapList:
        if (gapStart - nongapStart) % 3 != 0:
            raise StandardError('cannot translate: range %d:%d does not coincide with triplet boundaries' % (nongapStart, gapStart))
        if (gapEnd - gapStart) % 3 != 0:
            raise StandardError('cannot translate: range %d:%d does not coincide with triplet boundaries' % (gapStart, gapEnd))
        if nongapStart < gapStart:
            ungappedSeq = seq[nongapStart:gapStart]
            ungappedSeq.alphabet = ungappedAlphabet
            t = t + str(ungappedSeq.translate(table))
        t = t + seq.alphabet.gap_char * ((gapEnd - gapStart) / 3)
        nongapStart = gapEnd
    if nongapStart < len(seq):
        if (len(seq) - nongapStart) % 3 != 0:
            raise StandardError('cannot translate: range %d:%d does not coincide with triplet boundaries' % (nongapStart, len(seq)))
        ungappedSeq = seq[nongapStart:]
        ungappedSeq.alphabet = ungappedAlphabet
        t = t + str(ungappedSeq.translate(table))
    return Bio.Seq.Seq(t, alphabet=Bio.Alphabet.Gapped(Bio.Alphabet.IUPAC.protein))


def ascendingRange(rangeStart, rangeEnd):
    """Get start and end into ascending order by switching them if necessary.
@param rangeStart: the start of the range
@type rangeStart: C{int}
@param rangeEnd: the end of the range
@type rangeEnd: C{int}
@return: range in ascending order
@rtype: C{tuple} of length 2
"""
    if rangeStart > rangeEnd:
        return rangeEnd, rangeStart
    else:
        return rangeStart, rangeEnd


class BlastAlignment(object):

    def __init__(self, rname, blastAlignment):
        raise StandardError, 'obsolete'
        #self.query = query
        self.rname = rname
        self.qname = blastAlignment.hit_id
        self.mapq = blastAlignment.hsps[0].score
        #self.qnameSet = set([alignment.hit_id for alignment in blastAlignment.alignments])
        self.expectValue = blastAlignment.hsps[0].expect

    def isMapped(self):
        return True


class SamAlignment(object):

    """Class to represent a SAM record.
This class follows the naming and definitions of the SAMv1 spec. It is incomplete
to provide fields required for Hyb-Seq analysis only.

@ivar qname: SAM query name (C{QNAME}), read or read pair ID
@type qname: C{str}
@ivar rname: SAM reference name (C{RNAME})
@type rname: C{str}
@ivar flag: SAM flag (C{FLAG})
@type flag: C{int}
@ivar pos: SAM mapping position (C{POS})
@type pos: C{int}
@ivar mapq: SAM mapping quality (C{MAPQ})
@type mapq: C{int}
@ivar cigar: SAM CIGAR string (unexpanded) (C{CIGAR})
@type mapq: C{str}
@ivar seq: SAM query (read) sequence (C{SEQ})
@type seq: C{str}
"""

    cigarElementRe = re.compile('([0-9]+)([MIDNSHP=X])')

    def __init__(self, samLine):
        if samLine[-1] == '\n':
            samLine = samLine[:-1]
        w = samLine.split('\t')
        self.qname = w[0]
        self.flag = int(w[1])
        self.rname = w[2]
        self.pos = int(w[3])
        self.mapq = int(w[4])
        self.cigar = w[5]
        self.seq = w[9]

    def isMapped(self):
        return self.flag & 4 == 0

    def getMatchLength(self):
        e = self.expandedCigar()
        return e.count('M') + e.count('D')

    def getEndpos(self):
        return self.pos + self.getMatchLength()

    def expandedCigar(self):
        if self.cigar is None:
            return None
        e = ''
        c = self.cigar
        while c != '':
            m = self.cigarElementRe.match(c)
            if m is None:
                raise StandardError('malformed CIGAR "%s" (stuck at "%s")' % (self.cigar, c))
            e = e + (m.group(2) * int(m.group(1)))
            c = c[len(m.group()):]
        return e

    def numCigarMatches(self):
        e = self.expandedCigar()
        if e is None:
            return None
        if e.count('=') > 0:
            logger.warning('found sequence match ("=") characters, unimplemented')
        if e.count('X') > 0:
            logger.warning('found sequence mismatch ("X") characters, unimplemented')
        return e.count('M')


class ExonerateResult(object):
    """Hold results from running C{exonerate}.

Many attributes correspond to conversions provided by C{exonerate}'s
roll your own (ryo) formatting facility.

@ivar querySeq: query sequence (first free parameter to C{exonerate})
@type querySeq: C{Bio.SeqRecord.SeqRecord}
@ivar targetSeq: target sequence file (second free parameter to C{exonerate})
@type targetSeq: C{Bio.SeqRecord.SeqRecord}, or C{None} if adding of target sequences was not requested
@ivar targetFname: name of a FASTA file of sequences to be scanned
@type targetFname: C{str}
@ivar exonerateModel: alignment model, currently only C{protein2genome:local} is really supported
@type exonerateModel: C{str}
@ivar queryId: query ID (C{%qi} conversion)
@type queryId: C{str}
@ivar queryDef: query definition (C{%qd} conversion)
@type queryDef: C{str}
@ivar queryStrand: query strand (C{%qS} conversion)
@type queryStrand: C{str}
@ivar queryAlignmentStart: query alignment start (C{%qab} conversion)
@type queryAlignmentStart: C{int}
@ivar queryAlignmentEnd: query alignment end (C{%qae} conversion)
@type queryAlignmentEnd: C{int}
@ivar queryAlignmentLength: query alignment length (C{%qal} conversion)
@type queryAlignmentLength: C{int}
@ivar queryCdsStart: query coding sequence start (C{%qcb} conversion)
@type queryCdsStart: C{int}
@ivar queryCdsEnd: query coding sequence end (C{%qce} conversion)
@type queryCdsEnd: C{int}
@ivar queryCdsLength: query coding sequence length (C{%qcl} conversion)
@type queryCdsLength: C{int}
@ivar targetId: target ID (C{%ti} conversion)
@type targetId: C{str}
@ivar targetDef: target definition (C{%td} conversion)
@type targetDef: C{str}
@ivar targetStrand: target strand (C{%tS} conversion)
@type targetStrand: C{str}
@ivar targetAlignmentStart: target alignment start (C{%tab} conversion)
@type targetAlignmentStart: C{int}
@ivar targetAlignmentEnd: target alignment end (C{%tae} conversion)
@type targetAlignmentEnd: C{int}
@ivar targetAlignmentLength: target alignment length (C{%tal} conversion)
@type targetAlignmentLength: C{int}
@ivar targetCdsStart: target coding sequence start (C{%tcb} conversion)
@type targetCdsStart: C{int}
@ivar targetCdsEnd: target coding sequence end (C{%tce} conversion)
@type targetCdsEnd: C{int}
@ivar targetCdsLength: target coding sequence length (C{%tcl} conversion)
@type targetCdsLength: C{int}
@ivar rawScore: raw score (C{%s} conversion)
@type rawScore: C{int}
@ivar percentIdentity: percent identity (C{%pi} conversion)
@type percentIdentity: C{float}
@ivar percentSimilarity: percent similarity (C{%ps} conversion)
@type percentSimilarity: C{float}
@ivar equivalencedTotal: equivalenced total (C{%et} conversion)
@type equivalencedTotal: C{int}
@ivar equivalencedIdentity: equivalenced identity (C{%ei} conversion)
@type equivalencedIdentity: C{int}
@ivar equivalencedSimilarity: equivalenced similarity (C{%es} conversion)
@type equivalencedSimilarity: C{int}
@ivar equivalencedMismatches: equivalenced mismatched (C{%em} conversion)
@type equivalencedMismatches: C{int}
@ivar queryAlignmentSeq: query alignment sequence (C{%qas} conversion)
@type queryAlignmentSeq: C{Bio.SeqRecord.SeqRecord}
@ivar queryCdsSeq: query coding sequence (C{%qcs} conversion)
@type queryCdsSeq: C{Bio.SeqRecord.SeqRecord}
@ivar targetAlignmentSeq: target alignment sequence (C{%tas} conversion)
@type targetAlignmentSeq: C{Bio.SeqRecord.SeqRecord}
@ivar targetCdsSeq: target coding sequence (C{%tcs} conversion)
@type targetCdsSeq: C{Bio.SeqRecord.SeqRecord}
    """
    queryRe = re.compile('Query: ([^ ]+)( (.*))?')
    targetRe = re.compile('Target: ([^ ]+)( (.*))?')
    queryRangeRe = re.compile('Query range: ([0-9]+) -> ([0-9]+)')
    targetRangeRe = re.compile('Target range: ([0-9]+) -> ([0-9]+)')
    rawScoreRe = re.compile('rawScore: ([0-9]+)')
    seqlineRe = re.compile('[ACGT]*')

    def __init__(self, querySeq, targetFname):
        self.querySeq = querySeq
        self.targetFname = targetFname
        self.targetSeq = None
        self.exonerateModel = None
        self.queryId = None
        self.queryDef = None
        self.queryStrand = None
        self.queryAlignmentStart = None
        self.queryAlignmentEnd = None
        self.queryAlignmentLength = None
        self.queryCdsStart = None
        self.queryCdsEnd = None
        self.queryCdsLength = None
        self.targetId = None
        self.targetDef = None
        self.targetStrand = None
        self.targetAlignmentStart = None
        self.targetAlignmentEnd = None
        self.targetAlignmentLength = None
        self.targetCdsStart = None
        self.targetCdsEnd = None
        self.targetCdsLength = None
        self.rawScore = None
        self.percentIdentity = None
        self.percentSimilarity = None
        self.equivalencedTotal = None
        self.equivalencedIdentity = None
        self.equivalencedSimilarity = None
        self.equivalencedMismatches = None
        self.queryAlignmentSeq = None
        self.queryCdsSeq = None
        self.targetAlignmentSeq = None
        self.targetCdsSeq = None

    def reverseComplementTarget(self):
        """Reverse complement target sequences, switching them from reverse to forward direction or vice versa.

Ranges and strand orientations are changed accordingly.
"""
        if self.targetStrand not in '+-':
            raise StandardError('cannot reverse complement strand with orientation "%s"' % self.targetStrand)
        if self.targetAlignmentSeq is not None:
            self.targetAlignmentSeq = self.targetAlignmentSeq.reverse_complement(id=True, name=True, description=True)
            self.targetAlignmentStart, self.targetAlignmentEnd = self.targetAlignmentEnd, self.targetAlignmentStart
        if self.targetCdsSeq is not None:
            self.targetCdsSeq = self.targetCdsSeq.reverse_complement(id=True, name=True, description=True)
            self.targetCdsStart, self.targetCdsEnd = self.targetCdsEnd, self.targetCdsStart
        self.targetStrand = '+' if self.targetStrand == '-' else '-'

    def queryAlignmentOverlap(self, other):
        """Find overlap between query alignment range in this and another exonerate result.

Ranges are canonicalised to be ascending, therefore returned ranges are ascending too.
@param other: the other exonerate result
@type other: L{ExonerateResult}
@return: the range of the overlap if any, or else C{None}
@rtype: C{tuple} of length 2
"""
        # FIXME: check whether self and other pertain to same query (also for containment)?
        selfQas, selfQae = ascendingRange(self.queryAlignmentStart, self.queryAlignmentEnd)
        otherQas, otherQae = ascendingRange(other.queryAlignmentStart, other.queryAlignmentEnd)
        if selfQas > otherQae or selfQae < otherQas:
            return None
        return max(selfQas, otherQas), min(selfQae, otherQae)

    def containsQueryAlignmentRange(self, other):
        """Determine whether the query alignment range in this result completely contains the range in another result.

@param other: the other exonerate result
@type other: L{ExonerateResult}
@return: C{True} if this result's alignment range fully contains that in C{other}, else C{False}
@rtype: C{bool}
"""
        selfQas, selfQae = ascendingRange(self.queryAlignmentStart, self.queryAlignmentEnd)
        otherQas, otherQae = ascendingRange(other.queryAlignmentStart, other.queryAlignmentEnd)
        return selfQas <= otherQas and selfQae >= otherQae

    def overlapsQueryAlignmentRange(self, other):
        """Determine whether the query alignment range in this result overlaps with the range in another result.

@param other: the other exonerate result
@type other: L{ExonerateResult}
@return: C{True} if this result's query alignment range overlaps with that in C{other}, else C{False}
@rtype: C{bool}
"""
        return self.queryAlignmentOverlap(other) is not None

    def nucleotideAlignment(self, appendFlanking=False):
        """Construct an alignment of the protein query and target sequences in this result.
@return: the protein sequence alignment
@rtype: C{Bio.Align.MultipleSeqAlignment}
"""
        #FIXME: take gapChar from alphabet or parameter?
        gapChar = '-'
        if self.exonerateModel != 'affine:local:dna2dna':
            raise StandardError('nucleotideAlignment is not supported for exonerate model "%s"' % self.exonerateModel)
        v = self.vulgar.split()
        qAln = ''
        tAln = ''
        vulgarLetterAnnotation = []
        qPos = 0
        tPos = 0
        qSeq = str(self.querySeq.seq)
        qLeftFlank = qSeq[:self.queryAlignmentStart].lower()
        qAligned = qSeq[self.queryAlignmentStart:self.queryAlignmentEnd].upper()
        if qAligned != str(self.queryAlignmentSeq.seq).upper():
            sys.stderr.write('%s:\nqseq: %s\nqaln: %s\n' % (self.queryId, qAligned, str(self.querySeq.seq)))
            raise StandardError, 'query sequence and query alignment sequence mismatch'
        qRightFlank = qSeq[self.queryAlignmentEnd:].lower()
        # sys.stderr.write('%s: tas = %d, tae = %d, strand = %s\n' % (self.targetId, self.targetAlignmentStart, self.targetAlignmentEnd, self.targetStrand))
        if self.targetStrand == '+':
            tSeq = str(self.targetSeq.seq)
            tLeftFlank = tSeq[:self.targetAlignmentStart].lower()
            tAligned = tSeq[self.targetAlignmentStart:self.targetAlignmentEnd].upper()
            # sys.stderr.write('targetAlignmentStart = %d, targetAlignmentEnd = %d, len(tSeq) = %d, len(tAligned) = %d, tAligned = tSeq[%d:%d]\n' % (self.targetAlignmentStart, self.targetAlignmentEnd, len(tSeq), len(tAligned), self.targetAlignmentStart, self.targetAlignmentEnd))
            # sys.stderr.write('%s:\ntseq: %s\ntaln: %s\nqaln: %s\n' % (self.targetId, tSeq, tAligned, qAligned))
            if tAligned != str(self.targetAlignmentSeq.seq).upper():
                sys.stderr.write('%s:\ntseq: %s\ntaln: %s\nqaln: %s\n' % (self.targetId, tSeq, tAligned, qAligned))
                raise StandardError, 'target sequence and target alignment sequence mismatch'
            tRightFlank = tSeq[self.targetAlignmentEnd:].lower()
        elif self.targetStrand == '-':
            tSeq = str(self.targetSeq.reverse_complement().seq)
            rcTas = len(tSeq) - self.targetAlignmentStart
            rcTae = len(tSeq) - self.targetAlignmentEnd
            tLeftFlank = tSeq[:rcTas].lower()
            tAligned = tSeq[rcTas:rcTae].upper()
            # tAligned = str(self.targetSeq[self.targetAlignmentEnd:self.targetAlignmentStart].reverse_complement().seq).upper()
            # sys.stderr.write('targetAlignmentStart = %d, targetAlignmentEnd = %d, rcTas = %d, rcTae = %d, len(tSeq) = %d, len(tAligned) = %d, tAligned = tSeq[%d:%d]\n' % (self.targetAlignmentStart, self.targetAlignmentEnd, rcTas, rcTae, len(tSeq), len(tAligned), self.targetAlignmentEnd, self.targetAlignmentStart))
            # sys.stderr.write('%s:\ntseq: %s\ntaln: %s\nqaln: %s\n' % (self.targetId, tSeq, tAligned, qAligned))
            if tAligned != str(self.targetAlignmentSeq.seq).upper():
                # sys.stderr.write('targetAlignmentStart = %d, targetAlignmentEnd = %d, rcTas = %d, rcTae = %d, len(tSeq) = %d, len(tAligned) = %d, tAligned = tSeq[%d:%d]\n' % (self.targetAlignmentStart, self.targetAlignmentEnd, rcTas, rcTae, len(tSeq), len(tAligned), self.targetAlignmentEnd, self.targetAlignmentStart))
                # sys.stderr.write('%s:\ntseq: %s\ntaln: %s\nqaln: %s\n' % (self.targetId, tSeq, tAligned, qAligned))
                raise StandardError, 'target sequence and target alignment sequence mismatch'
            tRightFlank = tSeq[rcTae:].lower()
        else:
            raise StandardError, 'need target strand orientation +/-'
        # logger.debug('queryAlignmentSeq: %d, targetAlignmentSeq: %d', len(self.queryAlignmentSeq), len(self.targetAlignmentSeq))
        # logger.debug('%s', str(self.targetAlignmentSeq.seq))
        for i in xrange(0, len(v), 3):
            vLabel = v[i]
            vQueryLength = int(v[i + 1])
            vTargetLength = int(v[i + 2])
            # logger.debug('qPos = %d, tPos = %d, vLabel = %s, vql = %d, vtl = %d', qPos, tPos, vLabel, vQueryLength, vTargetLength)
            if vLabel == 'M':
                qAln = qAln + qAligned[qPos:(qPos + vQueryLength)]
                tAln = tAln + tAligned[tPos:(tPos + vTargetLength)]
            elif vLabel == 'G':
                if vQueryLength == 0:
                    qAln = qAln + gapChar * vTargetLength
                    tAln = tAln + tAligned[tPos:(tPos + vTargetLength)]
                elif vTargetLength == 0:
                    qAln = qAln + qAligned[qPos:(qPos + vQueryLength)]
                    tAln = tAln + gapChar * vQueryLength
            else:
                raise StandardError, 'unsupported VULGAR label: %s' % vLabel
            vLength = max(vQueryLength, vTargetLength)
            vulgarLetterAnnotation.extend([vLabel] * vLength)
            qPos = qPos + vQueryLength
            tPos = tPos + vTargetLength
            # logger.debug('%d: qPos = %d, tPos = %d, vLabel = %s, vql = %d, vtl = %d', i, qPos, tPos, vLabel, vQueryLength, vTargetLength)
            # logger.debug('qAln: %s', qAln)
            # logger.debug('tAln: %s', tAln)
            # s = Bio.Seq.Seq(tAln, alphabet=Bio.Alphabet.Gapped(self.targetAlignmentSeq.seq.alphabet))
            # s3 = s[:(len(s) - len(s) % 3)]
            # logger.debug('tA_tr: %s%s', str(translateGapped(s3)), '' if len(s) == len(s3) else '.')
        #FIXME: should not unconditionally assume unambiguous_dna
        # sys.stderr.write('len(qAln) = %d, len(tAln) = %d, len(vulgarLetterAnnotation) = %d\n' % (len(qAln), len(tAln), len(vulgarLetterAnnotation)))
        if appendFlanking:
            lfDiff = len(qLeftFlank) - len(tLeftFlank)
            if lfDiff > 0:
                tLeftFlank = gapChar * lfDiff + tLeftFlank
            elif lfDiff < 0:
                qLeftFlank = gapChar * (-lfDiff) + qLeftFlank
            rfDiff = len(qRightFlank) - len(tRightFlank)
            if rfDiff > 0:
                tRightFlank = tRightFlank + gapChar * rfDiff
            elif rfDiff < 0:
                qRightFlank = qRightFlank + gapChar * (-rfDiff)
            qAln = qLeftFlank + qAln + qRightFlank
            tAln = tLeftFlank + tAln + tRightFlank
            vulgarLetterAnnotation = ([None] * len(qLeftFlank)) + vulgarLetterAnnotation + ([None] * len(qRightFlank))
        # logger.debug('len(qAln) = %d, len(tAln) = %d, len(vulgarLetterAnnotation) = %d', len(qAln), len(tAln), len(vulgarLetterAnnotation))
        # sys.stderr.write('len(qAln) = %d, len(tAln) = %d, len(vulgarLetterAnnotation) = %d\n' % (len(qAln), len(tAln), len(vulgarLetterAnnotation)))
        qAlnSeq = Bio.Seq.Seq(qAln, alphabet=Bio.Alphabet.Gapped(Bio.Alphabet.IUPAC.unambiguous_dna))
        tAlnSeq = Bio.Seq.Seq(tAln, alphabet=Bio.Alphabet.Gapped(self.targetAlignmentSeq.seq.alphabet))
        # FIXME: assuming standard translation table -- check whether exonerate supports setting table?
        return Bio.Align.MultipleSeqAlignment([Bio.SeqRecord.SeqRecord(qAlnSeq, id=self.queryAlignmentSeq.id, letter_annotations = {'vulgar': vulgarLetterAnnotation[:]}), Bio.SeqRecord.SeqRecord(tAlnSeq, id=self.targetAlignmentSeq.id, letter_annotations = {'vulgar': vulgarLetterAnnotation[:]})])

    def proteinAlignment(self):
        """Construct an alignment of the protein query and target sequences in this result.
@return: the protein sequence alignment
@rtype: C{Bio.Align.MultipleSeqAlignment}
"""
        gapChar = '-'
        if self.exonerateModel != 'protein2genome:local':
            raise StandardError('proteinAlignment is not supported for exonerate model "%s"' % self.exonerateModel)
        v = self.vulgar.split()
        qAln = ''
        tAln = ''
        qPos = 0
        tPos = 0
        # logger.debug('queryAlignmentSeq: %d, targetAlignmentSeq: %d', len(self.queryAlignmentSeq), len(self.targetAlignmentSeq))
        # logger.debug('%s', str(self.targetAlignmentSeq.seq))
        for i in xrange(0, len(v), 3):
            vLabel = v[i]
            vQueryLength = int(v[i + 1])
            vTargetLength = int(v[i + 2])
            # logger.debug('qPos = %d, tPos = %d, vLabel = %s, vql = %d, vtl = %d', qPos, tPos, vLabel, vQueryLength, vTargetLength)
            if vLabel == 'M' or vLabel == 'S':
                qAln = qAln + str(self.queryAlignmentSeq[qPos:(qPos + vQueryLength)].seq)
                tAln = tAln + str(self.targetAlignmentSeq[tPos:(tPos + vTargetLength)].seq)
            elif vLabel == 'G':
                if vQueryLength == 0:
                    if vTargetLength % 3 != 0:
                        raise StandardError('cannot process nucleotide gaps with length not a multiple of 3')
                    qAln = qAln + gapChar * (vTargetLength / 3)
                    tAln = tAln + str(self.targetAlignmentSeq[tPos:(tPos + vTargetLength)].seq)
                elif vTargetLength == 0:
                    qAln = qAln + str(self.queryAlignmentSeq[qPos:(qPos + vQueryLength)].seq)
                    tAln = tAln + gapChar * (vQueryLength * 3)
            elif vLabel == '5' or vLabel == '3' or vLabel == 'I' or vLabel == 'F':
                pass
            qPos = qPos + vQueryLength
            tPos = tPos + vTargetLength
            # logger.debug('%d: qPos = %d, tPos = %d, vLabel = %s, vql = %d, vtl = %d', i, qPos, tPos, vLabel, vQueryLength, vTargetLength)
            # logger.debug('qAln: %s', qAln)
            # logger.debug('tAln: %s', tAln)
            # s = Bio.Seq.Seq(tAln, alphabet=Bio.Alphabet.Gapped(self.targetAlignmentSeq.seq.alphabet))
            # s3 = s[:(len(s) - len(s) % 3)]
            # logger.debug('tA_tr: %s%s', str(translateGapped(s3)), '' if len(s) == len(s3) else '.')
        qAlnSeq = Bio.Seq.Seq(qAln, alphabet=Bio.Alphabet.Gapped(Bio.Alphabet.IUPAC.protein))
        tAlnSeq = Bio.Seq.Seq(tAln, alphabet=Bio.Alphabet.Gapped(self.targetAlignmentSeq.seq.alphabet))
        # FIXME: assuming standard translation table -- check whether exonerate supports setting table?
        tAlnProt = translateGapped(tAlnSeq)
        return Bio.Align.MultipleSeqAlignment([Bio.SeqRecord.SeqRecord(qAlnSeq, id=self.queryAlignmentSeq.id), Bio.SeqRecord.SeqRecord(tAlnProt, id='%s_pep' % self.targetAlignmentSeq.id)])

    def __str__(self):
        """String representation of this C{ExonerateResult} instance.

@return: string representation
@rtype: C{str}
"""
        return '%s, query: %s %s(%d -> %d), target: %s:%s %s(%d -> %d), target CDS %dnt' % (self.exonerateModel, self.queryId, self.queryStrand, self.queryAlignmentStart, self.queryAlignmentEnd, self.targetFname, self.targetId, self.targetStrand, self.targetAlignmentStart, self.targetAlignmentEnd, len(self.targetCdsSeq))

    def isEmpty(self):
        """Check whether this result has actually been populated with material from running C{exonerate}.

The the C{query} and C{targetFname} attributes may be populated and
such an empty instance can be used to indicate that scanning the
target file with the query resulted in no hits. This is useful for
recording the fact that no hits were found in a CSV file written with
a L{ExonerateCsvDictWriter}.

@return: C{True} if this result is populated
@rtype: C{bool}

        """
        return self.queryId is None

    def writeAllSeqs(self, f):
        """Write all C{SeqRecord} type attributes to a FASTA file.

@param f: the file to write the sequences to
@type f: C{file} (or file-like that is suitable for C{Bio.SeqIO.write})
"""
        seqList = []
        if self.querySeq is not None:
            seqList.append(self.querySeq)
        if self.queryAlignmentSeq is not None:
            seqList.append(self.queryAlignmentSeq)
        if self.queryCdsSeq is not None:
            seqList.append(self.queryCdsSeq)
        if self.targetAlignmentSeq is not None:
            seqList.append(self.targetAlignmentSeq)
        if self.targetCdsSeq is not None:
            seqList.append(self.targetCdsSeq)
        Bio.SeqIO.write(seqList, f, 'fasta')

    def isQueryReversed(self):
        """Determine whether the query is in reverse orientation.
"""
        return self.queryAlignmentStart > self.queryAlignmentEnd

    def isTargetReversed(self):
        """Determine whether the target is in reverse orientation.
"""
        return self.targetAlignmentStart > self.targetAlignmentEnd


class ExonerateRunner(object):
    """Run C{exonerate} and construct L{ExonerateResult} instances on that basis.
"""

    labelledLineRe = re.compile('([A-Za-z][A-Za-z0-9_]*): (.*)')
    seqStartRe = re.compile('seqStart (.*)')

    def __init__(self):
        # FIXME: consider setting this from environment?
        self.scratchDir = '/tmp'
        self.scratchPrefix = 'exoneratequery'
        pass

    def nextLine(self, f):
        line = f.readline()
        # logger.debug('%s', line.strip())
        if line == '':
            raise StandardError('unexpected EOF')
        if line[-1] == '\n':
            line = line[:-1]
        return line

    def parseString(self, f, label):
        line = self.nextLine(f)
        m = self.labelledLineRe.match(line)
        if m is None:
            raise StandardError('malformed line (expected label %s): %s' % (label, line.strip()))
        if m.group(1) != label:
            raise StandardError('expected label %s but got %s' % (label, m.group(1)))
        return m.group(2).strip()

    def parseInt(self, f, label):
        s = self.parseString(f, label)
        if s == 'NA':
            return None
        return int(s)

    def parseFloat(self, f, label):
        s = self.parseString(f, label)
        if s == 'NA':
            return None
        return float(s)

    def parseSeq(self, f, label, alphabet, seqId):
        line = self.nextLine(f)
        m = self.seqStartRe.match(line)
        if m is None:
            raise StandardError('malformed line (expected seqStart): %s' % line.strip())
        if m.group(1) != label:
            raise StandardError('expected sequence label %s but got %s' % (label, m.group(1)))
        seq = ''
        s = self.nextLine(f)
        while s != 'seqEnd':
            seq = seq + s
            s = self.nextLine(f)
        return Bio.SeqRecord.SeqRecord(Bio.Seq.Seq(seq, alphabet), id=seqId)

    def makeSeqId(self, exonerateResult, seqType):
        return('%s_%s_%s' % (exonerateResult.queryId, exonerateResult.targetId, seqType))

    def parseExonerateResult(self, f, exonerateResult, targetSeqDict):
        line = f.readline()
        if line == '':
            return None
        if line.strip() != 'ryoStart':
            raise StandardError('malformed input: ryoStart missing, got %s instead' % line.strip())
        exonerateResult.exonerateModel = self.parseString(f, 'exonerateModel')
        exonerateResult.queryId = self.parseString(f, 'queryId')
        if exonerateResult.queryId != exonerateResult.querySeq.id:
            raise StandardError('result incompatible with query: querySeq.id = %s, exonerate queryId = %s' % (self.querySeq.id, exonerateResult.queryId))
        exonerateResult.queryDef = self.parseString(f, 'queryDef')
        exonerateResult.queryStrand = self.parseString(f, 'queryStrand')
        exonerateResult.queryAlignmentStart = self.parseInt(f, 'queryAlignmentStart')
        exonerateResult.queryAlignmentEnd = self.parseInt(f, 'queryAlignmentEnd')
        exonerateResult.queryAlignmentLength = self.parseInt(f, 'queryAlignmentLength')
        exonerateResult.queryCdsStart = self.parseInt(f, 'queryCdsStart')
        exonerateResult.queryCdsEnd = self.parseInt(f, 'queryCdsEnd')
        exonerateResult.queryCdsLength = self.parseInt(f, 'queryCdsLength')
        exonerateResult.targetId = self.parseString(f, 'targetId')
        if targetSeqDict is not None:
            if exonerateResult.targetId not in targetSeqDict:
                raise StandardError, 'found targetId %s but no corresponding sequence' % exonerateResult.targetId
            exonerateResult.targetSeq = targetSeqDict[exonerateResult.targetId]
        exonerateResult.targetDef = self.parseString(f, 'targetDef')
        exonerateResult.targetStrand = self.parseString(f, 'targetStrand')
        exonerateResult.targetAlignmentStart = self.parseInt(f, 'targetAlignmentStart')
        exonerateResult.targetAlignmentEnd = self.parseInt(f, 'targetAlignmentEnd')
        exonerateResult.targetAlignmentLength = self.parseInt(f, 'targetAlignmentLength')
        exonerateResult.targetCdsStart = self.parseInt(f, 'targetCdsStart')
        exonerateResult.targetCdsEnd = self.parseInt(f, 'targetCdsEnd')
        exonerateResult.targetCdsLength = self.parseInt(f, 'targetCdsLength')
        exonerateResult.rawScore = self.parseInt(f, 'rawScore')
        exonerateResult.percentIdentity = self.parseFloat(f, 'percentIdentity')
        exonerateResult.percentSimilarity = self.parseFloat(f, 'percentSimilarity')
        exonerateResult.equivalencedTotal = self.parseInt(f, 'equivalencedTotal')
        exonerateResult.equivalencedIdentity = self.parseInt(f, 'equivalencedIdentity')
        exonerateResult.equivalencedSimilarity = self.parseInt(f, 'equivalencedSimilarity')
        exonerateResult.equivalencedMismatches = self.parseInt(f, 'equivalencedMismatches')
        exonerateResult.vulgar = self.parseString(f, 'vulgar')
        exonerateResult.queryCdsSeq = self.parseSeq(f, 'queryCds', Bio.Alphabet.IUPAC.ambiguous_dna, self.makeSeqId(exonerateResult, 'qcds'))
        exonerateResult.queryAlignmentSeq = self.parseSeq(f, 'queryAlignment', Bio.Alphabet.IUPAC.protein, self.makeSeqId(exonerateResult, 'qaln'))
        exonerateResult.targetCdsSeq = self.parseSeq(f, 'targetCds', Bio.Alphabet.IUPAC.ambiguous_dna, self.makeSeqId(exonerateResult, 'tcds'))
        exonerateResult.targetAlignmentSeq = self.parseSeq(f, 'targetAlignment', Bio.Alphabet.IUPAC.ambiguous_dna, self.makeSeqId(exonerateResult, 'taln'))
        if exonerateResult.exonerateModel == 'protein2genome:local':
            # FIXME: kludge -- throwing away rubbish output of exonerate, would be much better not to generate it in the first place
            # FIXME: would seem better to not include fields that don't apply and that exonerate populates with garbage
            exonerateResult.queryCdsSeq = None
        elif exonerateResult.exonerateModel == 'affine:local:dna2dna':
            exonerateResult.targetCdsStart = None
            exonerateResult.targetCdsEnd = None
            exonerateResult.targetCdsLength = None
            exonerateResult.queryCdsSeq = None
            exonerateResult.targetCdsSeq = None
        else:
            raise StandardError('unsupported exonerate model: %s' % exonerateResult.exonerateModel)
        line = self.nextLine(f)
        if line.strip() != 'ryoEnd':
            raise StandardError('malformed input: ryoEnd missing')
        return exonerateResult

    def parse(self, querySeq, targetFname, exonerateModel, bestn=None, minPercentIdentity=None, addRawTargetSeqs=False):
        """Run C{exonerate} and return a C{list} of C{ExonerateResult}s.

@param querySeq: the query sequence
@type querySeq: C{Bio.SeqRecord.SeqRecord}
@param targetFname: the name of the FASTA sequence file containing the targets
@type targetFname: C{str}
@param exonerateModel: the alignment model
@type exonerateModel: C{str}
@param bestn: max. number of hits
@type bestn: C{int}, or C{None} to use default
@param minPercentIdentity: minimum percent identity threshold for reporting alignments
@type minPercentIdentity: C{float}, or C{None} to use default
@param addRawTargetSeqs: whether to add C{SeqRecord}s of target sequences to L{ExonerateResult} instances
@type addRawTargetSeqs: C{bool}
@return: list of results
@rtype: C{list} of L{ExonerateResult}
"""
        # FIXME: hard-coded to generate query scratch file in cwd
        targetSeqDict = None
        if addRawTargetSeqs:
            targetSeqDict = Bio.SeqIO.to_dict(Bio.SeqIO.parse(targetFname, 'fasta'))
        queryScratchFd, queryScratchFname = tempfile.mkstemp('.fasta', self.scratchPrefix, self.scratchDir)
        try:
            queryScratchFile = os.fdopen(queryScratchFd, 'w')
            Bio.SeqIO.write(querySeq, queryScratchFile, 'fasta')
            queryScratchFile.close()
            exonerateArgv = ['exonerate', '--model', exonerateModel, '--verbose', '0', '--showalignment', 'no', '--showvulgar', 'no']
            if bestn is not None:
                exonerateArgv.extend(['--bestn', '%d' % bestn])
            if minPercentIdentity is not None:
                exonerateArgv.extend(['--percent', '%f' % minPercentIdentity])
            exonerateArgv.extend(['--ryo', 'ryoStart\\nexonerateModel: %m\\nqueryId: %qi\\nqueryDef: %qd\\nqueryStrand: %qS\\nqueryAlignmentStart: %qab\\nqueryAlignmentEnd: %qae\\nqueryAlignmentLength: %qal\\nqueryCdsStart: NA\\nqueryCdsEnd: NA\\nqueryCdsLength: NA\\ntargetId: %ti\\ntargetDef: %td\\ntargetStrand: %tS\\ntargetAlignmentStart: %tab\\ntargetAlignmentEnd: %tae\\ntargetAlignmentLength: %tal\\ntargetCdsStart: %tcb\\ntargetCdsEnd: %tce\\ntargetCdsLength: %tcl\\nrawScore: %s\\npercentIdentity: %pi\\npercentSimilarity: %ps\\nequivalencedTotal: %et\\nequivalencedIdentity: %ei\\nequivalencedSimilarity: %es\\nequivalencedMismatches: %em\\nvulgar: %V\\nseqStart queryCds\\n%qcsseqEnd\\nseqStart queryAlignment\\n%qasseqEnd\\nseqStart targetCds\\n%tcsseqEnd\\nseqStart targetAlignment\\n%tasseqEnd\\nryoEnd\\n', queryScratchFname, targetFname])
            logger.debug('%s', ' '.join(exonerateArgv))
            p = subprocess.Popen(exonerateArgv, stdin=subprocess.PIPE, stdout=subprocess.PIPE)
            pid = os.fork()
            if pid == 0:
                p.stdout.close()
                # exonerate seems to be one of those stupid programs that don't read from stdin...
                # p.stdin.write('%s\n' % querySr.format('fasta'))
                p.stdin.close()
                os._exit(0)
            p.stdin.close()
            exonerateResultList = []
            exonerateResult = self.parseExonerateResult(p.stdout, ExonerateResult(querySeq, targetFname), targetSeqDict)
            while exonerateResult is not None:
                exonerateResultList.append(exonerateResult)
                exonerateResult = self.parseExonerateResult(p.stdout, ExonerateResult(querySeq, targetFname), targetSeqDict)
            p.stdout.close()
            wPid, wExit = os.waitpid(pid, 0)
            if pid != wPid:
                raise StandardError('wait returned pid %s (expected %d)' % (wPid, pid))
            if wExit != 0:
                raise StandardError('wait on forked process returned %d' % wExit)
            r = p.wait()
            if r != 0:
                raise StandardError('exonerate process exited with %d' % r)
        finally:
            if paftol.keepTmp:
                logger.warning('not deleting query scratch file %s', queryScratchFname)
            else:
                os.unlink(queryScratchFname)
        return exonerateResultList


class ExonerateCsvDictWriter(object):
    """Write CSV files containing stats from L{ExonerateResult} instances.
"""

    csvFieldnames = [
        'querySeqId', 'targetFname', 'genomeName', 'querySeqLength',
        'queryId',
        'queryAlignmentStart', 'queryAlignmentEnd', 'queryAlignmentLength',
        'queryCdsStart', 'queryCdsEnd', 'queryCdsLength', 'queryStrand',
        'targetId',
        'targetAlignmentStart', 'targetAlignmentEnd', 'targetAlignmentLength',
        'targetCdsStart', 'targetCdsEnd', 'targetCdsLength', 'targetStrand',
        'rawScore','percentIdentity','percentSimilarity',
        'equivalencedTotal', 'equivalencedIdentity', 'equivalencedSimilarity', 'equivalencedMismatches',
        'targetCdsSeqLength']

    def __init__(self, csvfile):
        """Create an L{ExonerateCsvDictWriter}.

The C{csvfile} parameter may either be a string containing the name of
the CSV file to be written, or a file like object. If it's a string,
the file will be opened for writing, thus erasing any previous
content.

@param csvfile: The file to write to
@type csvfile: C{str} or file like
"""
        if isinstance(csvfile, types.StringType):
            self.csvfile = open(csvfile, 'w')
            self.csvDictWriter = csv.DictWriter(self.csvfile, self.csvFieldnames)
        else:
            self.csvfile = None
            self.csvDictWriter = csv.DictWriter(csvfile, csvFieldnames)
        self.csvDictWriter.writeheader()

    def close(self):
        """Close this writer.

If the underlying file was opened at construction of this instance, it
is closed. If a file like object was passed in at construction, it's
the caller's responsibility to close it.
"""
        if self.csvfile is not None:
            self.csvfile.close()
        self.csvfile = None
        self.csvDictWriter = None

    def writeCsvStats(self, exonerateResult):
        """Write a CSV row containing data from an C{ExonerateResult} instance.

@param exonerateResult: the instance from which to take the row's content
"""
        if self.csvDictWriter is None:
            raise StandardError('illegal state: no DictWriter (close called previously?)')
        d = {}
        d['querySeqId'] = None if exonerateResult.querySeq is None else exonerateResult.querySeq.id
        d['targetFname'] = exonerateResult.targetFname
        d['querySeqLength'] = None if exonerateResult.querySeq is None else len(exonerateResult.querySeq)
        d['queryId'] = exonerateResult.queryId
        d['queryAlignmentStart'] = exonerateResult.queryAlignmentStart
        d['queryAlignmentEnd'] = exonerateResult.queryAlignmentEnd
        d['queryAlignmentLength'] = exonerateResult.queryAlignmentLength
        d['queryCdsStart'] = exonerateResult.queryCdsStart
        d['queryCdsEnd'] = exonerateResult.queryCdsEnd
        d['queryCdsLength'] = exonerateResult.queryCdsLength
        d['queryStrand'] = exonerateResult.queryStrand
        d['targetId'] = exonerateResult.targetId
        d['targetAlignmentStart'] = exonerateResult.targetAlignmentStart
        d['targetAlignmentEnd'] = exonerateResult.targetAlignmentEnd
        d['targetAlignmentLength'] = exonerateResult.targetAlignmentLength
        d['targetCdsStart'] = exonerateResult.targetCdsStart
        d['targetCdsEnd'] = exonerateResult.targetCdsEnd
        d['targetCdsLength'] = exonerateResult.targetCdsLength
        d['targetStrand'] = exonerateResult.targetStrand
        d['rawScore'] = exonerateResult.rawScore
        d['percentIdentity'] = exonerateResult.percentIdentity
        d['percentSimilarity'] = exonerateResult.percentSimilarity
        d['equivalencedTotal'] = exonerateResult.equivalencedTotal
        d['equivalencedIdentity'] = exonerateResult.equivalencedIdentity
        d['equivalencedSimilarity'] = exonerateResult.equivalencedSimilarity
        d['equivalencedMismatches'] = exonerateResult.equivalencedMismatches
        d['targetCdsSeqLength'] = None if exonerateResult.targetCdsSeq is None else len(exonerateResult.targetCdsSeq)
        self.csvDictWriter.writerow(d)


def alignMerge(pairwiseAlignmentList):
    gapChar = '-'
    alphabet = pairwiseAlignmentList[0][0].seq.alphabet
    refId  = pairwiseAlignmentList[0][0].id
    refseqList = []
    otherseqList = []
    refSymbolSeq = []
    otherSymbolSeqList = []
    vulgarLetterAnnotationList = []
    staralignLetterAnnotationList = []
    for i in xrange(len(pairwiseAlignmentList)):
        if 'vulgar' in pairwiseAlignmentList[i][1].letter_annotations:
            vulgarLetterAnnotationList.append([])
        else:
            vulgarLetterAnnotationList.append(None)
        staralignLetterAnnotationList.append([])
    for pairwiseAlignment in pairwiseAlignmentList:
        refseqList.append(str(pairwiseAlignment[0].seq))
        otherseqList.append(str(pairwiseAlignment[1].seq))
        otherSymbolSeqList.append([])
    i = [0] * len(pairwiseAlignmentList)
    while max([len(refseqList[n]) - i[n] for n in xrange(len(pairwiseAlignmentList))]) > 0:
        refGap = False
        refSym = None
        refSymAlignIndex = None
        for n in xrange(len(pairwiseAlignmentList)):
            # sys.stderr.write('n = %d, i[n] = %d, l = %d\n' % (n, i[n], len(refseqList[n])))
            if i[n] == len(refseqList[n]):
                r  = gapChar
                # Bio.AlignIO.write(pairwiseAlignmentList[n], sys.stderr, 'fasta')
            else:
                r = refseqList[n][i[n]].lower()
            if r == gapChar:
                refGap = True
            elif refSym is None:
                refSym = r
                refSymAlignIndex = n
            else:
                if r != refSym:
                    raise StandardError, 'reference sequences inconsistent, found symbols %s (%s:%d) and %s (%s:%d)' % (refSym, pairwiseAlignmentList[refSymAlignIndex][1].id, refSymAlignIndex, r, pairwiseAlignmentList[n][1].id, n)
        if refGap:
            refSymbolSeq.append(gapChar)
            for n in xrange(len(pairwiseAlignmentList)):
                if i[n] ==len(otherseqList[n]):
                    otherSymbolSeqList[n].append(gapChar)
                    vulgarLetterAnnotationList[n].append(None)
                    staralignLetterAnnotationList[n].append(None)
                else:
                    if refseqList[n][i[n]] == gapChar:
                        otherSymbolSeqList[n].append(otherseqList[n][i[n]])
                        if 'vulgar' in pairwiseAlignmentList[n][1].letter_annotations:
                            vulgarLetterAnnotationList[n].append(pairwiseAlignmentList[n][1].letter_annotations['vulgar'][i[n]])
                        staralignLetterAnnotationList[n].append(i[n])
                        i[n] = i[n] + 1
                    elif refseqList[n][i[n]].lower() == refSym:
                        otherSymbolSeqList[n].append(gapChar)
                        if 'vulgar' in pairwiseAlignmentList[n][1].letter_annotations:
                            vulgarLetterAnnotationList[n].append(None)
                        staralignLetterAnnotationList[n].append(None)
                    else:
                        raise StandardError, 'weird: refSym = %s, refseqsym = %s, n = %d, i[n] = %d, id = %s' % (refSym, refseqList[n][i[n]], n, i[n], pairwiseAlignmentList[n][1].id)
        else:
            refSymbolSeq.append(refSym)
            for n in xrange(len(pairwiseAlignmentList)):
                otherSymbolSeqList[n].append(otherseqList[n][i[n]])
                if 'vulgar' in pairwiseAlignmentList[n][1].letter_annotations:
                    vulgarLetterAnnotationList[n].append(pairwiseAlignmentList[n][1].letter_annotations['vulgar'][i[n]])
                staralignLetterAnnotationList[n].append(i[n])
                i[n] = i[n] + 1
    refSr = Bio.SeqRecord.SeqRecord(Bio.Seq.Seq(''.join(refSymbolSeq), alphabet=alphabet), id='%s_ref' % refId)
    srList = [refSr]
    # sys.stderr.write('ref. %s: length: %d\n' % (refSr.id, len(refSr)))
    for n in xrange(len(pairwiseAlignmentList)):
        letter_annotations = {}
        if vulgarLetterAnnotationList[n] is not None:
            letter_annotations['vulgar'] = vulgarLetterAnnotationList[n]
        letter_annotations['staralign'] = staralignLetterAnnotationList[n]
        sr = Bio.SeqRecord.SeqRecord(Bio.Seq.Seq(''.join(otherSymbolSeqList[n]), alphabet=alphabet), id='%s_%05d' % (pairwiseAlignmentList[n][1].id, n), letter_annotations = letter_annotations)
        # sys.stderr.write('%s: length %d\n' % (sr.id, len(sr)))
        srList.append(sr)
    # Bio.SeqIO.write(srList, sys.stderr, 'fasta')
    return Bio.Align.MultipleSeqAlignment(srList)


class ExonerateStarAlignment(object):

    def __init__(self, reference, fastaFname):
        self.reference = reference
        self.fastaFname = fastaFname
        self.xstarAlignment = None
        self.makeStarAlignment()

    def makeExonerateResult(self, seqRecord, exonerateRunner):
        tmpFastaFd, tmpFastaFname = tempfile.mkstemp('.fasta', 'xstar', '.')
        try:
            with os.fdopen(tmpFastaFd, 'w') as tmpFastaFile:
                Bio.SeqIO.write([seqRecord], tmpFastaFile, 'fasta')
            # only best alignment for now, extending to multiple local alignments will require name fiddling...
            exonerateResultList = exonerateRunner.parse(self.reference, tmpFastaFname, 'affine:local', bestn=1, addRawTargetSeqs=True)
        except StandardError as e:
            raise e
        finally:
            # FIXME: check keepTmp when moving this into pypaftol
            os.unlink(tmpFastaFname)
            pass
        return exonerateResultList

    def makeStarAlignment(self):
        exonerateRunner = paftol.tools.ExonerateRunner()
        exonerateResultList = []
        for seqRecord in Bio.SeqIO.parse(self.fastaFname, 'fasta'):
            exonerateResultList.extend(self.makeExonerateResult(seqRecord, exonerateRunner))
        # sys.stderr.write('got %d exonerate results\n')
        exonerateResultList.sort(lambda e1, e2: cmp(e1.queryAlignmentStart, e2.queryAlignmentStart))
        pairwiseAlignmentList = [er.nucleotideAlignment(appendFlanking=True) for er in exonerateResultList]
        self.xstarAlignment = alignMerge(pairwiseAlignmentList)

    def epsSketchSr(self, sr, epsFile, y, symbolWidth, symbolHeight, symbolRgb, gapRgb, unalignedRgb):
        gapChar = '-'
        s = str(sr.seq)
        v = None
        if 'vulgar' in sr.letter_annotations:
            v = sr.letter_annotations['vulgar']
        rgbList = [symbolRgb, gapRgb, unalignedRgb]
        lastRgbIndex = None
        rgbIndex = None
        lastRgbStart = None
        epsFile.write('%% sequence %s\n' % sr.id)
        for i in xrange(len(sr)):
            if s[i] == gapChar:
                if v is not None and v[i] is None:
                    rgbIndex = 2
                else:
                    rgbIndex = 1
            else:
                if v is not None and v[i] is None:
                    rgbIndex = 2
                else:
                    rgbIndex = 0
            if rgbIndex != lastRgbIndex or i == len(sr) - 1:
                if lastRgbIndex is not None:
                    epsFile.write('%f %f %f setrgbcolor\n' % rgbList[lastRgbIndex])
                    epsFile.write('%f %f moveto %f 0 rlineto stroke\n' % (lastRgbStart * symbolWidth , y, (i - lastRgbStart) * symbolWidth))
                lastRgbIndex = rgbIndex
                lastRgbStart = i

    def epsSketch(self, epsFile, width=None, height=None):
        if self.xstarAlignment is None:
            raise StandardError, 'no star alignment to sketch'
        if width is None:
            width = self.xstarAlignment.get_alignment_length()
        if height is None:
            height = len(self.xstarAlignment)
        symbolWidth = width / self.xstarAlignment.get_alignment_length()
        lineHeight = height / len(self.xstarAlignment)
        symbolHeight = lineHeight * 0.8
        y = height - symbolHeight
        epsFile.write('%!PS-Adobe-3.0 EPSF-3.0\n')
        epsFile.write('%%%%BoundingBox: 0 0 %f %f\n' % (width, height))
        epsFile.write('%%Creator: exonerateStarAlignment.epsSketch\n')
        epsFile.write('%%Title: star alignment sketch\n')
        epsFile.write('%%EndComments\n')
        epsFile.write('%s setlinewidth 0 setlinecap\n' % symbolHeight)
        self.epsSketchSr(self.xstarAlignment[0], epsFile, height - 0.5 * symbolHeight, symbolWidth, symbolHeight, (0.0, 0.0, 1.0, ), (0.5, 0.5, 1.0, ), (1.0, 0.0, 0.0, ))
        y = lineHeight * (len(self.xstarAlignment) - 1.5)
        for sr in self.xstarAlignment[1:]:
            self.epsSketchSr(sr, epsFile, y, symbolWidth, symbolHeight, (0.0, 0.0, 0.0, ), (0.5, 0.5, 0.5, ), (0.9, 0.9, 0.9, ))
            y = y - lineHeight
        epsFile.write('%%EOF\n')


class BwaRunner(object):
    """Wrapper class for running C{bwa}.

Instance variables are passed to C{bwa} via its CLI. C{None} means no
command line option / parameter to be specified, which usually results
in C{bwa} using a default (which may depend on the C{bwa} version).

@ivar numThreads: BWA number of threads (C{-t} option)
@type numThreads: C{int}, or C{None}
@ivar minSeedLength: BWA minimum seed length (C{-k} option)
@type minSeedLength: C{int}, or C{None}
@ivar scoreThreshold: BWA score threshold for recording reads as mapped (C{-T} option)
@type scoreThreshold: C{int}, or C{None}
@ivar reseedTrigger: BWA re-seed trigger (C{-r} option)
@type reseedTrigger: C{float}, or C{None}
"""
    def __init__(self, numThreads=None, minSeedLength=None, scoreThreshold=None, reseedTrigger=None, workingDirectory=None):
        """Constructor.

Parameters correspond to instance variables, see their documentation.
"""
        self.numThreads = numThreads
        self.minSeedLength = minSeedLength
        self.scoreThreshold = scoreThreshold
        self.reseedTrigger = reseedTrigger
        self.workingDirectory = workingDirectory

    def indexReferenceArgv(self, referenceFname):
        return ['bwa', 'index', referenceFname]

    def mappingMemArgv(self, referenceFname, forwardReadsFname, reverseReadsFname=None):
        argv = ['bwa', 'mem', '-M']
        if self.minSeedLength is not None:
            argv.extend(['-k', '%d' % self.minSeedLength])
        if self.reseedTrigger is not None:
            argv.extend(['-r', '%f' % self.reseedTrigger])
        if self.scoreThreshold is not None:
            argv.extend(['-T', '%d' % self.scoreThreshold])
        if self.numThreads is not None:
            argv.extend(['-t', '%d' % self.numThreads])
        argv.append(referenceFname)
        argv.append(forwardReadsFname)
        if reverseReadsFname is not None:
            argv.append(reverseReadsFname)
        return argv

    def indexReference(self, referenceFname):
        """Index a reference sequence (using C{bwa index}).

Notice that the index files created by BWA will be generated as a side
effect. It is the responsibility of clients to tidy these up, if
necessary.

@param referenceFname: the name of the reference sequence FASTA file which is to be indexed
@type referenceFname: C{str}

        """
        bwaIndexArgv = self.indexReferenceArgv(referenceFname)
        logger.debug('%s', ' '.join(bwaIndexArgv))
        subprocess.check_call(bwaIndexArgv)

    def processBwa(self, samAlignmentProcessor, referenceFname, forwardReadsFname, reverseReadsFname=None):
        """Process reads mapped to to reference sequences.

This method runs C{bwa} with the reference sequence, forwards reads
file and reverse reads file as specified by the respective parameters.
The C{samAlignmentProcessor} must implement a C{processSamAlignment}
method which is called for each SAM alignment emitted by C{bwa}.

In other words, by implementing a C{processSamAlignment} method, a
class is of a suitable "duck type" to be used as a
C{samAlignmentProcessor}.

@param samAlignmentProcessor: the object for processing SAM alignments
@type samAlignmentProcessor: object of suitable "duck type"
@param referenceFname: name of the reference sequence file (FASTA format)
@type referenceFname: C{str}
@param forwardReadsFname: name of the forward reads file (FASTQ format)
@type forwardReadsFname: C{str}
@param reverseReadsFname: name of the reverse reads file (FASTQ format)
@type reverseReadsFname: C{str}, or C{None}
"""
        sys.stderr.write('effective mapReadsBwa logging level: %d\n' % logger.getEffectiveLevel())
        logger.debug('mapping reads to gene sequences')
        self.indexReference(referenceFname)
        bwaArgv = self.mappingMemArgv(referenceFname, forwardReadsFname, reverseReadsFname)
        logger.debug('%s', ' '.join(bwaArgv))
        bwaProcess = subprocess.Popen(bwaArgv, stdout=subprocess.PIPE, cwd=self.workingDirectory)
        # samtoolsArgv = ['samtools', 'view', '-h', '-S', '-F', '4', '-']
        # logger.debug('%s', ' '.join(samtoolsArgv))
        # samtoolsProcess = subprocess.Popen(samtoolsArgv, stdin=bwaProcess.stdout.fileno(), stdout=subprocess.PIPE, cwd = self.workingDirectory)
        samLine = bwaProcess.stdout.readline()
        while samLine != '':
            # logger.debug(samLine)
            if samLine[0] != '@':
                samAlignment = SamAlignment(samLine)
                samAlignmentProcessor.processSamAlignment(samAlignment)
            samLine = bwaProcess.stdout.readline()
        bwaProcess.stdout.close()
        # samtoolsProcess.stdout.close()
        bwaReturncode = bwaProcess.wait()
        # samtoolsReturncode = samtoolsProcess.wait()
        if bwaReturncode != 0:
            raise StandardError('process "%s" returned %d' % (' '.join(bwaArgv), bwaReturncode))
        # if samtoolsReturncode != 0:
        #     raise StandardError('process "%s" returned %d' % (' '.join(samtoolsArgv), samtoolsReturncode))


class BlastAlignmentProcessor(object):
    """Simple BLAST alignment processor.

This class provides a dictionary in its C{blastAlignmentDict} instance
variable. Keys are query names and values are lists of
L{BlastAlignment} instances. The C{processBlastAlignment} method
appends BLAST alignment to the appropriate list, adding an empty list
to the C{blastAlignmentDict} when necessary.

@ivar blastAlignmentDict: dictionary of lists of BLAST alignments
@type blastAlignmentDict: C{dict}
"""
    def __init__(self):
        """Constructor.
"""
        self.blastAlignmentDict = {}

    def processBlastAlignment(self, query, blastAlignment):
        """Process a BLAST alignment.
"""
        if query not in self.blastAlignmentDict:
            self.blastAlignmentDict[query] = []
        self.blastAlignmentDict[query].append(blastAlignment)


class BlastRunner(object):
    """Wrapper class for running BLAST programs.

This is a base class for runners that wrap specific BLAST programs
(C{blastn}, C{tblastn} etc.), and that should therefore be considered
abstract, i.e. it should not be instantiated.

Instance variables correspond to command line options that are the
same for all BLAST programs. Setting an instance variable to C{None}
means no command line options are generated, usually resulting in a
default to be used.

@ivar numThreads: number of threads
@type numThreads: C{int}, or C{None}
@ivar gapOpen: gap open penalty
@type gapOpen: C{int}, or C{None}
@ivar gapExtend: gap extension penalty
@type gapExtend: C{int}, or C{None}
@ivar maxTargetSeqs: maximum number of target sequences
@type maxTargetSeqs: C{int}, or C{None}
@ivar numAlignments: number of sequences to show alignments for
@type numAlignments: C{int}, or C{None}
@ivar maxHsps: maximum number of HSPs per subject sequence
@type maxHsps: C{int}, or C{None}
@ivar evalue: E-value threshold for reporting hits
@type evalue: C{float}, or C{None}
@ivar windowSize: multiple hits window size
@type windowSize: C{int}, or C{None}
"""
    def __init__(self, numThreads, gapOpen, gapExtend, maxTargetSeqs, numAlignments, maxHsps, evalue, windowSize):
        self.numThreads = numThreads
        self.gapOpen = gapOpen
        self.gapExtend = gapExtend
        self.maxTargetSeqs = maxTargetSeqs
        self.numAlignments = numAlignments
        self.maxHsps = maxHsps
        self.evalue = evalue
        self.windowSize = windowSize

    def indexDatabase(self, databaseFname, dbtype):
        makeblastdbArgv = ['makeblastdb', '-dbtype', dbtype, '-in', databaseFname, '-parse_seqids']
        logger.debug('%s', ' '.join(makeblastdbArgv))
        makeblastdbProcess = subprocess.check_call(makeblastdbArgv)

    def makeBlastArgv(self, blastProgram, databaseFname):
        blastArgv = [blastProgram]
        if self.numThreads is not None:
            blastArgv.extend(['-num_threads', '%d' % self.numThreads])
        if self.gapOpen is not None:
            blastArgv.extend(['-gapopen', '%d' % self.gapOpen])
        if self.gapExtend is not None:
            blastArgv.extend(['-gapextend', '%d' % self.gapExtend])
        if self.maxTargetSeqs is not None:
            blastArgv.extend(['-max_target_seqs', '%d' % self.maxTargetSeqs])
        if self.numAlignments is not None:
            blastArgv.extend(['-num_alignments', '%d' % self.numAlignments])
        if self.maxHsps is not None:
            blastArgv.extend(['-max_hsps', '%d' % self.maxHsps])
        if self.evalue is not None:
            blastArgv.extend(['-evalue', '%1.12g' % self.evalue])
        if self.windowSize is not None:
            blastArgv.extend(['-window_size', '%d' % self.windowSize])
        blastArgv.extend(['-db', databaseFname, '-outfmt', '5'])
        return blastArgv

    def processBlast(self, blastProgram, blastAlignmentProcessor, databaseFname, queryList):
        """Run a BLAST program and process the alignments output by it.

The C{blastAlignmentProcessor} object needs to provide a
C{processBlastAlignment} method that takes the name of the query and a
corresponding BLAST alignment as a parameter.

In other words, by implementing a C{processBlastAlignment} method, a
class is of a suitable "duck type" to be used as a
C{blastAlignmentProcessor}.

@param blastProgram: the BLAST program to be used, provided by the subclass
@type blastProgram: C{str}
@param blastAlignmentProcessor: the BLAST alignment processor
@type blastAlignmentProcessor: suitable duck type
@param databaseFname: the database to be used, needs to be suitable indexed by C{makeblastdb}
@type databaseFname: C{str}
@param queryList: List of query sequences
@type queryList: C{list} of C{Bio.SeqRecordSeqRecord} instances
"""
        blastArgv = self.makeBlastArgv(blastProgram, databaseFname)
        logger.debug('%s', ' '.join(blastArgv))
        blastProcess = subprocess.Popen(blastArgv, stdin=subprocess.PIPE, stdout=subprocess.PIPE)
        pid = os.fork()
        if pid == 0:
            blastProcess.stdout.close()
            for query in queryList:
                blastProcess.stdin.write(query.format('fasta'))
            blastProcess.stdin.close()
            os._exit(0)
        blastProcess.stdin.close()
        for blastRecord in Bio.Blast.NCBIXML.parse(blastProcess.stdout):
            for alignment in blastRecord.alignments:
                blastAlignmentProcessor.processBlastAlignment(blastRecord.query, alignment)
        blastProcess.stdout.close()
        wPid, wExit = os.waitpid(pid, 0)
        if pid != wPid:
            raise StandardError, 'wait returned pid %s (expected %d)' % (wPid, pid)
        if wExit != 0:
            raise StandardError, 'wait on forked process returned %d' % wExit
        blastReturncode = blastProcess.wait()
        if blastReturncode != 0:
            raise StandardError, 'process "%s" returned %d' % (' '.join(blastArgv), blastReturncode)


class BlastnRunner(BlastRunner):
    """Runner for C{blastn}.
"""

    def __init__(self, numThreads=None, gapOpen=None, gapExtend=None, maxTargetSeqs=None, numAlignments=None, maxHsps=None, evalue=None, windowSize=None):
        super(BlastnRunner, self).__init__(numThreads, gapOpen, gapExtend, maxTargetSeqs, numAlignments, maxHsps, evalue, windowSize)

    def indexDatabase(self, databaseFname):
        super(BlastnRunner, self).indexDatabase(databaseFname, 'nucl')

    def processBlast(self, blastAlignmentProcessor, databaseFname, queryList):
        logger.debug('BlastnRunner.processBlast')
        super(BlastnRunner, self).processBlast('blastn', blastAlignmentProcessor, databaseFname, queryList)


class TblastnRunner(BlastRunner):
    """Runner for C{tblastn}.
"""

    def __init__(self, numThreads=None, gapOpen=None, gapExtend=None, maxTargetSeqs=None, numAlignments=None, maxHsps=None, evalue=None, windowSize=None):
        super(TblastnRunner, self).__init__(numThreads, gapOpen, gapExtend, maxTargetSeqs, numAlignments, maxHsps, evalue, windowSize)

    def indexDatabase(self, databaseFname):
        super(TblastnRunner, self).indexDatabase(databaseFname, 'nucl')

    def processTblastn(self, blastAlignmentProcessor, databaseFname, queryList):
        super(TblastnRunner, self).processBlast('tblastn', blastAlignmentProcessor, databaseFname, queryList)


def selectLongestReads(readsList, numReads):
    """Select the longest reads from a list of reads.
@param readsList: list of reads
@type readsList: C{list} of C{Bio.SeqRecord.SeqRecord}
@param numReads: the number of longest reads to return
@type numReads: C{int}
"""
    # maximally naive implementation would be:
    # return readsList[:numReads]
    l = readsList[:]
    l.sort(cmp=lambda x, y: cmp(len(x), len(y)), reverse=True)
    return l[:numReads]


class SpadesRunner(object):
    """

"""

    SINGLE = 1
    INTERLACED = 2
    PAIRED = 3

    def __init__(self, numThreads=None, covCutoff=None, kvalList=None):
        self.numThreads = numThreads
        self.covCutoff = covCutoff
        if kvalList is None:
            self.kvalList = None
        else:
            self.kvalList = kvalList[:]

    def assemble(self, readsFname, libraryType, outputDirname, workDirname=None):
        if libraryType == self.INTERLACED:
            spadesInputArgs = ['--12', readsFname]
        elif libraryType == self.SINGLE:
            spadesInputArgs = ['-s', readsFname]
        else:
            raise StandardError, 'library type %d unknown / unsupported' % libraryType
        spadesArgv = ['spades.py', '--only-assembler']
        if self.numThreads is not None:
            spadesArgv.extend(['--threads', '%d' % self.numThreads])
        if self.covCutoff is not None:
            spadesArgv.extend(['--cov-cutoff', '%s' % self.covCutoff])
        if self.kvalList is not None:
            spadesArgv.extend(['-k', ','.join(['%d' % k for k in self.spadesKvalList])])
        spadesArgv.extend(spadesInputArgs)
        spadesArgv.extend(['-o', outputDirname])
        logger.debug('%s', ' '.join(spadesArgv))
        spadesProcess = subprocess.Popen(spadesArgv, cwd=workDirname)
        returncode = spadesProcess.wait()
        if returncode != 0:
            # raise StandardError('spades process "%s" exited with %d' % (' '.join(spadesArgv), returncode))
            logger.warning('spades process "%s" exited with %d', ' '.join(spadesArgv), returncode)
        if workDirname is None:
            contigFname = os.path.join(outputDirname, 'contigs.fasta')
        else:
            contigFname = os.path.join(workDirname, outputDirname, 'contigs.fasta')
        # logger.debug('contigFname: %s', spadesContigFname)
        if os.path.exists(contigFname):
            contigList = list(Bio.SeqIO.parse(contigFname, 'fasta'))
            # logger.debug('contigFname: %s, %d contigs', contigFname, len(contigList))
        else:
            contigList = None
            # logger.debug('contigFname: %s, no contigs', contigFname)
        return contigList


class MeanAndStddev(object):

    def __init__(self, l):
        self.l = l
        self.mean = sum(l) / float(len(l))
        sdList = []
        for num in l:
            sdList.append((self.mean - num) ** 2)
        self.std = math.sqrt(sum(sdList) / (len(sdList) - 1))


class DataFrame(object):

    def __init__(self, columnHeaderList):
        self.columnHeaderList = columnHeaderList[:]
        self.rowDictList = []

    def addRow(self, rowDict):
        if set(rowDict.keys()) != set(self.columnHeaderList):
            raise StandardError, 'key set %s is not compatible with column headers %s' % (', '.join([str(k) for k in rowDict.keys()]), ', '.join(self.columnHeaderList))
        self.rowDictList.append(copy.copy(rowDict))

    def nrow(self):
        return len(self.rowDictList)

    def getRowDict(self, rowIndex):
        return self.rowDictList[rowIndex]

    def writeCsv(self, f):
        csvDictWriter = csv.DictWriter(f, self.columnHeaderList)
        csvDictWriter.writeheader()
        for rowDict in self.rowDictList:
            csvDictWriter.writerow(rowDict)

    def getColumn(self, columnName):
        return [rowDict[columnName] for rowDict in self.rowDictList]

    def colMeanAndStddev(self, columnName):
        return MeanAndStddev(self.getColumn(columnName))

    def colStats(self, columnName):
        d = {}
        for columnName in self.columnHeaderList:
            d[columnHeader] = self.colMeanAndStddev(columnName)
        return d


def numIdenticalSymbols(sr1, sr2, ignoreCase=True):
    if len(sr1) != len(sr2):
        raise StandardError, 'sequences %s and %s differ in length: %d != %d' % (sr1.id, sr2.id, len(sr1), len(sr2))
    gapChar = None
    if isinstance(sr1.seq.alphabet, Bio.Alphabet.Gapped) and isinstance(sr2.seq.alphabet, Bio.Alphabet.Gapped):
        if sr1.seq.alphabet.gap_char == sr2.seq.alphabet.gap_char:
            gapChar = sr1.seq.alphabet.gap_char
    s1 = str(sr1.seq)
    s2 = str(sr2.seq)
    if ignoreCase:
        s1 = s1.lower()
        s2 = s2.lower()
    n = 0
    for i in xrange(len(s1)):
        if s1[i] == s2[i]:
            if gapChar is None or s1[i] != gapChar:
                # logger.debug('s1[%d] = %s, s2[%d] = %s', i, s1[i], i, s2[i])
                n = n + 1
    # logger.debug('%s / %s, length: %d, n: %d, gapChar: %s', sr1.id, sr2.id, len(sr1), n, str(gapChar))
    return n


def addGapClassAnnotation(sr):
    """Add letter annotation classifying gaps as internal or terminal.

C{letter_annotations['gapClass']} is set to a list containing C{'t'}
for terminal gaps and C{'i'} for internal gaps. Non-gap symbols are
annotated with C{None}.

@param sr: sequence record to be annotated
@type sr: C{Bio.SeqRecord.SeqRecord}
@return: C{list} of gap class annotations
@rtype: C{list}
    """
    if not isinstance(sr.seq.alphabet, Bio.Alphabet.Gapped):
        raise StandardError, 'sequence is not gapped (not an aligned sequence?)'
    s = str(sr.seq)
    gapClass = [None] * len(sr)
    for i in xrange(len(s)):
        if s[i] == sr.seq.alphabet.gap_char:
            gapClass[i] = 'i'
    i = 0
    while i < len(sr) and gapClass[i] == 'i':
        gapClass[i] = 't'
        i = i + 1
    if i < len(sr):
        i = len(sr) - 1
        while s[i] == sr.seq.alphabet.gap_char:
            gapClass[i] = 't'
            i = i - 1
    sr.letter_annotations['gapClass'] = gapClass
    return gapClass


def pairwiseAlignmentStatsRowDict(alignment):
    rowDict = {}
    if alignment is None:
        rowDict['alignmentLength'] = None
        rowDict['numIdentity'] = None
        rowDict['terminalGapLength1'] = None
        rowDict['terminalGapLength2'] = None
        rowDict['internalGapLength1'] = None
        rowDict['internalGapLength2'] = None
        rowDict['numInternalGaps1'] = None
        rowDict['numInternalGaps2'] = None
    else:
        a1 = alignment[0]
        a2 = alignment[1]
        a1gc = addGapClassAnnotation(a1)
        a2gc = addGapClassAnnotation(a2)
        rowDict['alignmentLength'] = alignment.get_alignment_length()
        rowDict['numIdentity'] = numIdenticalSymbols(a1, a2)
        rowDict['terminalGapLength1'] = sum([1 if gc == 't' else 0 for gc in a1gc])
        rowDict['terminalGapLength2'] = sum([1 if gc == 't' else 0 for gc in a2gc])
        rowDict['internalGapLength1'] = sum([1 if gc == 'i' else 0 for gc in a1gc])
        rowDict['internalGapLength2'] = sum([1 if gc == 'i' else 0 for gc in a2gc])
        rowDict['numInternalGaps1'] = sum([1 if a1gc[i] != 'i' and a1gc[i + 1] == 'i' else 0 for i in xrange(len(a1gc) - 1)])
        rowDict['numInternalGaps2'] = sum([1 if a2gc[i] != 'i' and a2gc[i + 1] == 'i' else 0 for i in xrange(len(a2gc) - 1)])
    return rowDict


def pairwiseAlignmentStats(sr1Dict, sr2Dict, alignmentRunner, alignmentFastaFname=None):
    alignmentStatsFrame = DataFrame(['seqKey', 'seqId1', 'seqId2', 'seqLength1', 'seqLength2', 'alignmentLength', 'numIdentity', 'terminalGapLength1', 'terminalGapLength2', 'internalGapLength1', 'internalGapLength2', 'numInternalGaps1', 'numInternalGaps2'])
    alignmentSaveList = []
    keySet = set(sr1Dict.keys() + sr2Dict.keys())
    for k in keySet:
        sr1 = None
        sr2 = None
        if k in sr1Dict:
            sr1 = sr1Dict[k]
        if k in sr2Dict:
            sr2 = sr2Dict[k]
        if sr1 is not None and sr2 is not None:
            alignmentList = alignmentRunner.align(sr1, sr2)
            alignment = alignmentList[0]
            alignmentSaveList.append(alignment)
        else:
            alignment = None
        rowDict = pairwiseAlignmentStatsRowDict(alignment)
        if sr1 is not None:
            rowDict['seqId1'] = sr1.id
            rowDict['seqLength1'] = len(sr1)
        else:
            rowDict['seqId1'] = None
            rowDict['seqLength1'] = None
        if sr2 is not None:
            rowDict['seqId2'] = sr2.id
            rowDict['seqLength2'] = len(sr2)
        else:
            rowDict['seqId2'] = None
            rowDict['seqLength2'] = None
        rowDict['seqKey'] = k
        alignmentStatsFrame.addRow(rowDict)
    if alignmentFastaFname is not None:
        Bio.AlignIO.write(alignmentSaveList, alignmentFastaFname, 'fasta')
    return alignmentStatsFrame


class PairwiseAlignmentRunner(object):

    def __init__(self):
        pass

    def align(self, sra, srbList):
        raise StandardError, 'abstract method'


class NeedleRunner(PairwiseAlignmentRunner):

    def __init__(self):
        super(NeedleRunner, self).__init__()

    def align(self, sra, srbList):
        logger.debug('starting')
        # FIXME: hard-coded to generate temporary bsequence file in cwd
        bsequenceFd, bsequenceFname = tempfile.mkstemp('.fasta', 'needle_b', '.')
        try:
            bsequenceFile = os.fdopen(bsequenceFd, 'w')
            Bio.SeqIO.write(srbList, bsequenceFile, 'fasta')
            bsequenceFile.close()
            # sys.stderr.write('wrote bsequence file %s\n' % bsequenceFname)
            needleArgv = ['needle', '-asequence', 'stdin', '-bsequence', bsequenceFname, '-outfile', 'stdout', '-auto']
            logger.debug('%s', ' '.join(needleArgv))
            needleProcess = subprocess.Popen(needleArgv, stdin=subprocess.PIPE, stdout=subprocess.PIPE)
            pid = os.fork()
            if pid == 0:
                needleProcess.stdout.close()
                Bio.SeqIO.write([sra], needleProcess.stdin, 'fasta')
                needleProcess.stdin.close()
                os._exit(0)
            needleProcess.stdin.close()
            alignmentList = list(Bio.AlignIO.parse(needleProcess.stdout, 'emboss', alphabet=Bio.Alphabet.Gapped(sra.seq.alphabet)))
            needleProcess.stdout.close()
            wPid, wExit = os.waitpid(pid, 0)
            if pid != wPid:
                raise StandardError('wait returned pid %s (expected %d)' % (wPid, pid))
            if wExit != 0:
                raise StandardError('wait on forked process returned %d' % wExit)
            r = needleProcess.wait()
            if r != 0:
                raise StandardError('needle process exited with %d' % r)
        finally:
            if paftol.keepTmp:
                logger.warning('not deleting needle bsequence file %s', bsequenceFname)
            else:
                os.unlink(bsequenceFname)
        return alignmentList


class WaterRunner(PairwiseAlignmentRunner):

    def __init__(self):
        super(WaterRunner, self).__init__()

    def align(self, sra, srbList):
        logger.debug('starting')
        # FIXME: hard-coded to generate temporary bsequence file in cwd
        bsequenceFd, bsequenceFname = tempfile.mkstemp('.fasta', 'water_b', '.')
        try:
            bsequenceFile = os.fdopen(bsequenceFd, 'w')
            Bio.SeqIO.write(srbList, bsequenceFile, 'fasta')
            bsequenceFile.close()
            # sys.stderr.write('wrote bsequence file %s\n' % bsequenceFname)
            waterArgv = ['water', '-asequence', 'stdin', '-bsequence', bsequenceFname, '-outfile', 'stdout', '-auto']
            logger.debug('%s', ' '.join(waterArgv))
            waterProcess = subprocess.Popen(waterArgv, stdin=subprocess.PIPE, stdout=subprocess.PIPE)
            pid = os.fork()
            if pid == 0:
                waterProcess.stdout.close()
                Bio.SeqIO.write([sra], waterProcess.stdin, 'fasta')
                waterProcess.stdin.close()
                os._exit(0)
            waterProcess.stdin.close()
            alignmentList = list(Bio.AlignIO.parse(waterProcess.stdout, 'emboss', alphabet=Bio.Alphabet.Gapped(sra.seq.alphabet)))
            waterProcess.stdout.close()
            wPid, wExit = os.waitpid(pid, 0)
            if pid != wPid:
                raise StandardError('wait returned pid %s (expected %d)' % (wPid, pid))
            if wExit != 0:
                raise StandardError('wait on forked process returned %d' % wExit)
            r = waterProcess.wait()
            if r != 0:
                raise StandardError('water process exited with %d' % r)
        finally:
            if paftol.keepTmp:
                logger.warning('not deleting water bsequence file %s', bsequenceFname)
            else:
                os.unlink(bsequenceFname)
        return alignmentList


class SemiglobalAlignmentRunner(PairwiseAlignmentRunner):

    def __init__(self):
        pass

    def align(self, sra, srbList):
        alignmentList = []
        sa = str(sra.seq)
        for srb in srbList:
            sb = str(srb.seq)
            aa, ab, alignmentScore = paftol.clib.align_semiglobal(sa, sb)
            asra = Bio.SeqRecord.SeqRecord(Bio.Seq.Seq(aa), id=sra.id, description='%s, aligned semiglobally, score %f' % (sra.description, alignmentScore))
            asrb = Bio.SeqRecord.SeqRecord(Bio.Seq.Seq(ab), id=srb.id, description='%s, aligned semiglobally, score %f' % (srb.description, alignmentScore))
            alignmentList.append(Bio.Align.MultipleSeqAlignment([asra, asrb]))
        return alignmentList


def findRelativeIdentity(alignment):
    n = 0
    for i in xrange(alignment.get_alignment_length()):
        if len(set(alignment[:, i])) == 1:
            n = n + 1
    return float(n) / float(alignment.get_alignment_length())


def findMaxRelativeIdentity(alignment, windowSize):
    maxNumIdentities = None
    numSymbolsList = [len(set(alignment[:, i])) for i in xrange(alignment.get_alignment_length())]
    constColumn = [1 if numSymbolsList[i] == 1 else 0 for i in xrange(alignment.get_alignment_length())]
    numIdentities = 0
    for i in xrange(windowSize - 1):
        numIdentities = numIdentities + constColumn[i]
    for i in xrange(alignment.get_alignment_length() - windowSize + 1):
        numIdentities = numIdentities + constColumn[i + windowSize - 1]
        # sys.stderr.write('i = %d, numIdentities = %d\n' % (i, numIdentities))
        if maxNumIdentities is None or maxNumIdentities < numIdentities:
            maxNumIdentities = numIdentities
        numIdentities = numIdentities - constColumn[i]
    return float(maxNumIdentities) / float(windowSize)


# FIXME: make proper test of this
def testMaxRelativeIdentity(alignment):
    Bio.AlignIO.write(alignment, sys.stderr, 'fasta')
    m = alignment.get_alignment_length() / 2
    for i in xrange(m):
        windowSize = i + 1
        sys.stderr.write('max. relative identity at window %d: %f\n' % (windowSize, findMaxRelativeIdentity(alignment, windowSize)))


class PositionedRead(object):

    def __init__(self, readSr, position, maxRelativeIdentity, coreLength, coreMatch):
        self.position = position
        self.maxRelativeIdentity = maxRelativeIdentity
        self.coreLength = coreLength
        self.coreMatch = coreMatch
        self.readSr = copy.deepcopy(readSr)
        self.readSr.description = 'pos: %d' % self.position

    def __cmp__(self, other):
        return cmp(self.position, other.position)


# FIXME: tidy up these ad-hoc designed helper functions?
def alignSemiglobal(sr0, sr1):
    s0 = str(sr0.seq)
    s1 = str(sr1.seq)
    a0, a1, alignmentScore = paftol.clib.align_semiglobal(s0, s1)
    asr0 = Bio.SeqRecord.SeqRecord(Bio.Seq.Seq(a0), id=sr0.id, description='%s, aligned semiglobally, score %f' % (sr0.description, alignmentScore))
    asr1 = Bio.SeqRecord.SeqRecord(Bio.Seq.Seq(a1), id=sr1.id, description='%s, aligned semiglobally, score %f' % (sr1.description, alignmentScore))
    return Bio.Align.MultipleSeqAlignment([asr0, asr1])


def semiglobalOneVsAll(sr0, sr1List):
    alignmentList = []
    for sr1 in sr1List:
        alignmentList.append(alignSemiglobal(sr0, sr1))
    return alignmentList


def reverseComplementSeqRecordList(seqRecordList):
    rList = []
    for seqRecord in seqRecordList:
        r = seqRecord.reverse_complement()
        r.id = '%s-rc' % seqRecord.id
        r.description = '%s, reverse complement' % seqRecord.description
    return rList


def findFirstNongapPosition(seqRecord, gapChar='-'):
    s = str(seqRecord.seq)
    # logger.debug(s)
    p = 0
    while p < len(s) and s[p] == gapChar:
        p = p + 1
    if p == len(s):
        return None
    else:
        return p


def findLastNongapPosition(seqRecord, gapChar='-'):
    s = str(seqRecord.seq)
    p = -1
    while p > -len(s) and s[p] == gapChar:
        p = p - 1
    if p == -len(s) - 1:
        return None
    else:
        return p + len(s)


def findReadPosition(alignment):
    p = findFirstNongapPosition(alignment[1])
    if p > 0:
        return p
    else:
        return -findFirstNongapPosition(alignment[0])


def findOverlapAlignment(alignment):
    gapChar = '-'
    if len(alignment) != 2:
        raise StandardError, 'pairwise alignment required but this one has %d sequences' % len(alignment)
    l = 0
    if alignment[0, 0] == gapChar:
        l = findFirstNongapPosition(alignment[0], gapChar)
    elif alignment[1, 0] == gapChar:
        l = findFirstNongapPosition(alignment[1], gapChar)
    r = alignment.get_alignment_length() - 1
    if alignment[0, -1] == gapChar:
        r = findLastNongapPosition(alignment[0], gapChar)
    elif alignment[1, -1] == gapChar:
        r = findLastNongapPosition(alignment[1], gapChar)
    # logger.debug('l = %d, r = %d', l, r)
    return alignment[:, l:(r + 1)]


class ContigColumn(object):

    def __init__(self, numRows=0, symbol=None):
        self.symbolList = []
        for i in xrange(numRows):
            self.addRow(symbol)

    def getNumRows(self):
        return len(self.symbolList)

    def addRow(self, symbol):
        self.symbolList.append(symbol)

    def getSymbol(self, rowIndex):
        return self.symbolList[rowIndex]

    def setSymbol(self, rowIndex, symbol):
        self.symbolList[rowIndex] = symbol

    def getNumNongaps(self, gapChar):
        return sum([1 if symbol is not None and symbol != gapChar else 0 for symbol in self.symbolList])

    def getMostFrequentSymbolList(self):
        frequencyDict = {}
        for symbol in self.symbolList:
            if symbol is not None:
                if symbol not in frequencyDict:
                    frequencyDict[symbol] = 0
                frequencyDict[symbol] = frequencyDict[symbol] + 1
        maxFrequency = max(frequencyDict.values())
        mfSymbolList = []
        for symbol in frequencyDict.keys():
            if frequencyDict[symbol] == maxFrequency:
                mfSymbolList.append(symbol)
        return mfSymbolList


class Contig(object):

    def __init__(self, overlapLengthThreshold, overlapMatchThreshold, alignmentRunner, gapChar='-'):
        self.overlapLengthThreshold = overlapLengthThreshold
        self.overlapMatchThreshold = overlapMatchThreshold
        self.alignmentRunner = alignmentRunner
        self.gapChar = gapChar
        self.readList = []
        self.columnList = []

    def numRows(self):
        return len(self.readList)

    def numColumns(self):
        return len(self.columnList)

    def getSymbol(self, rowIndex, columnIndex):
        # logger.debug('rowIndex: %d / %d, columnIndex: %d / %d', rowIndex, self.numRows(), columnIndex, self.numColumns())
        return self.columnList[columnIndex].getSymbol(rowIndex)

    def setSymbol(self, rowIndex, columnIndex, symbol):
        self.columnList[columnIndex].setSymbol(rowIndex, symbol)

    def getSeqRecord(self, rowIndex, terminalGapChar):
        # logger.debug('terminalGapChar: %s', terminalGapChar)
        rawSymbolList = [c.getSymbol(rowIndex) for c in self.columnList]
        # logger.debug('rawSymbolList: %s', str(rawSymbolList))
        symbolList = [rawSymbol if rawSymbol is not None else terminalGapChar for rawSymbol in rawSymbolList]
        s = ''.join(symbolList)
        sr = Bio.SeqRecord.SeqRecord(Bio.Seq.Seq(s), id=self.readList[rowIndex].id, description='')
        return sr

    def getAlignment(self, terminalGapChar=None):
        # logger.debug('terminalGapChar: %s', str(terminalGapChar))
        srList = []
        if terminalGapChar is None:
            terminalGapChar = self.gapChar
        # logger.debug('terminalGapChar: %s', str(terminalGapChar))
        for r in xrange(self.getNumReads()):
            srList.append(self.getSeqRecord(r, terminalGapChar))
        return Bio.Align.MultipleSeqAlignment(srList)

    def getDepthProfile(self):
        return [column.getNumNongaps(self.gapChar) for column in self.columnList]

    def getMeanDepth(self):
        return float(sum(self.getDepthProfile())) / float(self.numColumns())

    def getConsensus(self):
        if len(self.columnList) == 0:
            return None
        s = ''
        depthProfile = []
        for column in self.columnList:
            symbolList = column.getMostFrequentSymbolList()
            symbol = symbolList[0]

            if symbol == self.gapChar and len(symbolList) > 1:
                symbol = symbolList[1]
            if symbol != self.gapChar:
                s = s + symbol
                depthProfile.append(column.getNumNongaps(self.gapChar))
        sr = Bio.SeqRecord.SeqRecord(Bio.Seq.Seq(s, alphabet=Bio.Alphabet.IUPAC.ambiguous_dna), id='contig', description='numReads=%s, meanDepth=%f' % (self.numRows(), float(sum(depthProfile)) / float(len(depthProfile))))
        sr.letter_annotations['depth'] = depthProfile
        return sr

    def getNumReads(self):
        if self.numColumns() == 0:
            return None
        else:
            return self.columnList[0].getNumRows()

    def getLastRead(self):
        return self.readList[-1]

    def findStartPosition(self, rowIndex):
        p = 0
        while self.columnList[p].getSymbol(rowIndex) == self.gapChar:
            # sys.stderr.write('p = %d\n' % p)
            p = p + 1
            if p >= self.numColumns():
                raise StandardError, 'reached end of column list'
        return p

    def insertGapColumn(self, columnIndex=None):
        newColumn = ContigColumn(self.numRows(), self.gapChar)
        if columnIndex is None:
            self.columnList.append(newColumn)
        else:
            self.columnList = self.columnList[:columnIndex] + [newColumn] + self.columnList[columnIndex:]

    def addRow(self):
        for column in self.columnList:
            column.addRow(self.gapChar)

    def removeTerminalGaps(self):
        for rowIndex in xrange(self.numRows()):
            # logger.debug('row: %d', rowIndex)
            columnIndex = 0
            while columnIndex < self.numColumns() and self.getSymbol(rowIndex, columnIndex) == self.gapChar:
                self.setSymbol(rowIndex, columnIndex, None)
                columnIndex = columnIndex + 1
            # logger.debug('final left columnIndex: %d', columnIndex)
            columnIndex = self.numColumns() - 1
            while columnIndex >= 0 and self.getSymbol(rowIndex, columnIndex) == self.gapChar:
                self.setSymbol(rowIndex, columnIndex, None)
                columnIndex = columnIndex - 1
            # logger.debug('final right columnIndex: %d', columnIndex)

    def addFirstRead(self, readSr):
        self.readList.append(readSr)
        for symbol in str(readSr.seq):
            self.columnList.append(ContigColumn(1, symbol))
        return True

    def addSubsequentRead(self, readSr):
        alignmentList = self.alignmentRunner.align(self.getLastRead(), [readSr])
        alignment = alignmentList[0]
        # sys.stderr.write('full alignment:\n')
        # Bio.AlignIO.write(alignment, sys.stderr, 'fasta')
        overlapAlignment = findOverlapAlignment(alignment)
        # Bio.AlignIO.write(overlapAlignment, sys.stderr, 'fasta')
        if overlapAlignment.get_alignment_length() == 0:
            logger.debug('overlapLength = 0, so not adding without further checks')
            return False
        logger.debug('overlapLength: %d, overlapMatch: %f', overlapAlignment.get_alignment_length(), overlapMatch)
        if overlapAlignment.get_alignment_length() >= self.overlapLengthThreshold and overlapMatch >= self.overlapMatchThreshold:
            self.readList.append(readSr)
            r0 = self.getNumReads() - 1
            r1 = r0 + 1
            self.addRow()
            cc = findFirstNongapPosition(alignment[0], self.gapChar)
            c = self.findStartPosition(r0)
            logger.debug('c = %d, cc = %d, adjusted c = %d', c, cc, c - cc)
            c = c - cc
            while c < 0:
                self.insertGapColumn(0)
                c = c + 1
            for i in xrange(alignment.get_alignment_length()):
                # logger.debug('c = %d', c)
                if alignment[0][i] == self.gapChar:
                    # sys.stderr.write('%d: gap in old read\n' % i)
                    if c == self.numColumns():
                        self.insertGapColumn()
                    if self.getSymbol(r0, c) != self.gapChar:
                        self.insertGapColumn(c)
                    self.setSymbol(r1, c, alignment[1][i])
                else:
                    # sys.stderr.write('%d: no gap in old read\n' % i)
                    while self.getSymbol(r0, c) == self.gapChar:
                        c = c + 1
                    self.setSymbol(r1, c, alignment[1][i])
                c = c + 1
            return True
        else:
            return False

    def addRead(self, readSr):
        # sys.stderr.write('adding read %s\n' % readSr.id)
        if self.numRows() == 0:
            isAdded = self.addFirstRead(readSr)
        else:
            isAdded = self.addSubsequentRead(readSr)
        # if isAdded:
            # sys.stderr.write('after adding %s:\n' % readSr.id)
            # Bio.AlignIO.write(self.getAlignment(), sys.stderr, 'fasta')
            # sys.stderr.write('\n')
        return isAdded
