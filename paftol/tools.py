#!/usr/bin/env python

import sys
import re
import os
import os.path
import subprocess
import csv
import tempfile

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


verbose = 0


def translateGapped(seq, table='Standard'):
    """Translate a gapped sequence.

All triplets in the sequence must either contain only non-gap symbols
or only gap symbols. Runs of non-gap triplets are translated as per
the specified translation table. Runs of gap symbol triplets are
translated into as many single gap symbols. An exception is raised
if the sequence contains triplets with both gap and non-gap symbols.
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
    # sys.stderr.write('gaplist: %s\n' % str(gapList))
    t = ''
    nongapStart = 0
    for gapStart, gapEnd in gapList:
        if (gapStart - nongapStart) % 3 != 0:
            raise StandardError, 'cannot translate: range %d:%d does not coincide with triplet boundaries' % (nongapStart, gapStart)
        if (gapEnd - gapStart) % 3 != 0:
            raise StandardError, 'cannot translate: range %d:%d does not coincide with triplet boundaries' % (gapStart, gapEnd)
        if nongapStart < gapStart:
            ungappedSeq = seq[nongapStart:gapStart]
            ungappedSeq.alphabet = ungappedAlphabet
            t = t + str(ungappedSeq.translate(table))
        t = t + seq.alphabet.gap_char * ((gapEnd - gapStart) / 3)
        nongapStart = gapEnd
    if nongapStart < len(seq):
        if (len(seq) - nongapStart) % 3 != 0:
            raise StandardError, 'cannot translate: range %d:%d does not coincide with triplet boundaries' % (nongapStart, len(seq))
        ungappedSeq = seq[nongapStart:]
        ungappedSeq.alphabet = ungappedAlphabet
        t = t + str(ungappedSeq.translate(table))
    return Bio.Seq.Seq(t, alphabet=Bio.Alphabet.Gapped(Bio.Alphabet.IUPAC.protein))
    

def ascendingRange(rangeStart, rangeEnd):
    if rangeStart > rangeEnd:
        return rangeEnd, rangeStart
    else:
        return rangeStart, rangeEnd
    
    
class ExonerateResult(object):

    queryRe = re.compile('Query: ([^ ]+)( (.*))?')
    targetRe = re.compile('Target: ([^ ]+)( (.*))?')
    queryRangeRe = re.compile('Query range: ([0-9]+) -> ([0-9]+)')
    targetRangeRe = re.compile('Target range: ([0-9]+) -> ([0-9]+)')
    rawScoreRe = re.compile('rawScore: ([0-9]+)')
    seqlineRe = re.compile('[ACGT]*')

    def __init__(self, querySeq, targetFname):
        self.querySeq = querySeq
        self.targetFname = targetFname
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
        self.genomicFragment = None
    
    def queryAlignmentOverlap(self, other):
        selfQas, selfQae = ascendingRange(self.queryAlignmentStart, self.queryAlignmentEnd)
        otherQas, otherQae = ascendingRange(other.queryAlignmentStart, other.queryAlignmentEnd)
        if selfQas > otherQae or selfQae < otherQas:
            return None
        return max(selfQas, otherQas), min(selfQae, otherQae)
        
    def containsQueryAlignmentRange(self, other):
        selfQas, selfQae = ascendingRange(self.queryAlignmentStart, self.queryAlignmentEnd)
        otherQas, otherQae = ascendingRange(other.queryAlignmentStart, other.queryAlignmentEnd)
        return selfQas <= otherQas and selfQae >= otherQae
    
    def overlapsQueryAlignmentRange(self, other):
        return self.queryAlignmentOverlap(other) is not None
    
    def proteinAlignment(self):
        if self.exonerateModel != 'protein2genome:local':
            raise StandardError, 'proteinAlignment is not supported for exonerate model "%s"' % self.exonerateModel
        v = self.vulgar.split()
        qAln = ''
        tAln = ''
        qPos = 0
        tPos = 0
        # sys.stderr.write('queryAlignmentSeq: %d, targetAlignmentSeq: %d\n' % (len(self.queryAlignmentSeq), len(self.targetAlignmentSeq)))
        # sys.stderr.write('%s\n' % str(self.targetAlignmentSeq.seq))
        for i in xrange(0, len(v), 3):
            vLabel = v[i]
            vQueryLength = int(v[i + 1])
            vTargetLength = int(v[i + 2])
            # sys.stderr.write('qPos = %d, tPos = %d, vLabel = %s, vql = %d, vtl = %d\n' % (qPos, tPos, vLabel, vQueryLength, vTargetLength))
            if vLabel == 'M' or vLabel == 'S':
                qAln = qAln + str(self.queryAlignmentSeq[qPos:(qPos + vQueryLength)].seq)
                tAln = tAln + str(self.targetAlignmentSeq[tPos:(tPos + vTargetLength)].seq)
            elif vLabel == 'G':
                if vQueryLength == 0:
                    if vTargetLength % 3 != 0:
                        raise StandardError, 'cannot process nucleotide gaps with length not a multiple of 3'
                    qAln = qAln + '-' * (vTargetLength / 3)
                    tAln = tAln + str(self.targetAlignmentSeq[tPos:(tPos + vTargetLength)].seq)
                elif vTargetLength == 0:
                    qAln = qAln + str(self.queryAlignmentSeq[qPos:(qPos + vQueryLength)].seq)
                    tAln = tAln + '-' * (vQueryLength * 3)
            elif vLabel == '5' or vLabel == '3' or vLabel == 'I' or vLabel == 'F':
                pass
            qPos = qPos + vQueryLength
            tPos = tPos + vTargetLength
            # sys.stderr.write('%d: qPos = %d, tPos = %d, vLabel = %s, vql = %d, vtl = %d\n' % (i, qPos, tPos, vLabel, vQueryLength, vTargetLength))
            # sys.stderr.write('qAln: %s\n' % qAln)
            # sys.stderr.write('tAln: %s\n' % tAln)
            # s = Bio.Seq.Seq(tAln, alphabet=Bio.Alphabet.Gapped(self.targetAlignmentSeq.seq.alphabet))
            # s3 = s[:(len(s) - len(s) % 3)]
            # sys.stderr.write('tA_tr: %s%s\n' % (str(translateGapped(s3)), '' if len(s) == len(s3) else '.'))
            # sys.stderr.write('\n')
        qAlnSeq = Bio.Seq.Seq(qAln, alphabet=Bio.Alphabet.Gapped(Bio.Alphabet.IUPAC.protein))
        tAlnSeq = Bio.Seq.Seq(tAln, alphabet=Bio.Alphabet.Gapped(self.targetAlignmentSeq.seq.alphabet))
        # FIXME: assuming standard translation table -- check whether exonerate supports setting table?
        tAlnProt = translateGapped(tAlnSeq)
        return Bio.Align.MultipleSeqAlignment([Bio.SeqRecord.SeqRecord(qAlnSeq, id=self.queryAlignmentSeq.id), Bio.SeqRecord.SeqRecord(tAlnProt, id='%s_pep' % self.targetAlignmentSeq.id)])
    
    def __str__(self):
        return '%s, query: %s %s(%s -> %s), target: %s:%s %s(%s -> %s), query fragment %dnt' % (self.exonerateModel, self.queryId, self.queryStrand, self.queryAlignmentStart, self.queryAlignmentEnd, self.targetFname, self.targetId, self.targetStrand, self.targetAlignmentStart, self.targetAlignmentEnd, len(self.targetCdsSeq))

    def isEmpty(self):
        return self.queryId is None
        
    def writeAllSeqs(self, f):
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

    def makeCsvDictWriter(self, csvfile):
        # FIXME: add columns for all attributes
        csvFieldnames = [
            'querySeqId', 'genomeName', 'querySeqLength',
            'queryId',
            'queryAlignmentStart', 'queryAlignmentEnd', 'queryAlignmentLength',
            'queryCdsStart', 'queryCdsEnd', 'queryCdsLength',
            'targetId',
            'targetAlignmentStart', 'targetAlignmentEnd', 'targetAlignmentLength',
            'targetCdsStart', 'targetCdsEnd', 'targetCdsLength',
            'rawScore','percentIdentity','percentSimilarity',
            'equivalencedTotal', 'equivalencedIdentity', 'equivalencedSimilarity', 'equivalencedMismatches', 
            'targetCdsSeqLength', 'genomicFragmentLength']
        csvDictWriter = csv.DictWriter(csvfile, csvFieldnames)
        csvDictWriter.writeheader()
        return csvDictWriter

    def writeCsvStats(self, csvDictWriter):
        d = {}
        d['querySeqId'] = None if self.querySeq is None else self.querySeq.id
        d['targetFname'] = self.targetFname
        d['querySeqLength'] = None if self.querySeq is None else len(self.querySeq)
        d['queryId'] = self.queryId
        d['queryAlignmentStart'] = self.queryAlignmentStart
        d['queryAlignmentEnd'] = self.queryAlignmentEnd
        d['queryAlignmentLength'] = self.queryAlignmentLength
        d['queryCdsStart'] = self.queryCdsStart
        d['queryCdsEnd'] = self.queryCdsEnd
        d['queryCdsLength'] = self.queryCdsLength
        d['targetId'] = self.targetId
        d['targetAlignmentStart'] = self.targetAlignmentStart
        d['targetAlignmentEnd'] = self.targetAlignmentEnd
        d['targetAlignmentLength'] = self.targetAlignmentLength
        d['targetCdsStart'] = self.targetCdsStart
        d['targetCdsEnd'] = self.targetCdsEnd
        d['targetCdsLength'] = self.targetCdsLength
        d['rawScore'] = self.rawScore
        d['percentIdentity'] = self.percentIdentity
        d['percentSimilarity'] = self.percentSimilarity
        d['equivalencedTotal'] = self.equivalencedTotal
        d['equivalencedIdentity'] = self.equivalencedIdentity
        d['equivalencedSimilarity'] = self.equivalencedSimilarity
        d['equivalencedMismatches'] = self.equivalencedMismatches
        d['targetCdsSeqLength'] = None if self.targetCdsSeq is None else len(self.targetCdsSeq)
        d['genomicFragmentLength'] = None if self.genomicFragment is None else len(self.genomicFragment)
        csvDictWriter.writerow(d)

    def isQueryReversed(self):
        return self.queryAlignmentStart > self.queryAlignmentEnd

    def isTargetReversed(self):
        return self.targetAlignmentStart > self.targetAlignmentEnd

    def makeSeqId(self, seqType):
        return('%s_%s_%s' % (self.queryId, self.targetId, seqType))


class ExonerateRunner(object):
    
    labelledLineRe = re.compile('([A-Za-z][A-Za-z0-9_]*): (.*)')
    seqStartRe = re.compile('seqStart (.*)')
    
    def __init__(self):
        pass
        
    def nextLine(self, f):
        line = f.readline()
        if verbose > 3:
            sys.stderr.write('nextLine: %s\n' % line.strip())
        if line == '':
            raise StandardError, 'unexpected EOF'
        if line[-1] == '\n':
            line = line[:-1]
        return line
        
    def parseString(self, f, label):
        line = self.nextLine(f)
        m = self.labelledLineRe.match(line)
        if m is None:
            raise StandardError, 'malformed line (expected label %s): %s' % (label, line.strip())
        if m.group(1) != label:
            raise StandardError, 'expected label %s but got %s' % (label, m.group(1))
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
            raise StandardError, 'malformed line (expected seqStart): %s' % line.strip()
        if m.group(1) != label:
            raise StandardError, 'expected sequence label %s but got %s' % (label, m.group(1))
        seq = ''
        s = self.nextLine(f)
        while s != 'seqEnd':
            seq = seq + s
            s = self.nextLine(f)
        return Bio.SeqRecord.SeqRecord(Bio.Seq.Seq(seq, alphabet), id=seqId)
    
    def parseExonerateResult(self, f, exonerateResult):
        line = f.readline()
        if line == '':
            return None
        if line.strip() != 'ryoStart':
            raise StandardError, 'malformed input: ryoStart missing, got %s instead' % line.strip()
        exonerateResult.exonerateModel = self.parseString(f, 'exonerateModel')
        exonerateResult.queryId = self.parseString(f, 'queryId')
        if exonerateResult.queryId != exonerateResult.querySeq.id:
            raise StandardError, 'result incompatible with query: querySeq.id = %s, exonerate queryId = %s' % (self.querySeq.id, exonerateResult.queryId)
        exonerateResult.queryDef = self.parseString(f, 'queryDef')
        exonerateResult.queryStrand = self.parseString(f, 'queryStrand')
        exonerateResult.queryAlignmentStart = self.parseInt(f, 'queryAlignmentStart')
        exonerateResult.queryAlignmentEnd = self.parseInt(f, 'queryAlignmentEnd')
        exonerateResult.queryAlignmentLength = self.parseInt(f, 'queryAlignmentLength')
        exonerateResult.queryCdsStart = self.parseInt(f, 'queryCdsStart')
        exonerateResult.queryCdsEnd = self.parseInt(f, 'queryCdsEnd')
        exonerateResult.queryCdsLength = self.parseInt(f, 'queryCdsLength')
        exonerateResult.targetId = self.parseString(f, 'targetId')
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
        if exonerateResult.exonerateModel == 'protein2genome:local':
            # exonerateResult.queryCdsSeq = self.parseSeq(f, 'queryCds', Bio.Alphabet.IUPAC.ambiguous_dna, exonerateResult.makeSeqId('qcds'))
            # FIXME: kludge -- throwing away rubbish output of exonerate, would be much better not to generate it in the first place
            self.parseSeq(f, 'queryCds', Bio.Alphabet.IUPAC.ambiguous_dna, exonerateResult.makeSeqId('qcds'))
            exonerateResult.queryCdsSeq = None
            exonerateResult.queryAlignmentSeq = self.parseSeq(f, 'queryAlignment', Bio.Alphabet.IUPAC.protein, exonerateResult.makeSeqId('qaln'))
            exonerateResult.targetCdsSeq = self.parseSeq(f, 'targetCds', Bio.Alphabet.IUPAC.ambiguous_dna, exonerateResult.makeSeqId('tcds'))
            exonerateResult.targetAlignmentSeq = self.parseSeq(f, 'targetAlignment', Bio.Alphabet.IUPAC.ambiguous_dna, exonerateResult.makeSeqId('taln'))
        else:
            raise StandardError, 'unsupported exonerate model: %s' % exonerateResult.exonerateModel
        line = self.nextLine(f)
        if line.strip() != 'ryoEnd':
            raise StandardError, 'malformed input: ryoEnd missing'
        return exonerateResult

    def parse(self, querySeq, targetFname, exonerateModel, bestn):
        queryScratchFd, queryScratchFname = tempfile.mkstemp('.fasta', 'scratch', '.')
        try:
            queryScratchFile = os.fdopen(queryScratchFd, 'w')
            Bio.SeqIO.write(querySeq, queryScratchFile, 'fasta')
            queryScratchFile.close()
            exonerateArgv = [
                'exonerate', '--model', exonerateModel, '--bestn', '%d' % bestn, '--verbose', '0', '--showalignment',
                'no', '--showvulgar', 'no', '--ryo', 'ryoStart\\nexonerateModel: %m\\nqueryId: %qi\\nqueryDef: %qd\\nqueryStrand: %qS\\nqueryAlignmentStart: %qab\\nqueryAlignmentEnd: %qae\\nqueryAlignmentLength: %qal\\nqueryCdsStart: NA\\nqueryCdsEnd: NA\\nqueryCdsLength: NA\\ntargetId: %ti\\ntargetDef: %td\\ntargetStrand: %tS\\ntargetAlignmentStart: %tab\\ntargetAlignmentEnd: %tae\\ntargetAlignmentLength: %tal\\ntargetCdsStart: %tcb\\ntargetCdsEnd: %tce\\ntargetCdsLength: %tcl\\nrawScore: %s\\npercentIdentity: %pi\\npercentSimilarity: %ps\\nequivalencedTotal: %et\\nequivalencedIdentity: %ei\\nequivalencedSimilarity: %es\\nequivalencedMismatches: %em\\nvulgar: %V\\nseqStart queryCds\\n%qcsseqEnd\\nseqStart queryAlignment\\n%qasseqEnd\\nseqStart targetCds\\n%tcsseqEnd\\nseqStart targetAlignment\\n%tasseqEnd\\nryoEnd\\n', queryScratchFname, targetFname]
            if verbose > 0:
                sys.stderr.write('%s\n' % ' '.join(exonerateArgv))
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
            exonerateResult = self.parseExonerateResult(p.stdout, ExonerateResult(querySeq, targetFname))
            while exonerateResult is not None:
                exonerateResultList.append(exonerateResult)
                exonerateResult = self.parseExonerateResult(p.stdout, ExonerateResult(querySeq, targetFname))
            p.stdout.close()
            wPid, wExit = os.waitpid(pid, 0)
            if pid != wPid:
                raise StandardError, 'wait returned pid %s (expected %d)' % (wPid, pid)
            if wExit != 0:
                raise StandardError, 'wait on forked process returned %d' % wExit
            r = p.wait()
            if r != 0:
                raise StandardError, 'exonerate process exited with %d' % r
        finally:
            if paftol.keepTmp:
                sys.stderr.write('not deleting query scratch file %s' % queryScratchFname)
            else:
                os.unlink(queryScratchFname)
        return exonerateResultList
