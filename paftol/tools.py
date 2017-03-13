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
import Bio.SeqIO
import Bio.AlignIO
import Bio.Data.CodonTable
import Bio.Blast
import Bio.Blast.NCBIXML


verbose = 0


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
        self.queryId = None
        self.queryDef = None
        self.queryAlignmentStart = None
        self.queryAlignmentEnd = None
        self.queryAlignmentLength = None
        self.queryCdsStart = None
        self.queryCdsEnd = None
        self.queryCdsLength = None
        self.targetId = None
        self.targetDef = None
        self.targetAlignmentStart = None
        self.targetAlignmentEnd = None
        self.targetAlignmentLength = None
        self.targetCdsStart = None
        self.targetCdsEnd = None
        self.targetCdsLength = None
        self.rawScore = None
        self.queryAlignmentSeq = None
        self.queryCdsSeq = None
        self.targetAlignmentSeq = None
        self.targetCdsSeq = None
        self.genomicFragment = None
        self.waterRelativeIdentity = None

    def proteinAlignment(self):
        v = self.vulgar.split()
        qAseq = ''
        tAseq = ''
        qPos = self.queryAlignmentStart
        tPos = self.targetAlignmentStart
        sys.stderr.write('queryAlignmentSeq: %d, targetAlignmentSeq: %d\n' % (len(self.queryAlignmentSeq), len(self.targetAlignmentSeq)))
        for i in xrange(0, len(v), 3):
            vLabel = v[i]
            vQueryLength = int(v[1])
            vTargetLength = int(v[2])
            sys.stderr.write('qPos = %d, tPos = %d, vLabel = %s, vql = %d, vtl = %d\n' % (qPos, tPos, vLabel, vQueryLength, vTargetLength))
            if vLabel == 'M' or vLabel == 'S':
                qAseq = qAseq + str(self.queryAlignmentSeq[qPos:(qPos + vQueryLength)].seq)
                tAseq = tAseq + str(self.targetAlignmentSeq[tPos:(tPos + vTargetLength)].seq)
            elif vLabel == 'G':
                if vQueryLength == 0:
                    if vTargetLength % 3 != 0:
                        raise StandardError, 'cannot process nucleotide gaps with length not a multiple of 3'
                    qAseq = qAseq + '-' * (vTargetLength / 3)
                    tAseq = tAseq + str(self.targetAlignmentSeq[tPos:(tPos + vTargetLength)].seq)
                elif vTargetLength == 0:
                    qAseq = qAseq + str(self.queryAlignmentSeq[qPos:(qPos + vQueryLength)].seq)
                    tAseq = tAseq + '-' * (vQueryLength * 3)
            elif vLabel == '5' or vLabel == '3' or vLabel == 'I' or vLabel == 'F':
                pass
            qPos = qPos + vQueryLength
            tPos = tPos + vTargetLength
            sys.stderr.write('%d: %s\n' % (i, qAseq))
            sys.stderr.write('%d: %s\n' % (i, tAseq))
        
    
    def __str__(self):
        return 'query: %s (%s -> %s), target: %s:%s (%s -> %s), query fragment %dnt' % (self.queryId, self.queryAlignmentStart, self.queryAlignmentEnd, self.targetFname, self.targetId, self.targetAlignmentStart, self.targetAlignmentEnd, len(self.targetCdsSeq))

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
            'rawScore',
            'targetCdsSeqLength', 'genomicFragmentLength', 'waterRelativeIdentity']
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
        d['targetCdsSeqLength'] = None if self.targetCdsSeq is None else len(self.targetCdsSeq)
        d['genomicFragmentLength'] = None if self.genomicFragment is None else len(self.genomicFragment)
        d['waterRelativeIdentity'] = None if self.waterRelativeIdentity is None else self.waterRelativeIdentity
        csvDictWriter.writerow(d)

    def isQueryReversed(self):
        return self.queryAlignmentStart > self.queryAlignmentEnd

    def isTargetReversed(self):
        return self.targetAlignmentStart > self.targetAlignmentEnd

    def makeSeqId(self, seqType):
        return('%s_%s_%s' % (self.queryId, self.targetId, seqType))


class RyoParser(object):
    
    labelledLineRe = re.compile('([A-Za-z][A-Za-z0-9_]*): (.*)')
    seqStartRe = re.compile('seqStart (.*)')
    
    def __init__(self, infile, querySeq, targetFname, exonerateModel):
        self.infile = infile
        self.querySeq = querySeq
        self.targetFname = targetFname
        self.exonerateModel = exonerateModel

    def nextLine(self):
        line = self.infile.readline()
        if verbose > 3:
            sys.stderr.write('nextLine: %s\n' % line.strip())
        if line == '':
            raise StandardError, 'unexpected EOF'
        if line[-1] == '\n':
            line = line[:-1]
        return line
        
    def parseString(self, label):
        line = self.nextLine()
        m = self.labelledLineRe.match(line)
        if m is None:
            raise StandardError, 'malformed line (expected label %s): %s' % (label, line.strip())
        if m.group(1) != label:
            raise StandardError, 'expected label %s but got %s' % (label, m.group(1))
        return m.group(2).strip()
    
    def parseInt(self, label):
        s = self.parseString(label)
        if s == 'NA':
            return None
        return int(s)
    
    def parseSeq(self, label, alphabet, seqId):
        line = self.nextLine()
        m = self.seqStartRe.match(line)
        if m is None:
            raise StandardError, 'malformed line (expected seqStart): %s' % line.strip()
        if m.group(1) != label:
            raise StandardError, 'expected sequence label %s but got %s' % (label, m.group(1))
        seq = ''
        s = self.nextLine()
        while s != 'seqEnd':
            seq = seq + s
            s = self.nextLine()
        return Bio.SeqRecord.SeqRecord(Bio.Seq.Seq(seq, alphabet), id=seqId)
    
    def parse(self):
        line = self.infile.readline()
        if line == '':
            return None
        if line.strip() != 'ryoStart':
            raise StandardError, 'malformed input: ryoStart missing, got %s instead' % line.strip()
        exonerateResult = ExonerateResult(self.querySeq, self.targetFname)
        exonerateResult.queryId = self.parseString('queryId')
        if exonerateResult.queryId != exonerateResult.querySeq.id:
            raise StandardError, 'result incompatible with query: querySeq.id = %s, exonerate queryId = %s' % (self.querySeq.id, exonerateResult.queryId)
        exonerateResult.queryDef = self.parseString('queryDef')
        exonerateResult.queryAlignmentStart = self.parseInt('queryAlignmentStart')
        exonerateResult.queryAlignmentEnd = self.parseInt('queryAlignmentEnd')
        exonerateResult.queryAlignmentLength = self.parseInt('queryAlignmentLength')
        exonerateResult.queryCdsStart = self.parseInt('queryCdsStart')
        exonerateResult.queryCdsEnd = self.parseInt('queryCdsEnd')
        exonerateResult.queryCdsLength = self.parseInt('queryCdsLength')
        exonerateResult.targetId = self.parseString('targetId')
        exonerateResult.targetDef = self.parseString('targetDef')
        exonerateResult.targetAlignmentStart = self.parseInt('targetAlignmentStart')
        exonerateResult.targetAlignmentEnd = self.parseInt('targetAlignmentEnd')
        exonerateResult.targetAlignmentLength = self.parseInt('targetAlignmentLength')
        exonerateResult.targetCdsStart = self.parseInt('targetCdsStart')
        exonerateResult.targetCdsEnd = self.parseInt('targetCdsEnd')
        exonerateResult.targetCdsLength = self.parseInt('targetCdsLength')
        exonerateResult.rawScore = self.parseInt('rawScore')
        exonerateResult.vulgar = self.parseString('vulgar')
        if self.exonerateModel == 'protein2genome':
            # exonerateResult.queryCdsSeq = self.parseSeq('queryCds', Bio.Alphabet.IUPAC.ambiguous_dna, exonerateResult.makeSeqId('qcds'))
            # FIXME: kludge -- throwing away rubbish output of exonerate, would be much better not to generate it in the first place
            self.parseSeq('queryCds', Bio.Alphabet.IUPAC.ambiguous_dna, exonerateResult.makeSeqId('qcds'))
            exonerateResult.queryCdsSeq = None
            exonerateResult.queryAlignmentSeq = self.parseSeq('queryAlignment', Bio.Alphabet.IUPAC.ambiguous_dna, exonerateResult.makeSeqId('qaln'))
            exonerateResult.targetCdsSeq = self.parseSeq('targetCds', Bio.Alphabet.IUPAC.ambiguous_dna, exonerateResult.makeSeqId('tcds'))
            exonerateResult.targetAlignmentSeq = self.parseSeq('targetAlignment', Bio.Alphabet.IUPAC.ambiguous_dna, exonerateResult.makeSeqId('taln'))
        else:
            raise StandardError, 'unsupported exonerate model: %s' % exonerateModel
        line = self.nextLine()
        if line.strip() != 'ryoEnd':
            raise StandardError, 'malformed input: ryoEnd missing'
        return exonerateResult

def findExonerateResultList(query, targetFname, exonerateModel, bestn, keepTmpFile=False):
    exonerateModelList = ['protein2genome']
    if exonerateModel not in exonerateModelList:
        raise StandardError, 'unknown exonerate alignment model: %s' % exonerateModel
    exonerateResultList = []
    queryScratchFd, queryScratchFname = tempfile.mkstemp('.fasta', 'scratch', '.')
    try:
        queryScratchFile = os.fdopen(queryScratchFd, 'w')
        Bio.SeqIO.write(query, queryScratchFile, 'fasta')
        queryScratchFile.close()
        # FIXME: support multiple alignment models, currently protein2genome is hard-coded
        exonerateArgv = [
            'exonerate', '--model', exonerateModel, '--bestn', '%d' % bestn, '--verbose', '0', '--showalignment',
            'no', '--showvulgar', 'no', '--ryo', 'ryoStart\\nqueryId: %qi\\nqueryDef: %qd\\nqueryAlignmentStart: %qab\\nqueryAlignmentEnd: %qae\\nqueryAlignmentLength: %qal\\nqueryCdsStart: NA\\nqueryCdsEnd: NA\\nqueryCdsLength: NA\\ntargetId: %ti\\ntargetDef: %td\\ntargetAlignmentStart: %tab\\ntargetAlignmentEnd: %tae\\ntargetAlignmentLength: %tal\\ntargetCdsStart: %tcb\\ntargetCdsEnd: %tce\\ntargetCdsLength: %tcl\\nrawScore: %s\\nvulgar: %V\\nseqStart queryCds\\n%qcsseqEnd\\nseqStart queryAlignment\\n%qasseqEnd\\nseqStart targetCds\\n%tcsseqEnd\\nseqStart targetAlignment\\n%tasseqEnd\\nryoEnd\\n', queryScratchFname, targetFname]
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
        ryoParser = RyoParser(p.stdout, query, targetFname, exonerateModel)
        exonerateResult = ryoParser.parse()
        while exonerateResult is not None :
            exonerateResultList.append(exonerateResult)
            exonerateResult = ryoParser.parse()
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
        if keepTmpFile:
            sys.stderr.write('not deleting query scratch file %s' % queryScratchFname)
        else:
            os.unlink(queryScratchFname)
    return exonerateResultList
