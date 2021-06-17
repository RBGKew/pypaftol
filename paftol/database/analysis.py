#!/usr/bin/env python

# Copyright (c) 2020 The Board of Trustees of the Royal Botanic Gardens, Kew

import sys
import getopt
import re
import unicodedata

import mysql.connector

import paftol.database


class AnnotatedGenome(object):

    def __init__(self, idSequencing=None, accessionId=None, speciesLatinName=None, commonName=None, source=None, genomeVersion=None):
        self.id = None
        self.idSequencing = idSequencing
        self.accessionId = accessionId
        self.speciesLatinName = speciesLatinName
        self.commonName = commonName
        self.source = source
        self.genomeVersion = genomeVersion
        # one-to-many
        # fk_InputSequence_AnnotatedGenomeId: InputSequence.annotatedGenomeId REFERENCES AnnotatedGenome(annotatedGenomeId)
        self.inputSequenceAnnotatedGenomeList = []

    def insertIntoDatabase(self, cursor):
        sqlCmd = 'INSERT INTO `AnnotatedGenome` (`idSequencing`, `accessionId`, `speciesLatinName`, `commonName`, `source`, `genomeVersion`) VALUES (%s, %s, %s, %s, %s, %s)'
        l = []
        l.append(self.idSequencing)
        l.append(self.accessionId)
        l.append(self.speciesLatinName)
        l.append(self.commonName)
        l.append(self.source)
        l.append(self.genomeVersion)
        cursor.execute(sqlCmd, tuple(l))


class ContigRecovery(object):

    def __init__(self, fwdFastq=None, revFastq=None, fwdTrimmedFastqStats=None, revTrimmedFastqStats=None, contigFastaFileName=None, contigFastaFilePathName=None, contigFastaFileMd5sum=None, referenceTarget=None, numMappedReads=None, numUnmappedReads=None, softwareVersion=None, cmdLine=None, numRecoveredContigsCheck=None):
        self.id = None
        self.fwdFastq = fwdFastq
        self.revFastq = revFastq
        self.fwdTrimmedFastqStats = fwdTrimmedFastqStats
        self.revTrimmedFastqStats = revTrimmedFastqStats
        self.contigFastaFileName = contigFastaFileName
        self.contigFastaFilePathName = contigFastaFilePathName
        self.contigFastaFileMd5sum = contigFastaFileMd5sum
        self.referenceTarget = referenceTarget
        self.numMappedReads = numMappedReads
        self.numUnmappedReads = numUnmappedReads
        self.softwareVersion = softwareVersion
        self.cmdLine = cmdLine
        self.numRecoveredContigsCheck = numRecoveredContigsCheck
        # one-to-many
        # fk_ContigRecoveryDataRelease_contigRecoveryId: ContigRecoveryDataRelease.contigRecoveryId REFERENCES ContigRecovery(contigRecoveryId)
        self.contigRecoveryDataReleaseContigRecoveryList = []
        # fk_RecoveredContig_contigRecoveryId: RecoveredContig.contigRecoveryId REFERENCES ContigRecovery(contigRecoveryId)
        self.recoveredContigContigRecoveryList = []

    def insertIntoDatabase(self, cursor):
        sqlCmd = 'INSERT INTO `ContigRecovery` (`fwdFastqId`, `revFastqId`, `fwdTrimmedFastqStatsId`, `revTrimmedFastqStatsId`, `contigFastaFileName`, `contigFastaFilePathName`, `contigFastaFileMd5sum`, `referenceTargetId`, `numMappedReads`, `numUnmappedReads`, `softwareVersion`, `cmdLine`, `numRecoveredContigsCheck`) VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s)'
        l = []
        l.append(None if self.fwdFastq is None else self.fwdFastq.id)
        l.append(None if self.revFastq is None else self.revFastq.id)
        l.append(None if self.fwdTrimmedFastqStats is None else self.fwdTrimmedFastqStats.id)
        l.append(None if self.revTrimmedFastqStats is None else self.revTrimmedFastqStats.id)
        l.append(self.contigFastaFileName)
        l.append(self.contigFastaFilePathName)
        l.append(self.contigFastaFileMd5sum)
        l.append(None if self.referenceTarget is None else self.referenceTarget.id)
        l.append(self.numMappedReads)
        l.append(self.numUnmappedReads)
        l.append(self.softwareVersion)
        l.append(self.cmdLine)
        l.append(self.numRecoveredContigsCheck)
        cursor.execute(sqlCmd, tuple(l))


class ContigRecoveryDataRelease(object):

    def __init__(self, contigRecovery=None, dataRelease=None):
        self.id = None
        self.contigRecovery = contigRecovery
        self.dataRelease = dataRelease
        # one-to-many

    def insertIntoDatabase(self, cursor):
        sqlCmd = 'INSERT INTO `ContigRecoveryDataRelease` (`contigRecoveryId`, `dataReleaseId`) VALUES (%s, %s)'
        l = []
        l.append(None if self.contigRecovery is None else self.contigRecovery.id)
        l.append(None if self.dataRelease is None else self.dataRelease.idDataRelease)
        cursor.execute(sqlCmd, tuple(l))


class DataOrigin(object):

    def __init__(self, dataOriginName=None, acronym=None):
        self.id = None
        self.dataOriginName = dataOriginName
        self.acronym = acronym
        # one-to-many
        # fk_InputSequence_dataOriginId: InputSequence.dataOriginId REFERENCES DataOrigin(dataOriginId)
        self.inputSequenceDataOriginList = []

    def insertIntoDatabase(self, cursor):
        sqlCmd = 'INSERT INTO `DataOrigin` (`dataOriginName`, `acronym`) VALUES (%s, %s)'
        l = []
        l.append(self.dataOriginName)
        l.append(self.acronym)
        cursor.execute(sqlCmd, tuple(l))


class DataRelease(object):

    def __init__(self, releaseNumber=None, dataRelease=None):
        self.idDataRelease = None
        self.releaseNumber = releaseNumber
        self.dataRelease = dataRelease
        # one-to-many
        # fk_ContigRecoveryDataRelease_dataReleaseId: ContigRecoveryDataRelease.dataReleaseId REFERENCES DataRelease(dataReleaseId)
        self.contigRecoveryDataReleaseDataReleaseList = []
        # fk_GeneTreeDataRelease_dataReleaseId: GeneTreeDataRelease.dataReleaseId REFERENCES DataRelease(dataReleaseId)
        self.geneTreeDataReleaseDataReleaseList = []
        # fk_OneKP_SequenceDataRelease_dataReleaseId: OneKP_SequenceDataRelease.dataReleaseId REFERENCES DataRelease(dataReleaseId)
        self.oneKP_SequenceDataReleaseDataReleaseList = []
        # fk_SRA_RunSequenceDataRelease_dataReleaseId: SRA_RunSequenceDataRelease.dataReleaseId REFERENCES DataRelease(dataReleaseId)
        self.srA_RunSequenceDataReleaseDataReleaseList = []

    def insertIntoDatabase(self, cursor):
        sqlCmd = 'INSERT INTO `DataRelease` (`ReleaseNumber`, `DataRelease`) VALUES (%s, %s)'
        l = []
        l.append(self.releaseNumber)
        l.append(self.dataRelease)
        cursor.execute(sqlCmd, tuple(l))


class ENA_Accession(object):

    def __init__(self, accessionId=None):
        self.id = None
        self.accessionId = accessionId
        # one-to-many
        # fk_SRA_RunSequence_enaAccessionId: SRA_RunSequence.enaAccessionId REFERENCES ENA_Accession(enaAccessionId)
        self.srA_RunSequenceEnaAccessionList = []

    def insertIntoDatabase(self, cursor):
        sqlCmd = 'INSERT INTO `ENA_Accession` (`accessionId`) VALUES (%s)'
        l = []
        l.append(self.accessionId)
        cursor.execute(sqlCmd, tuple(l))


class ExemplarGene(object):

    def __init__(self, ac=None, gn=None, de=None, os=None, url=None):
        self.id = None
        self.ac = ac
        self.gn = gn
        self.de = de
        self.os = os
        self.url = url
        # one-to-many
        # fk_PaftolGene_exemplarGeneId: PaftolGene.exemplarGeneId REFERENCES ExemplarGene(exemplarGeneId)
        self.paftolGeneExemplarGeneList = []

    def insertIntoDatabase(self, cursor):
        sqlCmd = 'INSERT INTO `ExemplarGene` (`AC`, `GN`, `DE`, `OS`, `URL`) VALUES (%s, %s, %s, %s, %s)'
        l = []
        l.append(self.ac)
        l.append(self.gn)
        l.append(self.de)
        l.append(self.os)
        l.append(self.url)
        cursor.execute(sqlCmd, tuple(l))


class FastqStats(object):

    def __init__(self, numReads=None, numbrRecords=None, qual28=None, meanA=None, meanC=None, meanG=None, meanT=None, stddevA=None, stddevC=None, stddevG=None, stddevT=None, meanN=None, stddevN=None, meanAdapterContent=None, maxAdapterContent=None, sumLengthOfSeqs=None):
        self.id = None
        self.numReads = numReads
        self.numbrRecords = numbrRecords
        self.qual28 = qual28
        self.meanA = meanA
        self.meanC = meanC
        self.meanG = meanG
        self.meanT = meanT
        self.stddevA = stddevA
        self.stddevC = stddevC
        self.stddevG = stddevG
        self.stddevT = stddevT
        self.meanN = meanN
        self.stddevN = stddevN
        self.meanAdapterContent = meanAdapterContent
        self.maxAdapterContent = maxAdapterContent
        self.sumLengthOfSeqs = sumLengthOfSeqs
        # one-to-many
        # fk_ContigRecovery_fwdTrimmedFastqStatsId: ContigRecovery.fwdTrimmedFastqStatsId REFERENCES FastqStats(fwdTrimmedFastqStatsId)
        self.contigRecoveryFwdTrimmedFastqStatsList = []
        # fk_ContigRecovery_revTrimmedFastqStatsId: ContigRecovery.revTrimmedFastqStatsId REFERENCES FastqStats(revTrimmedFastqStatsId)
        self.contigRecoveryRevTrimmedFastqStatsList = []
        # fk_InputSequence_fastqStatsId: InputSequence.fastqStatsId REFERENCES FastqStats(fastqStatsId)
        self.inputSequenceFastqStatsList = []

    def insertIntoDatabase(self, cursor):
        sqlCmd = 'INSERT INTO `FastqStats` (`numReads`, `numbrRecords`, `qual28`, `meanA`, `meanC`, `meanG`, `meanT`, `stddevA`, `stddevC`, `stddevG`, `stddevT`, `meanN`, `stddevN`, `meanAdapterContent`, `maxAdapterContent`, `sumLengthOfSeqs`) VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s)'
        l = []
        l.append(self.numReads)
        l.append(self.numbrRecords)
        l.append(self.qual28)
        l.append(self.meanA)
        l.append(self.meanC)
        l.append(self.meanG)
        l.append(self.meanT)
        l.append(self.stddevA)
        l.append(self.stddevC)
        l.append(self.stddevG)
        l.append(self.stddevT)
        l.append(self.meanN)
        l.append(self.stddevN)
        l.append(self.meanAdapterContent)
        l.append(self.maxAdapterContent)
        l.append(self.sumLengthOfSeqs)
        cursor.execute(sqlCmd, tuple(l))


class GAP_Sequence(object):

    def __init__(self, idSequencing=None, sampleId=None):
        self.id = None
        self.idSequencing = idSequencing
        self.sampleId = sampleId
        # one-to-many
        # fk_InputSequence_GAP_SequenceId: InputSequence.GAP_SequenceId REFERENCES GAP_Sequence(GAP_SequenceId)
        self.inputSequenceGAP_SequenceList = []

    def insertIntoDatabase(self, cursor):
        sqlCmd = 'INSERT INTO `GAP_Sequence` (`idSequencing`, `sampleId`) VALUES (%s, %s)'
        l = []
        l.append(self.idSequencing)
        l.append(self.sampleId)
        cursor.execute(sqlCmd, tuple(l))


class GeneTree(object):

    def __init__(self, paftolGene=None, unAlnFastaFile=None, unAlnFastaFilePathName=None, alnFastaFile=None, alnFastaFilePathName=None, newickFile=None, newickFilePathName=None, speciesTree=None, cmdLine=None, softwareVersion=None):
        self.id = None
        self.paftolGene = paftolGene
        self.unAlnFastaFile = unAlnFastaFile
        self.unAlnFastaFilePathName = unAlnFastaFilePathName
        self.alnFastaFile = alnFastaFile
        self.alnFastaFilePathName = alnFastaFilePathName
        self.newickFile = newickFile
        self.newickFilePathName = newickFilePathName
        self.speciesTree = speciesTree
        self.cmdLine = cmdLine
        self.softwareVersion = softwareVersion
        # one-to-many
        # fk_GeneTreeDataRelease_geneTreeId: GeneTreeDataRelease.geneTreeId REFERENCES GeneTree(geneTreeId)
        self.geneTreeDataReleaseGeneTreeList = []

    def insertIntoDatabase(self, cursor):
        sqlCmd = 'INSERT INTO `GeneTree` (`paftolGeneId`, `unAlnFastaFile`, `unAlnFastaFilePathName`, `alnFastaFile`, `alnFastaFilePathName`, `newickFile`, `newickFilePathName`, `speciesTreeId`, `cmdLine`, `softwareVersion`) VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s, %s)'
        l = []
        l.append(None if self.paftolGene is None else self.paftolGene.id)
        l.append(self.unAlnFastaFile)
        l.append(self.unAlnFastaFilePathName)
        l.append(self.alnFastaFile)
        l.append(self.alnFastaFilePathName)
        l.append(self.newickFile)
        l.append(self.newickFilePathName)
        l.append(None if self.speciesTree is None else self.speciesTree.id)
        l.append(self.cmdLine)
        l.append(self.softwareVersion)
        cursor.execute(sqlCmd, tuple(l))


class GeneTreeDataRelease(object):

    def __init__(self, geneTree=None, dataRelease=None):
        self.id = None
        self.geneTree = geneTree
        self.dataRelease = dataRelease
        # one-to-many

    def insertIntoDatabase(self, cursor):
        sqlCmd = 'INSERT INTO `GeneTreeDataRelease` (`geneTreeId`, `dataReleaseId`) VALUES (%s, %s)'
        l = []
        l.append(None if self.geneTree is None else self.geneTree.id)
        l.append(None if self.dataRelease is None else self.dataRelease.idDataRelease)
        cursor.execute(sqlCmd, tuple(l))


class GeneType(object):

    def __init__(self, geneTypeName=None):
        self.id = None
        self.geneTypeName = geneTypeName
        # one-to-many
        # fk_PaftolGene_geneTypeId: PaftolGene.geneTypeId REFERENCES GeneType(geneTypeId)
        self.paftolGeneGeneTypeList = []

    def insertIntoDatabase(self, cursor):
        sqlCmd = 'INSERT INTO `GeneType` (`geneTypeName`) VALUES (%s)'
        l = []
        l.append(self.geneTypeName)
        cursor.execute(sqlCmd, tuple(l))


class InputSequence(object):

    def __init__(self, dataOrigin=None, sequenceType=None, filename=None, pathName=None, md5sum=None, fastqStats=None, paftolSequence=None, sraRunSequence=None, OneKP_Sequence=None, annotatedGenome=None, GAP_Sequence=None, UnannotatedGenome=None):
        self.id = None
        self.dataOrigin = dataOrigin
        self.sequenceType = sequenceType
        self.filename = filename
        self.pathName = pathName
        self.md5sum = md5sum
        self.fastqStats = fastqStats
        self.paftolSequence = paftolSequence
        self.sraRunSequence = sraRunSequence
        self.OneKP_Sequence = OneKP_Sequence
        self.annotatedGenome = annotatedGenome
        self.GAP_Sequence = GAP_Sequence
        self.UnannotatedGenome = UnannotatedGenome
        # one-to-many
        # fk_ContigRecovery_fwdFastqId: ContigRecovery.fwdFastqId REFERENCES InputSequence(fwdFastqId)
        self.contigRecoveryFwdFastqList = []
        # fk_ContigRecovery_revFastqId: ContigRecovery.revFastqId REFERENCES InputSequence(revFastqId)
        self.contigRecoveryRevFastqList = []

    def insertIntoDatabase(self, cursor):
        sqlCmd = 'INSERT INTO `InputSequence` (`dataOriginId`, `sequenceTypeId`, `filename`, `pathName`, `md5sum`, `fastqStatsId`, `paftolSequenceId`, `sraRunSequenceId`, `OneKP_SequenceId`, `annotatedGenomeId`, `GAP_SequenceId`, `UnannotatedGenomeId`) VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s)'
        l = []
        l.append(None if self.dataOrigin is None else self.dataOrigin.id)
        l.append(None if self.sequenceType is None else self.sequenceType.id)
        l.append(self.filename)
        l.append(self.pathName)
        l.append(self.md5sum)
        l.append(None if self.fastqStats is None else self.fastqStats.id)
        l.append(None if self.paftolSequence is None else self.paftolSequence.id)
        l.append(None if self.sraRunSequence is None else self.sraRunSequence.id)
        l.append(None if self.OneKP_Sequence is None else self.OneKP_Sequence.id)
        l.append(None if self.annotatedGenome is None else self.annotatedGenome.id)
        l.append(None if self.GAP_Sequence is None else self.GAP_Sequence.id)
        l.append(None if self.UnannotatedGenome is None else self.UnannotatedGenome.id)
        cursor.execute(sqlCmd, tuple(l))


class OneKP_Sequence(object):

    def __init__(self, idSequencing=None, sampleId=None):
        self.id = None
        self.idSequencing = idSequencing
        self.sampleId = sampleId
        # one-to-many
        # fk_InputSequence_OneKP_SequenceId: InputSequence.OneKP_SequenceId REFERENCES OneKP_Sequence(OneKP_SequenceId)
        self.inputSequenceOneKP_SequenceList = []
        # fk_OneKP_SequenceDataRelease_OneKP_SequenceId: OneKP_SequenceDataRelease.OneKP_SequenceId REFERENCES OneKP_Sequence(OneKP_SequenceId)
        self.oneKP_SequenceDataReleaseOneKP_SequenceList = []

    def insertIntoDatabase(self, cursor):
        sqlCmd = 'INSERT INTO `OneKP_Sequence` (`idSequencing`, `sampleId`) VALUES (%s, %s)'
        l = []
        l.append(self.idSequencing)
        l.append(self.sampleId)
        cursor.execute(sqlCmd, tuple(l))


class OneKP_SequenceDataRelease(object):

    def __init__(self, OneKP_Sequence=None, dataRelease=None):
        self.id = None
        self.OneKP_Sequence = OneKP_Sequence
        self.dataRelease = dataRelease
        # one-to-many

    def insertIntoDatabase(self, cursor):
        sqlCmd = 'INSERT INTO `OneKP_SequenceDataRelease` (`OneKP_SequenceId`, `dataReleaseId`) VALUES (%s, %s)'
        l = []
        l.append(None if self.OneKP_Sequence is None else self.OneKP_Sequence.id)
        l.append(None if self.dataRelease is None else self.dataRelease.idDataRelease)
        cursor.execute(sqlCmd, tuple(l))


class PaftolGene(object):

    def __init__(self, geneName=None, geneType=None, exemplarGene=None):
        self.id = None
        self.geneName = geneName
        self.geneType = geneType
        self.exemplarGene = exemplarGene
        # one-to-many
        # fk_GeneTree_paftolGeneId: GeneTree.paftolGeneId REFERENCES PaftolGene(paftolGeneId)
        self.geneTreePaftolGeneList = []
        # fk_RecoveredContig_paftolGeneId: RecoveredContig.paftolGeneId REFERENCES PaftolGene(paftolGeneId)
        self.recoveredContigPaftolGeneList = []
        # fk_ReferenceTarget_paftolGeneId: ReferenceTarget.paftolGeneId REFERENCES PaftolGene(paftolGeneId)
        self.referenceTargetPaftolGeneList = []

    def insertIntoDatabase(self, cursor):
        sqlCmd = 'INSERT INTO `PaftolGene` (`geneName`, `geneTypeId`, `exemplarGeneId`) VALUES (%s, %s, %s)'
        l = []
        l.append(self.geneName)
        l.append(None if self.geneType is None else self.geneType.id)
        l.append(None if self.exemplarGene is None else self.exemplarGene.id)
        cursor.execute(sqlCmd, tuple(l))


class PaftolSequence(object):

    def __init__(self, idSequencing=None, replicate=None):
        self.id = None
        self.idSequencing = idSequencing
        self.replicate = replicate
        # one-to-many
        # fk_InputSequence_paftolSequenceId: InputSequence.paftolSequenceId REFERENCES PaftolSequence(paftolSequenceId)
        self.inputSequencePaftolSequenceList = []

    def insertIntoDatabase(self, cursor):
        sqlCmd = 'INSERT INTO `PaftolSequence` (`idSequencing`, `replicateId`) VALUES (%s, %s)'
        l = []
        l.append(self.idSequencing)
        l.append(None if self.replicate is None else self.replicate.id)
        cursor.execute(sqlCmd, tuple(l))


class RecoveredContig(object):

    def __init__(self, contigRecovery=None, paftolGene=None, seqLength=None, md5sum=None, representativeReferenceTarget=None):
        self.id = None
        self.contigRecovery = contigRecovery
        self.paftolGene = paftolGene
        self.seqLength = seqLength
        self.md5sum = md5sum
        self.representativeReferenceTarget = representativeReferenceTarget
        # one-to-many

    def insertIntoDatabase(self, cursor):
        sqlCmd = 'INSERT INTO `RecoveredContig` (`contigRecoveryId`, `paftolGeneId`, `seqLength`, `md5sum`, `representativeReferenceTargetId`) VALUES (%s, %s, %s, %s, %s)'
        l = []
        l.append(None if self.contigRecovery is None else self.contigRecovery.id)
        l.append(None if self.paftolGene is None else self.paftolGene.id)
        l.append(self.seqLength)
        l.append(self.md5sum)
        l.append(None if self.representativeReferenceTarget is None else self.representativeReferenceTarget.id)
        cursor.execute(sqlCmd, tuple(l))


class ReferenceTarget(object):

    def __init__(self, paftolGene=None, paftolOrganism=None, paftolTargetLength=None, targetsFastaFile=None, targetsFastaFilePathName=None, numTargetSequences=None, md5sum=None):
        self.id = None
        self.paftolGene = paftolGene
        self.paftolOrganism = paftolOrganism
        self.paftolTargetLength = paftolTargetLength
        self.targetsFastaFile = targetsFastaFile
        self.targetsFastaFilePathName = targetsFastaFilePathName
        self.numTargetSequences = numTargetSequences
        self.md5sum = md5sum
        # one-to-many
        # fk_ContigRecovery_referenceTargetId: ContigRecovery.referenceTargetId REFERENCES ReferenceTarget(referenceTargetId)
        self.contigRecoveryReferenceTargetList = []
        # fk_RecoveredContig_representativeReferenceTargetId: RecoveredContig.representativeReferenceTargetId REFERENCES ReferenceTarget(representativeReferenceTargetId)
        self.recoveredContigRepresentativeReferenceTargetList = []

    def insertIntoDatabase(self, cursor):
        sqlCmd = 'INSERT INTO `ReferenceTarget` (`paftolGeneId`, `paftolOrganism`, `paftolTargetLength`, `targetsFastaFile`, `targetsFastaFilePathName`, `numTargetSequences`, `md5sum`) VALUES (%s, %s, %s, %s, %s, %s, %s)'
        l = []
        l.append(None if self.paftolGene is None else self.paftolGene.id)
        l.append(self.paftolOrganism)
        l.append(self.paftolTargetLength)
        l.append(self.targetsFastaFile)
        l.append(self.targetsFastaFilePathName)
        l.append(self.numTargetSequences)
        l.append(self.md5sum)
        cursor.execute(sqlCmd, tuple(l))


class ReplicateSequence(object):

    def __init__(self, replicateIdNumber=None):
        self.id = None
        self.replicateIdNumber = replicateIdNumber
        # one-to-many
        # fk_PaftolSequence_replicateId: PaftolSequence.replicateId REFERENCES ReplicateSequence(replicateId)
        self.paftolSequenceReplicateList = []
        # fk_SRA_RunSequence_replicateId: SRA_RunSequence.replicateId REFERENCES ReplicateSequence(replicateId)
        self.srA_RunSequenceReplicateList = []

    def insertIntoDatabase(self, cursor):
        sqlCmd = 'INSERT INTO `ReplicateSequence` (`replicateIdNumber`) VALUES (%s)'
        l = []
        l.append(self.replicateIdNumber)
        cursor.execute(sqlCmd, tuple(l))


class SRA_RunSequence(object):

    def __init__(self, idSequencing=None, accessionId=None, replicate=None, enaAccession=None):
        self.id = None
        self.idSequencing = idSequencing
        self.accessionId = accessionId
        self.replicate = replicate
        self.enaAccession = enaAccession
        # one-to-many
        # fk_InputSequence_sraRunSequenceId: InputSequence.sraRunSequenceId REFERENCES SRA_RunSequence(sraRunSequenceId)
        self.inputSequenceSraRunSequenceList = []
        # fk_SRA_RunSequenceDataRelease_sraRunSequenceId: SRA_RunSequenceDataRelease.sraRunSequenceId REFERENCES SRA_RunSequence(sraRunSequenceId)
        self.srA_RunSequenceDataReleaseSraRunSequenceList = []

    def insertIntoDatabase(self, cursor):
        sqlCmd = 'INSERT INTO `SRA_RunSequence` (`idSequencing`, `accessionId`, `replicateId`, `enaAccessionId`) VALUES (%s, %s, %s, %s)'
        l = []
        l.append(self.idSequencing)
        l.append(self.accessionId)
        l.append(None if self.replicate is None else self.replicate.id)
        l.append(None if self.enaAccession is None else self.enaAccession.id)
        cursor.execute(sqlCmd, tuple(l))


class SRA_RunSequenceDataRelease(object):

    def __init__(self, sraRunSequence=None, dataRelease=None):
        self.id = None
        self.sraRunSequence = sraRunSequence
        self.dataRelease = dataRelease
        # one-to-many

    def insertIntoDatabase(self, cursor):
        sqlCmd = 'INSERT INTO `SRA_RunSequenceDataRelease` (`sraRunSequenceId`, `dataReleaseId`) VALUES (%s, %s)'
        l = []
        l.append(None if self.sraRunSequence is None else self.sraRunSequence.id)
        l.append(None if self.dataRelease is None else self.dataRelease.idDataRelease)
        cursor.execute(sqlCmd, tuple(l))


class SequenceType(object):

    def __init__(self, sequenceType=None, acronym=None):
        self.id = None
        self.sequenceType = sequenceType
        self.acronym = acronym
        # one-to-many
        # fk_InputSequence_sequenceTypeId: InputSequence.sequenceTypeId REFERENCES SequenceType(sequenceTypeId)
        self.inputSequenceSequenceTypeList = []

    def insertIntoDatabase(self, cursor):
        sqlCmd = 'INSERT INTO `SequenceType` (`sequenceType`, `acronym`) VALUES (%s, %s)'
        l = []
        l.append(self.sequenceType)
        l.append(self.acronym)
        cursor.execute(sqlCmd, tuple(l))


class SpeciesTree(object):

    def __init__(self, cmdLine=None, softwareVersion=None, newickFile=None, newickFilePathName=None):
        self.id = None
        self.cmdLine = cmdLine
        self.softwareVersion = softwareVersion
        self.newickFile = newickFile
        self.newickFilePathName = newickFilePathName
        # one-to-many
        # fk_GeneTree_speciesTreeId: GeneTree.speciesTreeId REFERENCES SpeciesTree(speciesTreeId)
        self.geneTreeSpeciesTreeList = []

    def insertIntoDatabase(self, cursor):
        sqlCmd = 'INSERT INTO `SpeciesTree` (`cmdLine`, `softwareVersion`, `newickFile`, `newickFilePathName`) VALUES (%s, %s, %s, %s)'
        l = []
        l.append(self.cmdLine)
        l.append(self.softwareVersion)
        l.append(self.newickFile)
        l.append(self.newickFilePathName)
        cursor.execute(sqlCmd, tuple(l))


class UnannotatedGenome(object):

    def __init__(self, idSequencing=None, accessionId=None, speciesLatinName=None, commonName=None, source=None, genomeVersion=None):
        self.id = None
        self.idSequencing = idSequencing
        self.accessionId = accessionId
        self.speciesLatinName = speciesLatinName
        self.commonName = commonName
        self.source = source
        self.genomeVersion = genomeVersion
        # one-to-many
        # fk_unannotatedGenomeId: InputSequence.UnannotatedGenomeId REFERENCES UnannotatedGenome(UnannotatedGenomeId)
        self.inputSequenceUnannotatedGenomeList = []

    def insertIntoDatabase(self, cursor):
        sqlCmd = 'INSERT INTO `UnannotatedGenome` (`idSequencing`, `accessionId`, `speciesLatinName`, `commonName`, `source`, `genomeVersion`) VALUES (%s, %s, %s, %s, %s, %s)'
        l = []
        l.append(self.idSequencing)
        l.append(self.accessionId)
        l.append(self.speciesLatinName)
        l.append(self.commonName)
        l.append(self.source)
        l.append(self.genomeVersion)
        cursor.execute(sqlCmd, tuple(l))


def loadAnnotatedGenomeDict(connection, productionDatabase):
    cursor = connection.cursor()
    entityDict = {}
    sqlStatement = 'SELECT `id`, `idSequencing`, `accessionId`, `speciesLatinName`, `commonName`, `source`, `genomeVersion` FROM `AnnotatedGenome`'
    cursor.execute(sqlStatement)
    for row in cursor:
        entity = AnnotatedGenome()
        entity.id = paftol.database.intOrNone(row[0])
        entity.idSequencing = paftol.database.intOrNone(row[1])
        entity.accessionId = paftol.database.strOrNone(row[2])
        entity.speciesLatinName = paftol.database.strOrNone(row[3])
        entity.commonName = paftol.database.strOrNone(row[4])
        entity.source = paftol.database.strOrNone(row[5])
        entity.genomeVersion = paftol.database.strOrNone(row[6])
        entityDict[entity.id] = entity
    cursor.close()
    return entityDict


def loadContigRecoveryDict(connection, productionDatabase):
    cursor = connection.cursor()
    entityDict = {}
    sqlStatement = 'SELECT `id`, `fwdFastqId`, `revFastqId`, `fwdTrimmedFastqStatsId`, `revTrimmedFastqStatsId`, `contigFastaFileName`, `contigFastaFilePathName`, `contigFastaFileMd5sum`, `referenceTargetId`, `numMappedReads`, `numUnmappedReads`, `softwareVersion`, `cmdLine`, `numRecoveredContigsCheck` FROM `ContigRecovery`'
    cursor.execute(sqlStatement)
    for row in cursor:
        entity = ContigRecovery()
        entity.id = paftol.database.intOrNone(row[0])
        # many to one: fwdFastq
        entityId = paftol.database.intOrNone(row[1])
        if entityId is None:
            entity.fwdFastq = None
        elif entityId not in productionDatabase.inputSequenceDict:
            raise StandardError, 'no InputSequence entity with id = %d' % entityId
        else:
            entity.fwdFastq = productionDatabase.inputSequenceDict[entityId]
            # type: int, name: fwdFastqId, foreignTable: InputSequence, foreignColumn: id
            entity.fwdFastq.contigRecoveryFwdFastqList.append(entity)
        # many to one: revFastq
        entityId = paftol.database.intOrNone(row[2])
        if entityId is None:
            entity.revFastq = None
        elif entityId not in productionDatabase.inputSequenceDict:
            raise StandardError, 'no InputSequence entity with id = %d' % entityId
        else:
            entity.revFastq = productionDatabase.inputSequenceDict[entityId]
            # type: int, name: revFastqId, foreignTable: InputSequence, foreignColumn: id
            entity.revFastq.contigRecoveryRevFastqList.append(entity)
        # many to one: fwdTrimmedFastqStats
        entityId = paftol.database.intOrNone(row[3])
        if entityId is None:
            entity.fwdTrimmedFastqStats = None
        elif entityId not in productionDatabase.fastqStatsDict:
            raise StandardError, 'no FastqStats entity with id = %d' % entityId
        else:
            entity.fwdTrimmedFastqStats = productionDatabase.fastqStatsDict[entityId]
            # type: int, name: fwdTrimmedFastqStatsId, foreignTable: FastqStats, foreignColumn: id
            entity.fwdTrimmedFastqStats.contigRecoveryFwdTrimmedFastqStatsList.append(entity)
        # many to one: revTrimmedFastqStats
        entityId = paftol.database.intOrNone(row[4])
        if entityId is None:
            entity.revTrimmedFastqStats = None
        elif entityId not in productionDatabase.fastqStatsDict:
            raise StandardError, 'no FastqStats entity with id = %d' % entityId
        else:
            entity.revTrimmedFastqStats = productionDatabase.fastqStatsDict[entityId]
            # type: int, name: revTrimmedFastqStatsId, foreignTable: FastqStats, foreignColumn: id
            entity.revTrimmedFastqStats.contigRecoveryRevTrimmedFastqStatsList.append(entity)
        entity.contigFastaFileName = paftol.database.strOrNone(row[5])
        entity.contigFastaFilePathName = paftol.database.strOrNone(row[6])
        entity.contigFastaFileMd5sum = paftol.database.strOrNone(row[7])
        # many to one: referenceTarget
        entityId = paftol.database.intOrNone(row[8])
        if entityId is None:
            entity.referenceTarget = None
        elif entityId not in productionDatabase.referenceTargetDict:
            raise StandardError, 'no ReferenceTarget entity with id = %d' % entityId
        else:
            entity.referenceTarget = productionDatabase.referenceTargetDict[entityId]
            # type: int, name: referenceTargetId, foreignTable: ReferenceTarget, foreignColumn: id
            entity.referenceTarget.contigRecoveryReferenceTargetList.append(entity)
        entity.numMappedReads = paftol.database.intOrNone(row[9])
        entity.numUnmappedReads = paftol.database.intOrNone(row[10])
        entity.softwareVersion = paftol.database.strOrNone(row[11])
        entity.cmdLine = paftol.database.strOrNone(row[12])
        entity.numRecoveredContigsCheck = paftol.database.intOrNone(row[13])
        entityDict[entity.id] = entity
    cursor.close()
    return entityDict


def loadContigRecoveryDataReleaseDict(connection, productionDatabase):
    cursor = connection.cursor()
    entityDict = {}
    sqlStatement = 'SELECT `id`, `contigRecoveryId`, `dataReleaseId` FROM `ContigRecoveryDataRelease`'
    cursor.execute(sqlStatement)
    for row in cursor:
        entity = ContigRecoveryDataRelease()
        entity.id = paftol.database.intOrNone(row[0])
        # many to one: contigRecovery
        entityId = paftol.database.intOrNone(row[1])
        if entityId is None:
            entity.contigRecovery = None
        elif entityId not in productionDatabase.contigRecoveryDict:
            raise StandardError, 'no ContigRecovery entity with id = %d' % entityId
        else:
            entity.contigRecovery = productionDatabase.contigRecoveryDict[entityId]
            # type: int, name: contigRecoveryId, foreignTable: ContigRecovery, foreignColumn: id
            entity.contigRecovery.contigRecoveryDataReleaseContigRecoveryList.append(entity)
        # many to one: dataRelease
        entityId = paftol.database.intOrNone(row[2])
        if entityId is None:
            entity.dataRelease = None
        elif entityId not in productionDatabase.dataReleaseDict:
            raise StandardError, 'no DataRelease entity with idDataRelease = %d' % entityId
        else:
            entity.dataRelease = productionDatabase.dataReleaseDict[entityId]
            # type: int, name: dataReleaseId, foreignTable: DataRelease, foreignColumn: idDataRelease
            entity.dataRelease.contigRecoveryDataReleaseDataReleaseList.append(entity)
        entityDict[entity.id] = entity
    cursor.close()
    return entityDict


def loadDataOriginDict(connection, productionDatabase):
    cursor = connection.cursor()
    entityDict = {}
    sqlStatement = 'SELECT `id`, `dataOriginName`, `acronym` FROM `DataOrigin`'
    cursor.execute(sqlStatement)
    for row in cursor:
        entity = DataOrigin()
        entity.id = paftol.database.intOrNone(row[0])
        entity.dataOriginName = paftol.database.strOrNone(row[1])
        entity.acronym = paftol.database.strOrNone(row[2])
        entityDict[entity.id] = entity
    cursor.close()
    return entityDict


def loadDataReleaseDict(connection, productionDatabase):
    cursor = connection.cursor()
    entityDict = {}
    sqlStatement = 'SELECT `idDataRelease`, `ReleaseNumber`, `DataRelease` FROM `DataRelease`'
    cursor.execute(sqlStatement)
    for row in cursor:
        entity = DataRelease()
        entity.idDataRelease = paftol.database.intOrNone(row[0])
        entity.releaseNumber = paftol.database.floatOrNone(row[1])
        entity.dataRelease = paftol.database.strOrNone(row[2])
        entityDict[entity.idDataRelease] = entity
    cursor.close()
    return entityDict


def loadENA_AccessionDict(connection, productionDatabase):
    cursor = connection.cursor()
    entityDict = {}
    sqlStatement = 'SELECT `id`, `accessionId` FROM `ENA_Accession`'
    cursor.execute(sqlStatement)
    for row in cursor:
        entity = ENA_Accession()
        entity.id = paftol.database.intOrNone(row[0])
        entity.accessionId = paftol.database.strOrNone(row[1])
        entityDict[entity.id] = entity
    cursor.close()
    return entityDict


def loadExemplarGeneDict(connection, productionDatabase):
    cursor = connection.cursor()
    entityDict = {}
    sqlStatement = 'SELECT `id`, `AC`, `GN`, `DE`, `OS`, `URL` FROM `ExemplarGene`'
    cursor.execute(sqlStatement)
    for row in cursor:
        entity = ExemplarGene()
        entity.id = paftol.database.intOrNone(row[0])
        entity.ac = paftol.database.strOrNone(row[1])
        entity.gn = paftol.database.strOrNone(row[2])
        entity.de = paftol.database.strOrNone(row[3])
        entity.os = paftol.database.strOrNone(row[4])
        entity.url = paftol.database.strOrNone(row[5])
        entityDict[entity.id] = entity
    cursor.close()
    return entityDict


def loadFastqStatsDict(connection, productionDatabase):
    cursor = connection.cursor()
    entityDict = {}
    sqlStatement = 'SELECT `id`, `numReads`, `numbrRecords`, `qual28`, `meanA`, `meanC`, `meanG`, `meanT`, `stddevA`, `stddevC`, `stddevG`, `stddevT`, `meanN`, `stddevN`, `meanAdapterContent`, `maxAdapterContent`, `sumLengthOfSeqs` FROM `FastqStats`'
    cursor.execute(sqlStatement)
    for row in cursor:
        entity = FastqStats()
        entity.id = paftol.database.intOrNone(row[0])
        entity.numReads = paftol.database.intOrNone(row[1])
        entity.numbrRecords = paftol.database.intOrNone(row[2])
        entity.qual28 = paftol.database.intOrNone(row[3])
        entity.meanA = paftol.database.floatOrNone(row[4])
        entity.meanC = paftol.database.floatOrNone(row[5])
        entity.meanG = paftol.database.floatOrNone(row[6])
        entity.meanT = paftol.database.floatOrNone(row[7])
        entity.stddevA = paftol.database.floatOrNone(row[8])
        entity.stddevC = paftol.database.floatOrNone(row[9])
        entity.stddevG = paftol.database.floatOrNone(row[10])
        entity.stddevT = paftol.database.floatOrNone(row[11])
        entity.meanN = paftol.database.floatOrNone(row[12])
        entity.stddevN = paftol.database.floatOrNone(row[13])
        entity.meanAdapterContent = paftol.database.floatOrNone(row[14])
        entity.maxAdapterContent = paftol.database.floatOrNone(row[15])
        entity.sumLengthOfSeqs = paftol.database.intOrNone(row[16])
        entityDict[entity.id] = entity
    cursor.close()
    return entityDict


def loadGAP_SequenceDict(connection, productionDatabase):
    cursor = connection.cursor()
    entityDict = {}
    sqlStatement = 'SELECT `id`, `idSequencing`, `sampleId` FROM `GAP_Sequence`'
    cursor.execute(sqlStatement)
    for row in cursor:
        entity = GAP_Sequence()
        entity.id = paftol.database.intOrNone(row[0])
        entity.idSequencing = paftol.database.intOrNone(row[1])
        entity.sampleId = paftol.database.intOrNone(row[2])
        entityDict[entity.id] = entity
    cursor.close()
    return entityDict


def loadGeneTreeDict(connection, productionDatabase):
    cursor = connection.cursor()
    entityDict = {}
    sqlStatement = 'SELECT `id`, `paftolGeneId`, `unAlnFastaFile`, `unAlnFastaFilePathName`, `alnFastaFile`, `alnFastaFilePathName`, `newickFile`, `newickFilePathName`, `speciesTreeId`, `cmdLine`, `softwareVersion` FROM `GeneTree`'
    cursor.execute(sqlStatement)
    for row in cursor:
        entity = GeneTree()
        entity.id = paftol.database.intOrNone(row[0])
        # many to one: paftolGene
        entityId = paftol.database.intOrNone(row[1])
        if entityId is None:
            entity.paftolGene = None
        elif entityId not in productionDatabase.paftolGeneDict:
            raise StandardError, 'no PaftolGene entity with id = %d' % entityId
        else:
            entity.paftolGene = productionDatabase.paftolGeneDict[entityId]
            # type: int, name: paftolGeneId, foreignTable: PaftolGene, foreignColumn: id
            entity.paftolGene.geneTreePaftolGeneList.append(entity)
        entity.unAlnFastaFile = paftol.database.strOrNone(row[2])
        entity.unAlnFastaFilePathName = paftol.database.strOrNone(row[3])
        entity.alnFastaFile = paftol.database.strOrNone(row[4])
        entity.alnFastaFilePathName = paftol.database.strOrNone(row[5])
        entity.newickFile = paftol.database.strOrNone(row[6])
        entity.newickFilePathName = paftol.database.strOrNone(row[7])
        # many to one: speciesTree
        entityId = paftol.database.intOrNone(row[8])
        if entityId is None:
            entity.speciesTree = None
        elif entityId not in productionDatabase.speciesTreeDict:
            raise StandardError, 'no SpeciesTree entity with id = %d' % entityId
        else:
            entity.speciesTree = productionDatabase.speciesTreeDict[entityId]
            # type: int, name: speciesTreeId, foreignTable: SpeciesTree, foreignColumn: id
            entity.speciesTree.geneTreeSpeciesTreeList.append(entity)
        entity.cmdLine = paftol.database.strOrNone(row[9])
        entity.softwareVersion = paftol.database.strOrNone(row[10])
        entityDict[entity.id] = entity
    cursor.close()
    return entityDict


def loadGeneTreeDataReleaseDict(connection, productionDatabase):
    cursor = connection.cursor()
    entityDict = {}
    sqlStatement = 'SELECT `id`, `geneTreeId`, `dataReleaseId` FROM `GeneTreeDataRelease`'
    cursor.execute(sqlStatement)
    for row in cursor:
        entity = GeneTreeDataRelease()
        entity.id = paftol.database.intOrNone(row[0])
        # many to one: geneTree
        entityId = paftol.database.intOrNone(row[1])
        if entityId is None:
            entity.geneTree = None
        elif entityId not in productionDatabase.geneTreeDict:
            raise StandardError, 'no GeneTree entity with id = %d' % entityId
        else:
            entity.geneTree = productionDatabase.geneTreeDict[entityId]
            # type: int, name: geneTreeId, foreignTable: GeneTree, foreignColumn: id
            entity.geneTree.geneTreeDataReleaseGeneTreeList.append(entity)
        # many to one: dataRelease
        entityId = paftol.database.intOrNone(row[2])
        if entityId is None:
            entity.dataRelease = None
        elif entityId not in productionDatabase.dataReleaseDict:
            raise StandardError, 'no DataRelease entity with idDataRelease = %d' % entityId
        else:
            entity.dataRelease = productionDatabase.dataReleaseDict[entityId]
            # type: int, name: dataReleaseId, foreignTable: DataRelease, foreignColumn: idDataRelease
            entity.dataRelease.geneTreeDataReleaseDataReleaseList.append(entity)
        entityDict[entity.id] = entity
    cursor.close()
    return entityDict


def loadGeneTypeDict(connection, productionDatabase):
    cursor = connection.cursor()
    entityDict = {}
    sqlStatement = 'SELECT `id`, `geneTypeName` FROM `GeneType`'
    cursor.execute(sqlStatement)
    for row in cursor:
        entity = GeneType()
        entity.id = paftol.database.intOrNone(row[0])
        entity.geneTypeName = paftol.database.strOrNone(row[1])
        entityDict[entity.id] = entity
    cursor.close()
    return entityDict


def loadInputSequenceDict(connection, productionDatabase):
    cursor = connection.cursor()
    entityDict = {}
    sqlStatement = 'SELECT `id`, `dataOriginId`, `sequenceTypeId`, `filename`, `pathName`, `md5sum`, `fastqStatsId`, `paftolSequenceId`, `sraRunSequenceId`, `OneKP_SequenceId`, `annotatedGenomeId`, `GAP_SequenceId`, `UnannotatedGenomeId` FROM `InputSequence`'
    cursor.execute(sqlStatement)
    for row in cursor:
        entity = InputSequence()
        entity.id = paftol.database.intOrNone(row[0])
        # many to one: dataOrigin
        entityId = paftol.database.intOrNone(row[1])
        if entityId is None:
            entity.dataOrigin = None
        elif entityId not in productionDatabase.dataOriginDict:
            raise StandardError, 'no DataOrigin entity with id = %d' % entityId
        else:
            entity.dataOrigin = productionDatabase.dataOriginDict[entityId]
            # type: int, name: dataOriginId, foreignTable: DataOrigin, foreignColumn: id
            entity.dataOrigin.inputSequenceDataOriginList.append(entity)
        # many to one: sequenceType
        entityId = paftol.database.intOrNone(row[2])
        if entityId is None:
            entity.sequenceType = None
        elif entityId not in productionDatabase.sequenceTypeDict:
            raise StandardError, 'no SequenceType entity with id = %d' % entityId
        else:
            entity.sequenceType = productionDatabase.sequenceTypeDict[entityId]
            # type: int, name: sequenceTypeId, foreignTable: SequenceType, foreignColumn: id
            entity.sequenceType.inputSequenceSequenceTypeList.append(entity)
        entity.filename = paftol.database.strOrNone(row[3])
        entity.pathName = paftol.database.strOrNone(row[4])
        entity.md5sum = paftol.database.strOrNone(row[5])
        # many to one: fastqStats
        entityId = paftol.database.intOrNone(row[6])
        if entityId is None:
            entity.fastqStats = None
        elif entityId not in productionDatabase.fastqStatsDict:
            raise StandardError, 'no FastqStats entity with id = %d' % entityId
        else:
            entity.fastqStats = productionDatabase.fastqStatsDict[entityId]
            # type: int, name: fastqStatsId, foreignTable: FastqStats, foreignColumn: id
            entity.fastqStats.inputSequenceFastqStatsList.append(entity)
        # many to one: paftolSequence
        entityId = paftol.database.intOrNone(row[7])
        if entityId is None:
            entity.paftolSequence = None
        elif entityId not in productionDatabase.paftolSequenceDict:
            raise StandardError, 'no PaftolSequence entity with id = %d' % entityId
        else:
            entity.paftolSequence = productionDatabase.paftolSequenceDict[entityId]
            # type: int, name: paftolSequenceId, foreignTable: PaftolSequence, foreignColumn: id
            entity.paftolSequence.inputSequencePaftolSequenceList.append(entity)
        # many to one: sraRunSequence
        entityId = paftol.database.intOrNone(row[8])
        if entityId is None:
            entity.sraRunSequence = None
        elif entityId not in productionDatabase.sRA_RunSequenceDict:
            raise StandardError, 'no SRA_RunSequence entity with id = %d' % entityId
        else:
            entity.sraRunSequence = productionDatabase.sRA_RunSequenceDict[entityId]
            # type: int, name: sraRunSequenceId, foreignTable: SRA_RunSequence, foreignColumn: id
            entity.sraRunSequence.inputSequenceSraRunSequenceList.append(entity)
        # many to one: OneKP_Sequence
        entityId = paftol.database.intOrNone(row[9])
        if entityId is None:
            entity.OneKP_Sequence = None
        elif entityId not in productionDatabase.oneKP_SequenceDict:
            raise StandardError, 'no OneKP_Sequence entity with id = %d' % entityId
        else:
            entity.OneKP_Sequence = productionDatabase.oneKP_SequenceDict[entityId]
            # type: int, name: OneKP_SequenceId, foreignTable: OneKP_Sequence, foreignColumn: id
            entity.OneKP_Sequence.inputSequenceOneKP_SequenceList.append(entity)
        # many to one: annotatedGenome
        entityId = paftol.database.intOrNone(row[10])
        if entityId is None:
            entity.annotatedGenome = None
        elif entityId not in productionDatabase.annotatedGenomeDict:
            raise StandardError, 'no AnnotatedGenome entity with id = %d' % entityId
        else:
            entity.annotatedGenome = productionDatabase.annotatedGenomeDict[entityId]
            # type: int, name: annotatedGenomeId, foreignTable: AnnotatedGenome, foreignColumn: id
            entity.annotatedGenome.inputSequenceAnnotatedGenomeList.append(entity)
        # many to one: GAP_Sequence
        entityId = paftol.database.intOrNone(row[11])
        if entityId is None:
            entity.GAP_Sequence = None
        elif entityId not in productionDatabase.gAP_SequenceDict:
            raise StandardError, 'no GAP_Sequence entity with id = %d' % entityId
        else:
            entity.GAP_Sequence = productionDatabase.gAP_SequenceDict[entityId]
            # type: int, name: GAP_SequenceId, foreignTable: GAP_Sequence, foreignColumn: id
            entity.GAP_Sequence.inputSequenceGAP_SequenceList.append(entity)
        # many to one: UnannotatedGenome
        entityId = paftol.database.intOrNone(row[12])
        if entityId is None:
            entity.UnannotatedGenome = None
        elif entityId not in productionDatabase.unannotatedGenomeDict:
            raise StandardError, 'no UnannotatedGenome entity with id = %d' % entityId
        else:
            entity.UnannotatedGenome = productionDatabase.unannotatedGenomeDict[entityId]
            # type: int, name: UnannotatedGenomeId, foreignTable: UnannotatedGenome, foreignColumn: id
            entity.UnannotatedGenome.inputSequenceUnannotatedGenomeList.append(entity)
        entityDict[entity.id] = entity
    cursor.close()
    return entityDict


def loadOneKP_SequenceDict(connection, productionDatabase):
    cursor = connection.cursor()
    entityDict = {}
    sqlStatement = 'SELECT `id`, `idSequencing`, `sampleId` FROM `OneKP_Sequence`'
    cursor.execute(sqlStatement)
    for row in cursor:
        entity = OneKP_Sequence()
        entity.id = paftol.database.intOrNone(row[0])
        entity.idSequencing = paftol.database.intOrNone(row[1])
        entity.sampleId = paftol.database.strOrNone(row[2])
        entityDict[entity.id] = entity
    cursor.close()
    return entityDict


def loadOneKP_SequenceDataReleaseDict(connection, productionDatabase):
    cursor = connection.cursor()
    entityDict = {}
    sqlStatement = 'SELECT `id`, `OneKP_SequenceId`, `dataReleaseId` FROM `OneKP_SequenceDataRelease`'
    cursor.execute(sqlStatement)
    for row in cursor:
        entity = OneKP_SequenceDataRelease()
        entity.id = paftol.database.intOrNone(row[0])
        # many to one: OneKP_Sequence
        entityId = paftol.database.intOrNone(row[1])
        if entityId is None:
            entity.OneKP_Sequence = None
        elif entityId not in productionDatabase.oneKP_SequenceDict:
            raise StandardError, 'no OneKP_Sequence entity with id = %d' % entityId
        else:
            entity.OneKP_Sequence = productionDatabase.oneKP_SequenceDict[entityId]
            # type: int, name: OneKP_SequenceId, foreignTable: OneKP_Sequence, foreignColumn: id
            entity.OneKP_Sequence.oneKP_SequenceDataReleaseOneKP_SequenceList.append(entity)
        # many to one: dataRelease
        entityId = paftol.database.intOrNone(row[2])
        if entityId is None:
            entity.dataRelease = None
        elif entityId not in productionDatabase.dataReleaseDict:
            raise StandardError, 'no DataRelease entity with idDataRelease = %d' % entityId
        else:
            entity.dataRelease = productionDatabase.dataReleaseDict[entityId]
            # type: int, name: dataReleaseId, foreignTable: DataRelease, foreignColumn: idDataRelease
            entity.dataRelease.oneKP_SequenceDataReleaseDataReleaseList.append(entity)
        entityDict[entity.id] = entity
    cursor.close()
    return entityDict


def loadPaftolGeneDict(connection, productionDatabase):
    cursor = connection.cursor()
    entityDict = {}
    sqlStatement = 'SELECT `id`, `geneName`, `geneTypeId`, `exemplarGeneId` FROM `PaftolGene`'
    cursor.execute(sqlStatement)
    for row in cursor:
        entity = PaftolGene()
        entity.id = paftol.database.intOrNone(row[0])
        entity.geneName = paftol.database.strOrNone(row[1])
        # many to one: geneType
        entityId = paftol.database.intOrNone(row[2])
        if entityId is None:
            entity.geneType = None
        elif entityId not in productionDatabase.geneTypeDict:
            raise StandardError, 'no GeneType entity with id = %d' % entityId
        else:
            entity.geneType = productionDatabase.geneTypeDict[entityId]
            # type: int, name: geneTypeId, foreignTable: GeneType, foreignColumn: id
            entity.geneType.paftolGeneGeneTypeList.append(entity)
        # many to one: exemplarGene
        entityId = paftol.database.intOrNone(row[3])
        if entityId is None:
            entity.exemplarGene = None
        elif entityId not in productionDatabase.exemplarGeneDict:
            raise StandardError, 'no ExemplarGene entity with id = %d' % entityId
        else:
            entity.exemplarGene = productionDatabase.exemplarGeneDict[entityId]
            # type: int, name: exemplarGeneId, foreignTable: ExemplarGene, foreignColumn: id
            entity.exemplarGene.paftolGeneExemplarGeneList.append(entity)
        entityDict[entity.id] = entity
    cursor.close()
    return entityDict


def loadPaftolSequenceDict(connection, productionDatabase):
    cursor = connection.cursor()
    entityDict = {}
    sqlStatement = 'SELECT `id`, `idSequencing`, `replicateId` FROM `PaftolSequence`'
    cursor.execute(sqlStatement)
    for row in cursor:
        entity = PaftolSequence()
        entity.id = paftol.database.intOrNone(row[0])
        entity.idSequencing = paftol.database.intOrNone(row[1])
        # many to one: replicate
        entityId = paftol.database.intOrNone(row[2])
        if entityId is None:
            entity.replicate = None
        elif entityId not in productionDatabase.replicateSequenceDict:
            raise StandardError, 'no ReplicateSequence entity with id = %d' % entityId
        else:
            entity.replicate = productionDatabase.replicateSequenceDict[entityId]
            # type: int, name: replicateId, foreignTable: ReplicateSequence, foreignColumn: id
            entity.replicate.paftolSequenceReplicateList.append(entity)
        entityDict[entity.id] = entity
    cursor.close()
    return entityDict


def loadRecoveredContigDict(connection, productionDatabase):
    cursor = connection.cursor()
    entityDict = {}
    sqlStatement = 'SELECT `id`, `contigRecoveryId`, `paftolGeneId`, `seqLength`, `md5sum`, `representativeReferenceTargetId` FROM `RecoveredContig`'
    cursor.execute(sqlStatement)
    for row in cursor:
        entity = RecoveredContig()
        entity.id = paftol.database.intOrNone(row[0])
        # many to one: contigRecovery
        entityId = paftol.database.intOrNone(row[1])
        if entityId is None:
            entity.contigRecovery = None
        elif entityId not in productionDatabase.contigRecoveryDict:
            raise StandardError, 'no ContigRecovery entity with id = %d' % entityId
        else:
            entity.contigRecovery = productionDatabase.contigRecoveryDict[entityId]
            # type: int, name: contigRecoveryId, foreignTable: ContigRecovery, foreignColumn: id
            entity.contigRecovery.recoveredContigContigRecoveryList.append(entity)
        # many to one: paftolGene
        entityId = paftol.database.intOrNone(row[2])
        if entityId is None:
            entity.paftolGene = None
        elif entityId not in productionDatabase.paftolGeneDict:
            raise StandardError, 'no PaftolGene entity with id = %d' % entityId
        else:
            entity.paftolGene = productionDatabase.paftolGeneDict[entityId]
            # type: int, name: paftolGeneId, foreignTable: PaftolGene, foreignColumn: id
            entity.paftolGene.recoveredContigPaftolGeneList.append(entity)
        entity.seqLength = paftol.database.intOrNone(row[3])
        entity.md5sum = paftol.database.strOrNone(row[4])
        # many to one: representativeReferenceTarget
        entityId = paftol.database.intOrNone(row[5])
        if entityId is None:
            entity.representativeReferenceTarget = None
        elif entityId not in productionDatabase.referenceTargetDict:
            raise StandardError, 'no ReferenceTarget entity with id = %d' % entityId
        else:
            entity.representativeReferenceTarget = productionDatabase.referenceTargetDict[entityId]
            # type: int, name: representativeReferenceTargetId, foreignTable: ReferenceTarget, foreignColumn: id
            entity.representativeReferenceTarget.recoveredContigRepresentativeReferenceTargetList.append(entity)
        entityDict[entity.id] = entity
    cursor.close()
    return entityDict


def loadReferenceTargetDict(connection, productionDatabase):
    cursor = connection.cursor()
    entityDict = {}
    sqlStatement = 'SELECT `id`, `paftolGeneId`, `paftolOrganism`, `paftolTargetLength`, `targetsFastaFile`, `targetsFastaFilePathName`, `numTargetSequences`, `md5sum` FROM `ReferenceTarget`'
    cursor.execute(sqlStatement)
    for row in cursor:
        entity = ReferenceTarget()
        entity.id = paftol.database.intOrNone(row[0])
        # many to one: paftolGene
        entityId = paftol.database.intOrNone(row[1])
        if entityId is None:
            entity.paftolGene = None
        elif entityId not in productionDatabase.paftolGeneDict:
            raise StandardError, 'no PaftolGene entity with id = %d' % entityId
        else:
            entity.paftolGene = productionDatabase.paftolGeneDict[entityId]
            # type: int, name: paftolGeneId, foreignTable: PaftolGene, foreignColumn: id
            entity.paftolGene.referenceTargetPaftolGeneList.append(entity)
        entity.paftolOrganism = paftol.database.strOrNone(row[2])
        entity.paftolTargetLength = paftol.database.intOrNone(row[3])
        entity.targetsFastaFile = paftol.database.strOrNone(row[4])
        entity.targetsFastaFilePathName = paftol.database.strOrNone(row[5])
        entity.numTargetSequences = paftol.database.intOrNone(row[6])
        entity.md5sum = paftol.database.strOrNone(row[7])
        entityDict[entity.id] = entity
    cursor.close()
    return entityDict


def loadReplicateSequenceDict(connection, productionDatabase):
    cursor = connection.cursor()
    entityDict = {}
    sqlStatement = 'SELECT `id`, `replicateIdNumber` FROM `ReplicateSequence`'
    cursor.execute(sqlStatement)
    for row in cursor:
        entity = ReplicateSequence()
        entity.id = paftol.database.intOrNone(row[0])
        entity.replicateIdNumber = paftol.database.intOrNone(row[1])
        entityDict[entity.id] = entity
    cursor.close()
    return entityDict


def loadSRA_RunSequenceDict(connection, productionDatabase):
    cursor = connection.cursor()
    entityDict = {}
    sqlStatement = 'SELECT `id`, `idSequencing`, `accessionId`, `replicateId`, `enaAccessionId` FROM `SRA_RunSequence`'
    cursor.execute(sqlStatement)
    for row in cursor:
        entity = SRA_RunSequence()
        entity.id = paftol.database.intOrNone(row[0])
        entity.idSequencing = paftol.database.intOrNone(row[1])
        entity.accessionId = paftol.database.strOrNone(row[2])
        # many to one: replicate
        entityId = paftol.database.intOrNone(row[3])
        if entityId is None:
            entity.replicate = None
        elif entityId not in productionDatabase.replicateSequenceDict:
            raise StandardError, 'no ReplicateSequence entity with id = %d' % entityId
        else:
            entity.replicate = productionDatabase.replicateSequenceDict[entityId]
            # type: int, name: replicateId, foreignTable: ReplicateSequence, foreignColumn: id
            entity.replicate.srA_RunSequenceReplicateList.append(entity)
        # many to one: enaAccession
        entityId = paftol.database.intOrNone(row[4])
        if entityId is None:
            entity.enaAccession = None
        elif entityId not in productionDatabase.eNA_AccessionDict:
            raise StandardError, 'no ENA_Accession entity with id = %d' % entityId
        else:
            entity.enaAccession = productionDatabase.eNA_AccessionDict[entityId]
            # type: int, name: enaAccessionId, foreignTable: ENA_Accession, foreignColumn: id
            entity.enaAccession.srA_RunSequenceEnaAccessionList.append(entity)
        entityDict[entity.id] = entity
    cursor.close()
    return entityDict


def loadSRA_RunSequenceDataReleaseDict(connection, productionDatabase):
    cursor = connection.cursor()
    entityDict = {}
    sqlStatement = 'SELECT `id`, `sraRunSequenceId`, `dataReleaseId` FROM `SRA_RunSequenceDataRelease`'
    cursor.execute(sqlStatement)
    for row in cursor:
        entity = SRA_RunSequenceDataRelease()
        entity.id = paftol.database.intOrNone(row[0])
        # many to one: sraRunSequence
        entityId = paftol.database.intOrNone(row[1])
        if entityId is None:
            entity.sraRunSequence = None
        elif entityId not in productionDatabase.sRA_RunSequenceDict:
            raise StandardError, 'no SRA_RunSequence entity with id = %d' % entityId
        else:
            entity.sraRunSequence = productionDatabase.sRA_RunSequenceDict[entityId]
            # type: int, name: sraRunSequenceId, foreignTable: SRA_RunSequence, foreignColumn: id
            entity.sraRunSequence.srA_RunSequenceDataReleaseSraRunSequenceList.append(entity)
        # many to one: dataRelease
        entityId = paftol.database.intOrNone(row[2])
        if entityId is None:
            entity.dataRelease = None
        elif entityId not in productionDatabase.dataReleaseDict:
            raise StandardError, 'no DataRelease entity with idDataRelease = %d' % entityId
        else:
            entity.dataRelease = productionDatabase.dataReleaseDict[entityId]
            # type: int, name: dataReleaseId, foreignTable: DataRelease, foreignColumn: idDataRelease
            entity.dataRelease.srA_RunSequenceDataReleaseDataReleaseList.append(entity)
        entityDict[entity.id] = entity
    cursor.close()
    return entityDict


def loadSequenceTypeDict(connection, productionDatabase):
    cursor = connection.cursor()
    entityDict = {}
    sqlStatement = 'SELECT `id`, `sequenceType`, `acronym` FROM `SequenceType`'
    cursor.execute(sqlStatement)
    for row in cursor:
        entity = SequenceType()
        entity.id = paftol.database.intOrNone(row[0])
        entity.sequenceType = paftol.database.strOrNone(row[1])
        entity.acronym = paftol.database.strOrNone(row[2])
        entityDict[entity.id] = entity
    cursor.close()
    return entityDict


def loadSpeciesTreeDict(connection, productionDatabase):
    cursor = connection.cursor()
    entityDict = {}
    sqlStatement = 'SELECT `id`, `cmdLine`, `softwareVersion`, `newickFile`, `newickFilePathName` FROM `SpeciesTree`'
    cursor.execute(sqlStatement)
    for row in cursor:
        entity = SpeciesTree()
        entity.id = paftol.database.intOrNone(row[0])
        entity.cmdLine = paftol.database.strOrNone(row[1])
        entity.softwareVersion = paftol.database.strOrNone(row[2])
        entity.newickFile = paftol.database.strOrNone(row[3])
        entity.newickFilePathName = paftol.database.strOrNone(row[4])
        entityDict[entity.id] = entity
    cursor.close()
    return entityDict


def loadUnannotatedGenomeDict(connection, productionDatabase):
    cursor = connection.cursor()
    entityDict = {}
    sqlStatement = 'SELECT `id`, `idSequencing`, `accessionId`, `speciesLatinName`, `commonName`, `source`, `genomeVersion` FROM `UnannotatedGenome`'
    cursor.execute(sqlStatement)
    for row in cursor:
        entity = UnannotatedGenome()
        entity.id = paftol.database.intOrNone(row[0])
        entity.idSequencing = paftol.database.intOrNone(row[1])
        entity.accessionId = paftol.database.strOrNone(row[2])
        entity.speciesLatinName = paftol.database.strOrNone(row[3])
        entity.commonName = paftol.database.strOrNone(row[4])
        entity.source = paftol.database.strOrNone(row[5])
        entity.genomeVersion = paftol.database.strOrNone(row[6])
        entityDict[entity.id] = entity
    cursor.close()
    return entityDict


class AnalysisDatabase(object):

    def __init__(self, connection):
        self.annotatedGenomeDict = {}
        self.contigRecoveryDict = {}
        self.contigRecoveryDataReleaseDict = {}
        self.dataOriginDict = {}
        self.dataReleaseDict = {}
        self.eNA_AccessionDict = {}
        self.exemplarGeneDict = {}
        self.fastqStatsDict = {}
        self.gAP_SequenceDict = {}
        self.geneTreeDict = {}
        self.geneTreeDataReleaseDict = {}
        self.geneTypeDict = {}
        self.inputSequenceDict = {}
        self.oneKP_SequenceDict = {}
        self.oneKP_SequenceDataReleaseDict = {}
        self.paftolGeneDict = {}
        self.paftolSequenceDict = {}
        self.recoveredContigDict = {}
        self.referenceTargetDict = {}
        self.replicateSequenceDict = {}
        self.sRA_RunSequenceDict = {}
        self.sRA_RunSequenceDataReleaseDict = {}
        self.sequenceTypeDict = {}
        self.speciesTreeDict = {}
        self.unannotatedGenomeDict = {}
        self.annotatedGenomeDict = loadAnnotatedGenomeDict(connection, self)
        self.dataOriginDict = loadDataOriginDict(connection, self)
        self.sequenceTypeDict = loadSequenceTypeDict(connection, self)
        self.fastqStatsDict = loadFastqStatsDict(connection, self)
        self.replicateSequenceDict = loadReplicateSequenceDict(connection, self)
        self.paftolSequenceDict = loadPaftolSequenceDict(connection, self)
        self.eNA_AccessionDict = loadENA_AccessionDict(connection, self)
        self.sRA_RunSequenceDict = loadSRA_RunSequenceDict(connection, self)
        self.oneKP_SequenceDict = loadOneKP_SequenceDict(connection, self)
        self.gAP_SequenceDict = loadGAP_SequenceDict(connection, self)
        self.unannotatedGenomeDict = loadUnannotatedGenomeDict(connection, self)
        self.inputSequenceDict = loadInputSequenceDict(connection, self)
        self.geneTypeDict = loadGeneTypeDict(connection, self)
        self.exemplarGeneDict = loadExemplarGeneDict(connection, self)
        self.paftolGeneDict = loadPaftolGeneDict(connection, self)
        self.referenceTargetDict = loadReferenceTargetDict(connection, self)
        self.contigRecoveryDict = loadContigRecoveryDict(connection, self)
        self.dataReleaseDict = loadDataReleaseDict(connection, self)
        self.contigRecoveryDataReleaseDict = loadContigRecoveryDataReleaseDict(connection, self)
        self.speciesTreeDict = loadSpeciesTreeDict(connection, self)
        self.geneTreeDict = loadGeneTreeDict(connection, self)
        self.geneTreeDataReleaseDict = loadGeneTreeDataReleaseDict(connection, self)
        self.oneKP_SequenceDataReleaseDict = loadOneKP_SequenceDataReleaseDict(connection, self)
        self.recoveredContigDict = loadRecoveredContigDict(connection, self)
        self.sRA_RunSequenceDataReleaseDict = loadSRA_RunSequenceDataReleaseDict(connection, self)

    def __str__(self):
        s = ''
        s = s + 'annotatedGenome: %d\n' % len(self.annotatedGenomeDict)
        s = s + 'dataOrigin: %d\n' % len(self.dataOriginDict)
        s = s + 'sequenceType: %d\n' % len(self.sequenceTypeDict)
        s = s + 'fastqStats: %d\n' % len(self.fastqStatsDict)
        s = s + 'replicateSequence: %d\n' % len(self.replicateSequenceDict)
        s = s + 'paftolSequence: %d\n' % len(self.paftolSequenceDict)
        s = s + 'eNA_Accession: %d\n' % len(self.eNA_AccessionDict)
        s = s + 'sRA_RunSequence: %d\n' % len(self.sRA_RunSequenceDict)
        s = s + 'oneKP_Sequence: %d\n' % len(self.oneKP_SequenceDict)
        s = s + 'gAP_Sequence: %d\n' % len(self.gAP_SequenceDict)
        s = s + 'unannotatedGenome: %d\n' % len(self.unannotatedGenomeDict)
        s = s + 'inputSequence: %d\n' % len(self.inputSequenceDict)
        s = s + 'geneType: %d\n' % len(self.geneTypeDict)
        s = s + 'exemplarGene: %d\n' % len(self.exemplarGeneDict)
        s = s + 'paftolGene: %d\n' % len(self.paftolGeneDict)
        s = s + 'referenceTarget: %d\n' % len(self.referenceTargetDict)
        s = s + 'contigRecovery: %d\n' % len(self.contigRecoveryDict)
        s = s + 'dataRelease: %d\n' % len(self.dataReleaseDict)
        s = s + 'contigRecoveryDataRelease: %d\n' % len(self.contigRecoveryDataReleaseDict)
        s = s + 'speciesTree: %d\n' % len(self.speciesTreeDict)
        s = s + 'geneTree: %d\n' % len(self.geneTreeDict)
        s = s + 'geneTreeDataRelease: %d\n' % len(self.geneTreeDataReleaseDict)
        s = s + 'oneKP_SequenceDataRelease: %d\n' % len(self.oneKP_SequenceDataReleaseDict)
        s = s + 'recoveredContig: %d\n' % len(self.recoveredContigDict)
        s = s + 'sRA_RunSequenceDataRelease: %d\n' % len(self.sRA_RunSequenceDataReleaseDict)
        return s

