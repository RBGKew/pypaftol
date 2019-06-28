#!/usr/bin/env python
import sys
import getopt
import re
import unicodedata

import mysql.connector

import paftol.database


class ContigRecovery(object):

    def __init__(self, id=None, fwdFastq=None, revFastq=None, targetsFastaFile=None, numMappedReads=None, totNumUnmappedReads=None, cmdLine=None):
        self.id = id
        self.fwdFastq = fwdFastq
        self.revFastq = revFastq
        self.targetsFastaFile = targetsFastaFile
        self.numMappedReads = numMappedReads
        self.totNumUnmappedReads = totNumUnmappedReads
        self.cmdLine = cmdLine
        # one-to-many
        # fk_RecoveredContig_contigRecoveryId: RecoveredContig.contigRecoveryId REFERENCES ContigRecovery(contigRecoveryId)
        self.recoveredContigContigRecoveryList = []


class FastaFile(object):

    def __init__(self, id=None, filename=None, md5sum=None, description=None, numSequences=None):
        self.id = id
        self.filename = filename
        self.md5sum = md5sum
        self.description = description
        self.numSequences = numSequences
        # one-to-many
        # fk_ContigRecovery_targetsFastaFileId: ContigRecovery.targetsFastaFileId REFERENCES FastaFile(targetsFastaFileId)
        self.contigRecoveryTargetsFastaFileList = []
        # fk_ReferenceTarget_targetsFastaFileId: ReferenceTarget.targetsFastaFileId REFERENCES FastaFile(targetsFastaFileId)
        self.referenceTargetTargetsFastaFileList = []


class FastqFile(object):

    def __init__(self, id=None, filename=None, md5sum=None, enaAccession=None, numReads=None, qual28=None, description=None):
        self.id = id
        self.filename = filename
        self.md5sum = md5sum
        self.enaAccession = enaAccession
        self.numReads = numReads
        self.qual28 = qual28
        self.description = description
        # one-to-many
        # fk_ContigRecovery_fwdFastqId_FastqFile: ContigRecovery.fwdFastqId REFERENCES FastqFile(fwdFastqId)
        self.contigRecoveryFwdFastqList = []
        # fk_ContigRecovery_revFastqId_FastqFile: ContigRecovery.revFastqId REFERENCES FastqFile(revFastqId)
        self.contigRecoveryRevFastqList = []
        # fk_PaftolFastqFile_fastqFileId_FastqFile: PaftolFastqFile.fastqFileId REFERENCES FastqFile(fastqFileId)
        self.paftolFastqFileFastqFileList = []
        # fk_Trimming_rawFwdFastqId_FastqFile: Trimming.rawFwdFastqId REFERENCES FastqFile(rawFwdFastqId)
        self.trimmingRawFwdFastqList = []
        # fk_Trimming_rawRevFastqId_FastqFile: Trimming.rawRevFastqId REFERENCES FastqFile(rawRevFastqId)
        self.trimmingRawRevFastqList = []
        # fk_Trimming_trimmedFwdPairedFastqId_FastqFile: Trimming.trimmedFwdPairedFastqId REFERENCES FastqFile(trimmedFwdPairedFastqId)
        self.trimmingTrimmedFwdPairedFastqList = []
        # fk_Trimming_trimmedRevPairedFastqId_FastqFile: Trimming.trimmedRevPairedFastqId REFERENCES FastqFile(trimmedRevPairedFastqId)
        self.trimmingTrimmedRevPairedFastqList = []
        # fk_Trimming_trimmedFwdUnpairedFastqId_FastqFile: Trimming.trimmedFwdUnpairedFastqId REFERENCES FastqFile(trimmedFwdUnpairedFastqId)
        self.trimmingTrimmedFwdUnpairedFastqList = []
        # fk_Trimming_trimmedRevUnpairedFastqId_FastqFile: Trimming.trimmedRevUnpairedFastqId REFERENCES FastqFile(trimmedRevUnpairedFastqId)
        self.trimmingTrimmedRevUnpairedFastqList = []


class GeneType(object):

    def __init__(self, id=None, geneTypeName=None):
        self.id = id
        self.geneTypeName = geneTypeName
        # one-to-many
        # fk_PaftolGene_geneTypeId_GeneType_id: PaftolGene.geneTypeId REFERENCES GeneType(geneTypeId)
        self.paftolGeneGeneTypeList = []


class PaftolFastqFile(object):

    def __init__(self, id=None, idSequencing=None, fastqFile=None):
        self.id = id
        self.idSequencing = idSequencing
        self.fastqFile = fastqFile
        # one-to-many
        # fk_RecoveredContig_fwdPaftolFastqId: RecoveredContig.fwdPaftolFastqId REFERENCES PaftolFastqFile(fwdPaftolFastqId)
        self.recoveredContigFwdPaftolFastqList = []
        # fk_RecoveredContig_revPaftolFastqId: RecoveredContig.revPaftolFastqId REFERENCES PaftolFastqFile(revPaftolFastqId)
        self.recoveredContigRevPaftolFastqList = []


class PaftolGene(object):

    def __init__(self, id=None, geneName=None, geneType=None):
        self.id = id
        self.geneName = geneName
        self.geneType = geneType
        # one-to-many
        # fk_RecoveredContig_paftolGeneId: RecoveredContig.paftolGeneId REFERENCES PaftolGene(paftolGeneId)
        self.recoveredContigPaftolGeneList = []
        # fk_ReferenceTarget_paftolGeneId_PaftolGene: ReferenceTarget.paftolGeneId REFERENCES PaftolGene(paftolGeneId)
        self.referenceTargetPaftolGeneList = []


class RecoveredContig(object):

    def __init__(self, id=None, contigRecovery=None, paftolGene=None, seqLength=None, fwdPaftolFastq=None, revPaftolFastq=None, representativeReferenceTarget=None):
        self.id = id
        self.contigRecovery = contigRecovery
        self.paftolGene = paftolGene
        self.seqLength = seqLength
        self.fwdPaftolFastq = fwdPaftolFastq
        self.revPaftolFastq = revPaftolFastq
        self.representativeReferenceTarget = representativeReferenceTarget
        # one-to-many


class ReferenceTarget(object):

    def __init__(self, id=None, paftolGene=None, paftolOrganism=None, paftolTargetLength=None, targetsFastaFile=None):
        self.id = id
        self.paftolGene = paftolGene
        self.paftolOrganism = paftolOrganism
        self.paftolTargetLength = paftolTargetLength
        self.targetsFastaFile = targetsFastaFile
        # one-to-many
        # fk_RecoveredContig_representativeReferenceTargetId: RecoveredContig.representativeReferenceTargetId REFERENCES ReferenceTarget(representativeReferenceTargetId)
        self.recoveredContigRepresentativeReferenceTargetList = []


class Trimming(object):

    def __init__(self, id=None, rawFwdFastq=None, rawRevFastq=None, trimmedFwdPairedFastq=None, trimmedRevPairedFastq=None, trimmedFwdUnpairedFastq=None, trimmedRevUnpairedFastq=None, cmdLine=None):
        self.id = id
        self.rawFwdFastq = rawFwdFastq
        self.rawRevFastq = rawRevFastq
        self.trimmedFwdPairedFastq = trimmedFwdPairedFastq
        self.trimmedRevPairedFastq = trimmedRevPairedFastq
        self.trimmedFwdUnpairedFastq = trimmedFwdUnpairedFastq
        self.trimmedRevUnpairedFastq = trimmedRevUnpairedFastq
        self.cmdLine = cmdLine
        # one-to-many


def loadContigRecoveryDict(connection, productionDatabase):
    cursor = connection.cursor()
    entityDict = {}
    sqlStatement = 'SELECT `id`, `fwdFastqId`, `revFastqId`, `targetsFastaFileId`, `numMappedReads`, `totNumUnmappedReads`, `cmdLine` FROM `ContigRecovery`'
    cursor.execute(sqlStatement)
    for row in cursor:
        entity = ContigRecovery()
        entity.id = paftol.database.intOrNone(row[0])
        # many to one: fwdFastq
        entityId = paftol.database.intOrNone(row[1])
        if entityId is None:
            entity.fwdFastq = None
        elif entityId not in productionDatabase.fastqFileDict:
            raise StandardError, 'no FastqFile entity with id = %d' % entityId
        else:
            entity.fwdFastq = productionDatabase.fastqFileDict[entityId]
            # type: int, name: fwdFastqId, foreignTable: FastqFile, foreignColumn: id
            entity.fwdFastq.contigRecoveryFwdFastqList.append(entity)
        # many to one: revFastq
        entityId = paftol.database.intOrNone(row[2])
        if entityId is None:
            entity.revFastq = None
        elif entityId not in productionDatabase.fastqFileDict:
            raise StandardError, 'no FastqFile entity with id = %d' % entityId
        else:
            entity.revFastq = productionDatabase.fastqFileDict[entityId]
            # type: int, name: revFastqId, foreignTable: FastqFile, foreignColumn: id
            entity.revFastq.contigRecoveryRevFastqList.append(entity)
        # many to one: targetsFastaFile
        entityId = paftol.database.intOrNone(row[3])
        if entityId is None:
            entity.targetsFastaFile = None
        elif entityId not in productionDatabase.fastaFileDict:
            raise StandardError, 'no FastaFile entity with id = %d' % entityId
        else:
            entity.targetsFastaFile = productionDatabase.fastaFileDict[entityId]
            # type: int, name: targetsFastaFileId, foreignTable: FastaFile, foreignColumn: id
            entity.targetsFastaFile.contigRecoveryTargetsFastaFileList.append(entity)
        entity.numMappedReads = paftol.database.intOrNone(row[4])
        entity.totNumUnmappedReads = paftol.database.intOrNone(row[5])
        entity.cmdLine = paftol.database.strOrNone(row[6])
        entityDict[entity.id] = entity
    cursor.close()
    return entityDict


def loadFastaFileDict(connection, productionDatabase):
    cursor = connection.cursor()
    entityDict = {}
    sqlStatement = 'SELECT `id`, `filename`, `md5sum`, `description`, `numSequences` FROM `FastaFile`'
    cursor.execute(sqlStatement)
    for row in cursor:
        entity = FastaFile()
        entity.id = paftol.database.intOrNone(row[0])
        entity.filename = paftol.database.strOrNone(row[1])
        entity.md5sum = paftol.database.strOrNone(row[2])
        entity.description = paftol.database.strOrNone(row[3])
        entity.numSequences = paftol.database.intOrNone(row[4])
        entityDict[entity.id] = entity
    cursor.close()
    return entityDict


def loadFastqFileDict(connection, productionDatabase):
    cursor = connection.cursor()
    entityDict = {}
    sqlStatement = 'SELECT `id`, `filename`, `md5sum`, `enaAccession`, `numReads`, `qual28`, `description` FROM `FastqFile`'
    cursor.execute(sqlStatement)
    for row in cursor:
        entity = FastqFile()
        entity.id = paftol.database.intOrNone(row[0])
        entity.filename = paftol.database.strOrNone(row[1])
        entity.md5sum = paftol.database.strOrNone(row[2])
        entity.enaAccession = paftol.database.strOrNone(row[3])
        entity.numReads = paftol.database.intOrNone(row[4])
        entity.qual28 = paftol.database.floatOrNone(row[5])
        entity.description = paftol.database.strOrNone(row[6])
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


def loadPaftolFastqFileDict(connection, productionDatabase):
    cursor = connection.cursor()
    entityDict = {}
    sqlStatement = 'SELECT `id`, `idSequencing`, `fastqFileId` FROM `PaftolFastqFile`'
    cursor.execute(sqlStatement)
    for row in cursor:
        entity = PaftolFastqFile()
        entity.id = paftol.database.intOrNone(row[0])
        entity.idSequencing = paftol.database.intOrNone(row[1])
        # many to one: fastqFile
        entityId = paftol.database.intOrNone(row[2])
        if entityId is None:
            entity.fastqFile = None
        elif entityId not in productionDatabase.fastqFileDict:
            raise StandardError, 'no FastqFile entity with id = %d' % entityId
        else:
            entity.fastqFile = productionDatabase.fastqFileDict[entityId]
            # type: int, name: fastqFileId, foreignTable: FastqFile, foreignColumn: id
            entity.fastqFile.paftolFastqFileFastqFileList.append(entity)
        entityDict[entity.id] = entity
    cursor.close()
    return entityDict


def loadPaftolGeneDict(connection, productionDatabase):
    cursor = connection.cursor()
    entityDict = {}
    sqlStatement = 'SELECT `id`, `geneName`, `geneTypeId` FROM `PaftolGene`'
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
        entityDict[entity.id] = entity
    cursor.close()
    return entityDict


def loadRecoveredContigDict(connection, productionDatabase):
    cursor = connection.cursor()
    entityDict = {}
    sqlStatement = 'SELECT `id`, `contigRecoveryId`, `paftolGeneId`, `seqLength`, `fwdPaftolFastqId`, `revPaftolFastqId`, `representativeReferenceTargetId` FROM `RecoveredContig`'
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
        # many to one: fwdPaftolFastq
        entityId = paftol.database.intOrNone(row[4])
        if entityId is None:
            entity.fwdPaftolFastq = None
        elif entityId not in productionDatabase.paftolFastqFileDict:
            raise StandardError, 'no PaftolFastqFile entity with id = %d' % entityId
        else:
            entity.fwdPaftolFastq = productionDatabase.paftolFastqFileDict[entityId]
            # type: int, name: fwdPaftolFastqId, foreignTable: PaftolFastqFile, foreignColumn: id
            entity.fwdPaftolFastq.recoveredContigFwdPaftolFastqList.append(entity)
        # many to one: revPaftolFastq
        entityId = paftol.database.intOrNone(row[5])
        if entityId is None:
            entity.revPaftolFastq = None
        elif entityId not in productionDatabase.paftolFastqFileDict:
            raise StandardError, 'no PaftolFastqFile entity with id = %d' % entityId
        else:
            entity.revPaftolFastq = productionDatabase.paftolFastqFileDict[entityId]
            # type: int, name: revPaftolFastqId, foreignTable: PaftolFastqFile, foreignColumn: id
            entity.revPaftolFastq.recoveredContigRevPaftolFastqList.append(entity)
        # many to one: representativeReferenceTarget
        entityId = paftol.database.intOrNone(row[6])
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
    sqlStatement = 'SELECT `id`, `paftolGeneId`, `paftolOrganism`, `paftolTargetLength`, `targetsFastaFileId` FROM `ReferenceTarget`'
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
        # many to one: targetsFastaFile
        entityId = paftol.database.intOrNone(row[4])
        if entityId is None:
            entity.targetsFastaFile = None
        elif entityId not in productionDatabase.fastaFileDict:
            raise StandardError, 'no FastaFile entity with id = %d' % entityId
        else:
            entity.targetsFastaFile = productionDatabase.fastaFileDict[entityId]
            # type: int, name: targetsFastaFileId, foreignTable: FastaFile, foreignColumn: id
            entity.targetsFastaFile.referenceTargetTargetsFastaFileList.append(entity)
        entityDict[entity.id] = entity
    cursor.close()
    return entityDict


def loadTrimmingDict(connection, productionDatabase):
    cursor = connection.cursor()
    entityDict = {}
    sqlStatement = 'SELECT `id`, `rawFwdFastqId`, `rawRevFastqId`, `trimmedFwdPairedFastqId`, `trimmedRevPairedFastqId`, `trimmedFwdUnpairedFastqId`, `trimmedRevUnpairedFastqId`, `cmdLine` FROM `Trimming`'
    cursor.execute(sqlStatement)
    for row in cursor:
        entity = Trimming()
        entity.id = paftol.database.intOrNone(row[0])
        # many to one: rawFwdFastq
        entityId = paftol.database.intOrNone(row[1])
        if entityId is None:
            entity.rawFwdFastq = None
        elif entityId not in productionDatabase.fastqFileDict:
            raise StandardError, 'no FastqFile entity with id = %d' % entityId
        else:
            entity.rawFwdFastq = productionDatabase.fastqFileDict[entityId]
            # type: int, name: rawFwdFastqId, foreignTable: FastqFile, foreignColumn: id
            entity.rawFwdFastq.trimmingRawFwdFastqList.append(entity)
        # many to one: rawRevFastq
        entityId = paftol.database.intOrNone(row[2])
        if entityId is None:
            entity.rawRevFastq = None
        elif entityId not in productionDatabase.fastqFileDict:
            raise StandardError, 'no FastqFile entity with id = %d' % entityId
        else:
            entity.rawRevFastq = productionDatabase.fastqFileDict[entityId]
            # type: int, name: rawRevFastqId, foreignTable: FastqFile, foreignColumn: id
            entity.rawRevFastq.trimmingRawRevFastqList.append(entity)
        # many to one: trimmedFwdPairedFastq
        entityId = paftol.database.intOrNone(row[3])
        if entityId is None:
            entity.trimmedFwdPairedFastq = None
        elif entityId not in productionDatabase.fastqFileDict:
            raise StandardError, 'no FastqFile entity with id = %d' % entityId
        else:
            entity.trimmedFwdPairedFastq = productionDatabase.fastqFileDict[entityId]
            # type: int, name: trimmedFwdPairedFastqId, foreignTable: FastqFile, foreignColumn: id
            entity.trimmedFwdPairedFastq.trimmingTrimmedFwdPairedFastqList.append(entity)
        # many to one: trimmedRevPairedFastq
        entityId = paftol.database.intOrNone(row[4])
        if entityId is None:
            entity.trimmedRevPairedFastq = None
        elif entityId not in productionDatabase.fastqFileDict:
            raise StandardError, 'no FastqFile entity with id = %d' % entityId
        else:
            entity.trimmedRevPairedFastq = productionDatabase.fastqFileDict[entityId]
            # type: int, name: trimmedRevPairedFastqId, foreignTable: FastqFile, foreignColumn: id
            entity.trimmedRevPairedFastq.trimmingTrimmedRevPairedFastqList.append(entity)
        # many to one: trimmedFwdUnpairedFastq
        entityId = paftol.database.intOrNone(row[5])
        if entityId is None:
            entity.trimmedFwdUnpairedFastq = None
        elif entityId not in productionDatabase.fastqFileDict:
            raise StandardError, 'no FastqFile entity with id = %d' % entityId
        else:
            entity.trimmedFwdUnpairedFastq = productionDatabase.fastqFileDict[entityId]
            # type: int, name: trimmedFwdUnpairedFastqId, foreignTable: FastqFile, foreignColumn: id
            entity.trimmedFwdUnpairedFastq.trimmingTrimmedFwdUnpairedFastqList.append(entity)
        # many to one: trimmedRevUnpairedFastq
        entityId = paftol.database.intOrNone(row[6])
        if entityId is None:
            entity.trimmedRevUnpairedFastq = None
        elif entityId not in productionDatabase.fastqFileDict:
            raise StandardError, 'no FastqFile entity with id = %d' % entityId
        else:
            entity.trimmedRevUnpairedFastq = productionDatabase.fastqFileDict[entityId]
            # type: int, name: trimmedRevUnpairedFastqId, foreignTable: FastqFile, foreignColumn: id
            entity.trimmedRevUnpairedFastq.trimmingTrimmedRevUnpairedFastqList.append(entity)
        entity.cmdLine = paftol.database.strOrNone(row[7])
        entityDict[entity.id] = entity
    cursor.close()
    return entityDict


class AnalysisDatabase(object):

    def __init__(self, connection):
        self.contigRecoveryDict = {}
        self.fastaFileDict = {}
        self.fastqFileDict = {}
        self.geneTypeDict = {}
        self.paftolFastqFileDict = {}
        self.paftolGeneDict = {}
        self.recoveredContigDict = {}
        self.referenceTargetDict = {}
        self.trimmingDict = {}
        self.fastqFileDict = loadFastqFileDict(connection, self)
        self.fastaFileDict = loadFastaFileDict(connection, self)
        self.contigRecoveryDict = loadContigRecoveryDict(connection, self)
        self.geneTypeDict = loadGeneTypeDict(connection, self)
        self.paftolFastqFileDict = loadPaftolFastqFileDict(connection, self)
        self.paftolGeneDict = loadPaftolGeneDict(connection, self)
        self.referenceTargetDict = loadReferenceTargetDict(connection, self)
        self.recoveredContigDict = loadRecoveredContigDict(connection, self)
        self.trimmingDict = loadTrimmingDict(connection, self)

    def __str__(self):
        s = ''
        s = s + 'fastqFile: %d\n' % len(self.fastqFileDict)
        s = s + 'fastaFile: %d\n' % len(self.fastaFileDict)
        s = s + 'contigRecovery: %d\n' % len(self.contigRecoveryDict)
        s = s + 'geneType: %d\n' % len(self.geneTypeDict)
        s = s + 'paftolFastqFile: %d\n' % len(self.paftolFastqFileDict)
        s = s + 'paftolGene: %d\n' % len(self.paftolGeneDict)
        s = s + 'referenceTarget: %d\n' % len(self.referenceTargetDict)
        s = s + 'recoveredContig: %d\n' % len(self.recoveredContigDict)
        s = s + 'trimming: %d\n' % len(self.trimmingDict)
        return s

