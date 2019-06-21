#!/usr/bin/env python
import sys
import getopt
import re
import unicodedata

import mysql.connector

import paftol.database


class ContigRecovery(object):

    def __init__(self, id=None, fwdFastqId=None, revFastqId=None, targetsFastaFileId=None, numMappedReads=None, totNumUnmappedReads=None, cmdLine=None):
        self.id = id
        self.fwdFastqId = fwdFastqId
        self.revFastqId = revFastqId
        self.targetsFastaFileId = targetsFastaFileId
        self.numMappedReads = numMappedReads
        self.totNumUnmappedReads = totNumUnmappedReads
        self.cmdLine = cmdLine
        # one-to-many
        self.recoveredContigList = []


class FastaFile(object):

    def __init__(self, id=None, filename=None, md5sum=None, description=None, numSequences=None):
        self.id = id
        self.filename = filename
        self.md5sum = md5sum
        self.description = description
        self.numSequences = numSequences
        # one-to-many
        self.contigRecoveryList = []
        self.referenceTargetList = []


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
        self.contigRecoveryList = []
        self.contigRecoveryList = []
        self.paftolFastqFileList = []
        self.trimmingList = []
        self.trimmingList = []
        self.trimmingList = []
        self.trimmingList = []
        self.trimmingList = []
        self.trimmingList = []


class GeneType(object):

    def __init__(self, id=None, geneTypeName=None):
        self.id = id
        self.geneTypeName = geneTypeName
        # one-to-many
        self.paftolGeneList = []


class PaftolFastqFile(object):

    def __init__(self, id=None, idSequencing=None, fastqFileId=None):
        self.id = id
        self.idSequencing = idSequencing
        self.fastqFileId = fastqFileId
        # one-to-many
        self.recoveredContigList = []
        self.recoveredContigList = []


class PaftolGene(object):

    def __init__(self, id=None, geneName=None, geneTypeId=None):
        self.id = id
        self.geneName = geneName
        self.geneTypeId = geneTypeId
        # one-to-many
        self.recoveredContigList = []
        self.referenceTargetList = []


class RecoveredContig(object):

    def __init__(self, id=None, contigRecoveryId=None, paftolGeneId=None, seqLength=None, fwdPaftolFastqId=None, revPaftolFastqId=None, representativeReferenceTargetId=None):
        self.id = id
        self.contigRecoveryId = contigRecoveryId
        self.paftolGeneId = paftolGeneId
        self.seqLength = seqLength
        self.fwdPaftolFastqId = fwdPaftolFastqId
        self.revPaftolFastqId = revPaftolFastqId
        self.representativeReferenceTargetId = representativeReferenceTargetId
        # one-to-many


class ReferenceTarget(object):

    def __init__(self, id=None, paftolGeneId=None, paftolOrganism=None, paftolTargetLength=None, targetsFastaFileId=None):
        self.id = id
        self.paftolGeneId = paftolGeneId
        self.paftolOrganism = paftolOrganism
        self.paftolTargetLength = paftolTargetLength
        self.targetsFastaFileId = targetsFastaFileId
        # one-to-many
        self.recoveredContigList = []


class Trimming(object):

    def __init__(self, id=None, rawFwdFastqId=None, rawRevFastqId=None, trimmedFwdPairedFastqId=None, trimmedRevPairedFastqId=None, trimmedFwdUnpairedFastqId=None, trimmedRevUnpairedFastqId=None, cmdLine=None):
        self.id = id
        self.rawFwdFastqId = rawFwdFastqId
        self.rawRevFastqId = rawRevFastqId
        self.trimmedFwdPairedFastqId = trimmedFwdPairedFastqId
        self.trimmedRevPairedFastqId = trimmedRevPairedFastqId
        self.trimmedFwdUnpairedFastqId = trimmedFwdUnpairedFastqId
        self.trimmedRevUnpairedFastqId = trimmedRevUnpairedFastqId
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
        # many to one: fwdFastqId
        entityId = paftol.database.intOrNone(row[1])
        if entityId is None:
            entity.fwdFastqId = None
        elif entityId not in productionDatabase.fastqFileDict:
            raise StandardError, 'no FastqFile entity with id = %d' % entityId
        else:
            entity.fwdFastqId = productionDatabase.fastqFileDict[entityId]
            entity.fwdFastqId.contigRecoveryList.append(entity)
        # many to one: revFastqId
        entityId = paftol.database.intOrNone(row[2])
        if entityId is None:
            entity.revFastqId = None
        elif entityId not in productionDatabase.fastqFileDict:
            raise StandardError, 'no FastqFile entity with id = %d' % entityId
        else:
            entity.revFastqId = productionDatabase.fastqFileDict[entityId]
            entity.revFastqId.contigRecoveryList.append(entity)
        # many to one: targetsFastaFileId
        entityId = paftol.database.intOrNone(row[3])
        if entityId is None:
            entity.targetsFastaFileId = None
        elif entityId not in productionDatabase.fastaFileDict:
            raise StandardError, 'no FastaFile entity with id = %d' % entityId
        else:
            entity.targetsFastaFileId = productionDatabase.fastaFileDict[entityId]
            entity.targetsFastaFileId.contigRecoveryList.append(entity)
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
        # many to one: fastqFileId
        entityId = paftol.database.intOrNone(row[2])
        if entityId is None:
            entity.fastqFileId = None
        elif entityId not in productionDatabase.fastqFileDict:
            raise StandardError, 'no FastqFile entity with id = %d' % entityId
        else:
            entity.fastqFileId = productionDatabase.fastqFileDict[entityId]
            entity.fastqFileId.paftolFastqFileList.append(entity)
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
        # many to one: geneTypeId
        entityId = paftol.database.intOrNone(row[2])
        if entityId is None:
            entity.geneTypeId = None
        elif entityId not in productionDatabase.geneTypeDict:
            raise StandardError, 'no GeneType entity with id = %d' % entityId
        else:
            entity.geneTypeId = productionDatabase.geneTypeDict[entityId]
            entity.geneTypeId.paftolGeneList.append(entity)
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
        # many to one: contigRecoveryId
        entityId = paftol.database.intOrNone(row[1])
        if entityId is None:
            entity.contigRecoveryId = None
        elif entityId not in productionDatabase.contigRecoveryDict:
            raise StandardError, 'no ContigRecovery entity with id = %d' % entityId
        else:
            entity.contigRecoveryId = productionDatabase.contigRecoveryDict[entityId]
            entity.contigRecoveryId.recoveredContigList.append(entity)
        # many to one: paftolGeneId
        entityId = paftol.database.intOrNone(row[2])
        if entityId is None:
            entity.paftolGeneId = None
        elif entityId not in productionDatabase.paftolGeneDict:
            raise StandardError, 'no PaftolGene entity with id = %d' % entityId
        else:
            entity.paftolGeneId = productionDatabase.paftolGeneDict[entityId]
            entity.paftolGeneId.recoveredContigList.append(entity)
        entity.seqLength = paftol.database.intOrNone(row[3])
        # many to one: fwdPaftolFastqId
        entityId = paftol.database.intOrNone(row[4])
        if entityId is None:
            entity.fwdPaftolFastqId = None
        elif entityId not in productionDatabase.paftolFastqFileDict:
            raise StandardError, 'no PaftolFastqFile entity with id = %d' % entityId
        else:
            entity.fwdPaftolFastqId = productionDatabase.paftolFastqFileDict[entityId]
            entity.fwdPaftolFastqId.recoveredContigList.append(entity)
        # many to one: revPaftolFastqId
        entityId = paftol.database.intOrNone(row[5])
        if entityId is None:
            entity.revPaftolFastqId = None
        elif entityId not in productionDatabase.paftolFastqFileDict:
            raise StandardError, 'no PaftolFastqFile entity with id = %d' % entityId
        else:
            entity.revPaftolFastqId = productionDatabase.paftolFastqFileDict[entityId]
            entity.revPaftolFastqId.recoveredContigList.append(entity)
        # many to one: representativeReferenceTargetId
        entityId = paftol.database.intOrNone(row[6])
        if entityId is None:
            entity.representativeReferenceTargetId = None
        elif entityId not in productionDatabase.referenceTargetDict:
            raise StandardError, 'no ReferenceTarget entity with id = %d' % entityId
        else:
            entity.representativeReferenceTargetId = productionDatabase.referenceTargetDict[entityId]
            entity.representativeReferenceTargetId.recoveredContigList.append(entity)
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
        # many to one: paftolGeneId
        entityId = paftol.database.intOrNone(row[1])
        if entityId is None:
            entity.paftolGeneId = None
        elif entityId not in productionDatabase.paftolGeneDict:
            raise StandardError, 'no PaftolGene entity with id = %d' % entityId
        else:
            entity.paftolGeneId = productionDatabase.paftolGeneDict[entityId]
            entity.paftolGeneId.referenceTargetList.append(entity)
        entity.paftolOrganism = paftol.database.strOrNone(row[2])
        entity.paftolTargetLength = paftol.database.intOrNone(row[3])
        # many to one: targetsFastaFileId
        entityId = paftol.database.intOrNone(row[4])
        if entityId is None:
            entity.targetsFastaFileId = None
        elif entityId not in productionDatabase.fastaFileDict:
            raise StandardError, 'no FastaFile entity with id = %d' % entityId
        else:
            entity.targetsFastaFileId = productionDatabase.fastaFileDict[entityId]
            entity.targetsFastaFileId.referenceTargetList.append(entity)
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
        # many to one: rawFwdFastqId
        entityId = paftol.database.intOrNone(row[1])
        if entityId is None:
            entity.rawFwdFastqId = None
        elif entityId not in productionDatabase.fastqFileDict:
            raise StandardError, 'no FastqFile entity with id = %d' % entityId
        else:
            entity.rawFwdFastqId = productionDatabase.fastqFileDict[entityId]
            entity.rawFwdFastqId.trimmingList.append(entity)
        # many to one: rawRevFastqId
        entityId = paftol.database.intOrNone(row[2])
        if entityId is None:
            entity.rawRevFastqId = None
        elif entityId not in productionDatabase.fastqFileDict:
            raise StandardError, 'no FastqFile entity with id = %d' % entityId
        else:
            entity.rawRevFastqId = productionDatabase.fastqFileDict[entityId]
            entity.rawRevFastqId.trimmingList.append(entity)
        # many to one: trimmedFwdPairedFastqId
        entityId = paftol.database.intOrNone(row[3])
        if entityId is None:
            entity.trimmedFwdPairedFastqId = None
        elif entityId not in productionDatabase.fastqFileDict:
            raise StandardError, 'no FastqFile entity with id = %d' % entityId
        else:
            entity.trimmedFwdPairedFastqId = productionDatabase.fastqFileDict[entityId]
            entity.trimmedFwdPairedFastqId.trimmingList.append(entity)
        # many to one: trimmedRevPairedFastqId
        entityId = paftol.database.intOrNone(row[4])
        if entityId is None:
            entity.trimmedRevPairedFastqId = None
        elif entityId not in productionDatabase.fastqFileDict:
            raise StandardError, 'no FastqFile entity with id = %d' % entityId
        else:
            entity.trimmedRevPairedFastqId = productionDatabase.fastqFileDict[entityId]
            entity.trimmedRevPairedFastqId.trimmingList.append(entity)
        # many to one: trimmedFwdUnpairedFastqId
        entityId = paftol.database.intOrNone(row[5])
        if entityId is None:
            entity.trimmedFwdUnpairedFastqId = None
        elif entityId not in productionDatabase.fastqFileDict:
            raise StandardError, 'no FastqFile entity with id = %d' % entityId
        else:
            entity.trimmedFwdUnpairedFastqId = productionDatabase.fastqFileDict[entityId]
            entity.trimmedFwdUnpairedFastqId.trimmingList.append(entity)
        # many to one: trimmedRevUnpairedFastqId
        entityId = paftol.database.intOrNone(row[6])
        if entityId is None:
            entity.trimmedRevUnpairedFastqId = None
        elif entityId not in productionDatabase.fastqFileDict:
            raise StandardError, 'no FastqFile entity with id = %d' % entityId
        else:
            entity.trimmedRevUnpairedFastqId = productionDatabase.fastqFileDict[entityId]
            entity.trimmedRevUnpairedFastqId.trimmingList.append(entity)
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

