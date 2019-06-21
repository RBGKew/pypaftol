#!/usr/bin/env python
import sys
import getopt
import re
import unicodedata

import mysql.connector

import paftol.database


class Action(object):

    def __init__(self, idAction=None, action=None):
        self.idAction = idAction
        self.action = action
        # one-to-many
        self.sampleList = []


class Coordinates(object):

    def __init__(self, coordinate=None, sortOrder=None):
        self.coordinate = coordinate
        self.sortOrder = sortOrder
        # one-to-many
        self.libraryList = []


class DBVersion(object):

    def __init__(self, id=None, dbName=None, dbDescription=None, dbVersion=None):
        self.id = id
        self.dbName = dbName
        self.dbDescription = dbDescription
        self.dbVersion = dbVersion
        # one-to-many


class DNAVolume(object):

    def __init__(self, idDnaVolume=None, dnaVolume=None):
        self.idDnaVolume = idDnaVolume
        self.dnaVolume = dnaVolume
        # one-to-many
        self.sampleList = []


class ExtractionType(object):

    def __init__(self, idExtractionType=None, extractionType=None):
        self.idExtractionType = idExtractionType
        self.extractionType = extractionType
        # one-to-many
        self.sampleList = []


class Family(object):

    def __init__(self, idFamily=None, family=None, order=None):
        self.idFamily = idFamily
        self.family = family
        self.order = order
        # one-to-many
        self.genusList = []


class Genus(object):

    def __init__(self, idGenus=None, family=None, genus=None, source=None, status=None, acceptedId=None, subfamily=None, tribe=None, subtribe=None, description=None, ipniId=None):
        self.idGenus = idGenus
        self.family = family
        self.genus = genus
        self.source = source
        self.status = status
        self.acceptedId = acceptedId
        self.subfamily = subfamily
        self.tribe = tribe
        self.subtribe = subtribe
        self.description = description
        self.ipniId = ipniId
        # one-to-many
        self.projectLinkList = []
        self.specimenList = []


class Indexes(object):

    def __init__(self, indexes=None, indexNameFwd=None, seqFwdPlatform1=None, seqFwdPlatform2=None, indexNameRv=None, seqRv=None):
        self.indexes = indexes
        self.indexNameFwd = indexNameFwd
        self.seqFwdPlatform1 = seqFwdPlatform1
        self.seqFwdPlatform2 = seqFwdPlatform2
        self.indexNameRv = indexNameRv
        self.seqRv = seqRv
        # one-to-many
        self.libraryList = []


class Library(object):

    def __init__(self, idLibrary=None, sample=None, libConcentration=None, hybridisationPool=None, libQuality=None, remainingVolume=None, libTapeStation=None, sonication=None, plate_No=None, coordinate=None, description=None, status=None, indexes=None, generateLibrary=None):
        self.idLibrary = idLibrary
        self.sample = sample
        self.libConcentration = libConcentration
        self.hybridisationPool = hybridisationPool
        self.libQuality = libQuality
        self.remainingVolume = remainingVolume
        self.libTapeStation = libTapeStation
        self.sonication = sonication
        self.plate_No = plate_No
        self.coordinate = coordinate
        self.description = description
        self.status = status
        self.indexes = indexes
        self.generateLibrary = generateLibrary
        # one-to-many
        self.sequenceList = []


class Order(object):

    def __init__(self, idOrder=None, order=None):
        self.idOrder = idOrder
        self.order = order
        # one-to-many
        self.familyList = []


class Project(object):

    def __init__(self, idProject=None, project=None):
        self.idProject = idProject
        self.project = project
        # one-to-many
        self.projectLinkList = []
        self.projectUserLinkList = []
        self.specimenList = []


class ProjectLink(object):

    def __init__(self, project=None, genus=None, notes=None, idProjectLink=None):
        self.project = project
        self.genus = genus
        self.notes = notes
        self.idProjectLink = idProjectLink
        # one-to-many


class ProjectUser(object):

    def __init__(self, idUser=None, name=None):
        self.idUser = idUser
        self.name = name
        # one-to-many
        self.projectUserLinkList = []


class ProjectUserLink(object):

    def __init__(self, user=None, project=None, idProjectUserLink=None):
        self.user = user
        self.project = project
        self.idProjectUserLink = idProjectUserLink
        # one-to-many


class Quality(object):

    def __init__(self, idQuality=None, quality=None):
        self.idQuality = idQuality
        self.quality = quality
        # one-to-many
        self.sampleList = []


class Sample(object):

    def __init__(self, idSample=None, specimen=None, description=None, action=None, extractionType=None, quality=None, sampleConcentration=None, newSampleConcentration=None, gelImage=None, sampleTapeStation=None, dnaVolume=None, compliant=None, materialSource=None, ageOfMaterial=None):
        self.idSample = idSample
        self.specimen = specimen
        self.description = description
        self.action = action
        self.extractionType = extractionType
        self.quality = quality
        self.sampleConcentration = sampleConcentration
        self.newSampleConcentration = newSampleConcentration
        self.gelImage = gelImage
        self.sampleTapeStation = sampleTapeStation
        self.dnaVolume = dnaVolume
        self.compliant = compliant
        self.materialSource = materialSource
        self.ageOfMaterial = ageOfMaterial
        # one-to-many
        self.libraryList = []


class Sequence(object):

    def __init__(self, idSequencing=None, library=None, platform=None, location=None, sequencingRun=None, numInferredCds=None, medianHybpiperCdsLength=None, status=None, hybridisationPool=None, r2FastqFile=None, r1FastqFile=None):
        self.idSequencing = idSequencing
        self.library = library
        self.platform = platform
        self.location = location
        self.sequencingRun = sequencingRun
        self.numInferredCds = numInferredCds
        self.medianHybpiperCdsLength = medianHybpiperCdsLength
        self.status = status
        self.hybridisationPool = hybridisationPool
        self.r2FastqFile = r2FastqFile
        self.r1FastqFile = r1FastqFile
        # one-to-many


class Source(object):

    def __init__(self, idSource=None, source=None):
        self.idSource = idSource
        self.source = source
        # one-to-many
        self.genusList = []


class SourceSpecimen(object):

    def __init__(self, idSourceSpecimen=None, sourceSpecimen=None):
        self.idSourceSpecimen = idSourceSpecimen
        self.sourceSpecimen = sourceSpecimen
        # one-to-many
        self.specimenList = []


class Species(object):

    def __init__(self, idSpecies=None, species=None, source=None):
        self.idSpecies = idSpecies
        self.species = species
        self.source = source
        # one-to-many
        self.specimenList = []


class Specimen(object):

    def __init__(self, idSpecimen=None, genus=None, idPaftol=None, species=None, description=None, bankId=None, lcd=None, msb=None, collector=None, collectorNo=None, voucherNo=None, museumId=None, museumBarcode=None, oldSpeciesName=None, sourceSpecimen=None, project=None):
        self.idSpecimen = idSpecimen
        self.genus = genus
        self.idPaftol = idPaftol
        self.species = species
        self.description = description
        self.bankId = bankId
        self.lcd = lcd
        self.msb = msb
        self.collector = collector
        self.collectorNo = collectorNo
        self.voucherNo = voucherNo
        self.museumId = museumId
        self.museumBarcode = museumBarcode
        self.oldSpeciesName = oldSpeciesName
        self.sourceSpecimen = sourceSpecimen
        self.project = project
        # one-to-many
        self.sampleList = []
        # no python attribute: SourceSpecimen
        # no python attribute: Project


class Status(object):

    def __init__(self, idStatus=None, status=None):
        self.idStatus = idStatus
        self.status = status
        # one-to-many
        self.libraryList = []
        self.sequenceList = []


def loadActionDict(connection, productionDatabase):
    cursor = connection.cursor()
    entityDict = {}
    sqlStatement = 'SELECT `idAction`, `Action` FROM `Action`'
    cursor.execute(sqlStatement)
    for row in cursor:
        entity = Action()
        entity.idAction = paftol.database.intOrNone(row[0])
        entity.action = paftol.database.strOrNone(row[1])
        entityDict[entity.idAction] = entity
    cursor.close()
    return entityDict


def loadCoordinatesDict(connection, productionDatabase):
    cursor = connection.cursor()
    entityDict = {}
    sqlStatement = 'SELECT `Coordinate`, `SortOrder` FROM `Coordinates`'
    cursor.execute(sqlStatement)
    for row in cursor:
        entity = Coordinates()
        entity.coordinate = paftol.database.strOrNone(row[0])
        entity.sortOrder = paftol.database.intOrNone(row[1])
        entityDict[entity.coordinate] = entity
    cursor.close()
    return entityDict


def loadDBVersionDict(connection, productionDatabase):
    cursor = connection.cursor()
    entityDict = {}
    sqlStatement = 'SELECT `ID`, `DBName`, `DBDescription`, `DBVersion` FROM `DBVersion`'
    cursor.execute(sqlStatement)
    for row in cursor:
        entity = DBVersion()
        entity.id = paftol.database.intOrNone(row[0])
        entity.dbName = paftol.database.strOrNone(row[1])
        entity.dbDescription = paftol.database.strOrNone(row[2])
        entity.dbVersion = paftol.database.strOrNone(row[3])
        entityDict[entity.id] = entity
    cursor.close()
    return entityDict


def loadDNAVolumeDict(connection, productionDatabase):
    cursor = connection.cursor()
    entityDict = {}
    sqlStatement = 'SELECT `idDNAVolume`, `DNAVolume` FROM `DNAVolume`'
    cursor.execute(sqlStatement)
    for row in cursor:
        entity = DNAVolume()
        entity.idDnaVolume = paftol.database.intOrNone(row[0])
        entity.dnaVolume = paftol.database.strOrNone(row[1])
        entityDict[entity.idDnaVolume] = entity
    cursor.close()
    return entityDict


def loadExtractionTypeDict(connection, productionDatabase):
    cursor = connection.cursor()
    entityDict = {}
    sqlStatement = 'SELECT `idExtractionType`, `ExtractionType` FROM `ExtractionType`'
    cursor.execute(sqlStatement)
    for row in cursor:
        entity = ExtractionType()
        entity.idExtractionType = paftol.database.intOrNone(row[0])
        entity.extractionType = paftol.database.strOrNone(row[1])
        entityDict[entity.idExtractionType] = entity
    cursor.close()
    return entityDict


def loadFamilyDict(connection, productionDatabase):
    cursor = connection.cursor()
    entityDict = {}
    sqlStatement = 'SELECT `idFamily`, `Family`, `idOrder` FROM `Family`'
    cursor.execute(sqlStatement)
    for row in cursor:
        entity = Family()
        entity.idFamily = paftol.database.intOrNone(row[0])
        entity.family = paftol.database.strOrNone(row[1])
        # many to one: order
        entityId = paftol.database.intOrNone(row[2])
        if entityId is None:
            entity.order = None
        elif entityId not in productionDatabase.orderDict:
            raise StandardError, 'no Order entity with idOrder = %d' % entityId
        else:
            entity.order = productionDatabase.orderDict[entityId]
            entity.order.familyList.append(entity)
        entityDict[entity.idFamily] = entity
    cursor.close()
    return entityDict


def loadGenusDict(connection, productionDatabase):
    cursor = connection.cursor()
    entityDict = {}
    sqlStatement = 'SELECT `idGenus`, `idFamily`, `Genus`, `idSource`, `Status`, `AcceptedId`, `Subfamily`, `Tribe`, `Subtribe`, `Description`, `IPNIid` FROM `Genus`'
    cursor.execute(sqlStatement)
    for row in cursor:
        entity = Genus()
        entity.idGenus = paftol.database.intOrNone(row[0])
        # many to one: family
        entityId = paftol.database.intOrNone(row[1])
        if entityId is None:
            entity.family = None
        elif entityId not in productionDatabase.familyDict:
            raise StandardError, 'no Family entity with idFamily = %d' % entityId
        else:
            entity.family = productionDatabase.familyDict[entityId]
            entity.family.genusList.append(entity)
        entity.genus = paftol.database.strOrNone(row[2])
        # many to one: source
        entityId = paftol.database.intOrNone(row[3])
        if entityId is None:
            entity.source = None
        elif entityId not in productionDatabase.sourceDict:
            raise StandardError, 'no Source entity with idSource = %d' % entityId
        else:
            entity.source = productionDatabase.sourceDict[entityId]
            entity.source.genusList.append(entity)
        entity.status = paftol.database.strOrNone(row[4])
        entity.acceptedId = paftol.database.intOrNone(row[5])
        entity.subfamily = paftol.database.strOrNone(row[6])
        entity.tribe = paftol.database.strOrNone(row[7])
        entity.subtribe = paftol.database.strOrNone(row[8])
        entity.description = paftol.database.strOrNone(row[9])
        entity.ipniId = paftol.database.strOrNone(row[10])
        entityDict[entity.idGenus] = entity
    cursor.close()
    return entityDict


def loadIndexesDict(connection, productionDatabase):
    cursor = connection.cursor()
    entityDict = {}
    sqlStatement = 'SELECT `Indexes`, `IndexNameFwd`, `SeqFwdPlatform1`, `SeqFwdPlatform2`, `IndexNameRv`, `SeqRv` FROM `Indexes`'
    cursor.execute(sqlStatement)
    for row in cursor:
        entity = Indexes()
        entity.indexes = paftol.database.strOrNone(row[0])
        entity.indexNameFwd = paftol.database.strOrNone(row[1])
        entity.seqFwdPlatform1 = paftol.database.strOrNone(row[2])
        entity.seqFwdPlatform2 = paftol.database.strOrNone(row[3])
        entity.indexNameRv = paftol.database.strOrNone(row[4])
        entity.seqRv = paftol.database.strOrNone(row[5])
        entityDict[entity.indexes] = entity
    cursor.close()
    return entityDict


def loadLibraryDict(connection, productionDatabase):
    cursor = connection.cursor()
    entityDict = {}
    sqlStatement = 'SELECT `idLibrary`, `idSample`, `LibConcentration`, `HybridisationPool`, `LibQuality`, `RemainingVolume`, `LibTapeStation`, `Sonication`, `Plate No`, `Coordinate`, `Description`, `idStatus`, `Indexes`, `GenerateLibrary` FROM `Library`'
    cursor.execute(sqlStatement)
    for row in cursor:
        entity = Library()
        entity.idLibrary = paftol.database.intOrNone(row[0])
        # many to one: sample
        entityId = paftol.database.intOrNone(row[1])
        if entityId is None:
            entity.sample = None
        elif entityId not in productionDatabase.sampleDict:
            raise StandardError, 'no Sample entity with idSample = %d' % entityId
        else:
            entity.sample = productionDatabase.sampleDict[entityId]
            entity.sample.libraryList.append(entity)
        entity.libConcentration = paftol.database.floatOrNone(row[2])
        entity.hybridisationPool = paftol.database.strOrNone(row[3])
        entity.libQuality = paftol.database.strOrNone(row[4])
        entity.remainingVolume = paftol.database.floatOrNone(row[5])
        entity.libTapeStation = paftol.database.strOrNone(row[6])
        entity.sonication = paftol.database.intOrNone(row[7])
        entity.plate_No = paftol.database.strOrNone(row[8])
        # many to one: coordinate
        entityId = paftol.database.strOrNone(row[9])
        if entityId is None:
            entity.coordinate = None
        elif entityId not in productionDatabase.coordinatesDict:
            raise StandardError, 'no Coordinates entity with coordinate = %d' % entityId
        else:
            entity.coordinate = productionDatabase.coordinatesDict[entityId]
            entity.coordinate.libraryList.append(entity)
        entity.description = paftol.database.strOrNone(row[10])
        # many to one: status
        entityId = paftol.database.intOrNone(row[11])
        if entityId is None:
            entity.status = None
        elif entityId not in productionDatabase.statusDict:
            raise StandardError, 'no Status entity with idStatus = %d' % entityId
        else:
            entity.status = productionDatabase.statusDict[entityId]
            entity.status.libraryList.append(entity)
        # many to one: indexes
        entityId = paftol.database.strOrNone(row[12])
        if entityId is None:
            entity.indexes = None
        elif entityId not in productionDatabase.indexesDict:
            raise StandardError, 'no Indexes entity with indexes = %d' % entityId
        else:
            entity.indexes = productionDatabase.indexesDict[entityId]
            entity.indexes.libraryList.append(entity)
        entity.generateLibrary = paftol.database.intOrNone(row[13])
        entityDict[entity.idLibrary] = entity
    cursor.close()
    return entityDict


def loadOrderDict(connection, productionDatabase):
    cursor = connection.cursor()
    entityDict = {}
    sqlStatement = 'SELECT `idOrder`, `Order` FROM `Order`'
    cursor.execute(sqlStatement)
    for row in cursor:
        entity = Order()
        entity.idOrder = paftol.database.intOrNone(row[0])
        entity.order = paftol.database.strOrNone(row[1])
        entityDict[entity.idOrder] = entity
    cursor.close()
    return entityDict


def loadProjectDict(connection, productionDatabase):
    cursor = connection.cursor()
    entityDict = {}
    sqlStatement = 'SELECT `idProject`, `Project` FROM `Project`'
    cursor.execute(sqlStatement)
    for row in cursor:
        entity = Project()
        entity.idProject = paftol.database.intOrNone(row[0])
        entity.project = paftol.database.strOrNone(row[1])
        entityDict[entity.idProject] = entity
    cursor.close()
    return entityDict


def loadProjectLinkDict(connection, productionDatabase):
    cursor = connection.cursor()
    entityDict = {}
    sqlStatement = 'SELECT `idProject`, `idGenus`, `Notes`, `idProjectLink` FROM `ProjectLink`'
    cursor.execute(sqlStatement)
    for row in cursor:
        entity = ProjectLink()
        # many to one: project
        entityId = paftol.database.intOrNone(row[0])
        if entityId is None:
            entity.project = None
        elif entityId not in productionDatabase.projectDict:
            raise StandardError, 'no Project entity with idProject = %d' % entityId
        else:
            entity.project = productionDatabase.projectDict[entityId]
            entity.project.projectLinkList.append(entity)
        # many to one: genus
        entityId = paftol.database.intOrNone(row[1])
        if entityId is None:
            entity.genus = None
        elif entityId not in productionDatabase.genusDict:
            raise StandardError, 'no Genus entity with idGenus = %d' % entityId
        else:
            entity.genus = productionDatabase.genusDict[entityId]
            entity.genus.projectLinkList.append(entity)
        entity.notes = paftol.database.strOrNone(row[2])
        entity.idProjectLink = paftol.database.intOrNone(row[3])
        entityDict[entity.idProjectLink] = entity
    cursor.close()
    return entityDict


def loadProjectUserDict(connection, productionDatabase):
    cursor = connection.cursor()
    entityDict = {}
    sqlStatement = 'SELECT `idUser`, `Name` FROM `ProjectUser`'
    cursor.execute(sqlStatement)
    for row in cursor:
        entity = ProjectUser()
        entity.idUser = paftol.database.strOrNone(row[0])
        entity.name = paftol.database.strOrNone(row[1])
        entityDict[entity.idUser] = entity
    cursor.close()
    return entityDict


def loadProjectUserLinkDict(connection, productionDatabase):
    cursor = connection.cursor()
    entityDict = {}
    sqlStatement = 'SELECT `idUser`, `idProject`, `idProjectUserLink` FROM `ProjectUserLink`'
    cursor.execute(sqlStatement)
    for row in cursor:
        entity = ProjectUserLink()
        # many to one: user
        entityId = paftol.database.strOrNone(row[0])
        if entityId is None:
            entity.user = None
        elif entityId not in productionDatabase.projectUserDict:
            raise StandardError, 'no ProjectUser entity with idUser = %d' % entityId
        else:
            entity.user = productionDatabase.projectUserDict[entityId]
            entity.user.projectUserLinkList.append(entity)
        # many to one: project
        entityId = paftol.database.intOrNone(row[1])
        if entityId is None:
            entity.project = None
        elif entityId not in productionDatabase.projectDict:
            raise StandardError, 'no Project entity with idProject = %d' % entityId
        else:
            entity.project = productionDatabase.projectDict[entityId]
            entity.project.projectUserLinkList.append(entity)
        entity.idProjectUserLink = paftol.database.intOrNone(row[2])
        entityDict[entity.idProjectUserLink] = entity
    cursor.close()
    return entityDict


def loadQualityDict(connection, productionDatabase):
    cursor = connection.cursor()
    entityDict = {}
    sqlStatement = 'SELECT `idQuality`, `Quality` FROM `Quality`'
    cursor.execute(sqlStatement)
    for row in cursor:
        entity = Quality()
        entity.idQuality = paftol.database.intOrNone(row[0])
        entity.quality = paftol.database.strOrNone(row[1])
        entityDict[entity.idQuality] = entity
    cursor.close()
    return entityDict


def loadSampleDict(connection, productionDatabase):
    cursor = connection.cursor()
    entityDict = {}
    sqlStatement = 'SELECT `idSample`, `idSpecimen`, `Description`, `idAction`, `idExtractionType`, `idQuality`, `SampleConcentration`, `NewSampleConcentration`, `GelImage`, `SampleTapeStation`, `idDNAVolume`, `Compliant`, `MaterialSource`, `AgeOfMaterial` FROM `Sample`'
    cursor.execute(sqlStatement)
    for row in cursor:
        entity = Sample()
        entity.idSample = paftol.database.intOrNone(row[0])
        # many to one: specimen
        entityId = paftol.database.intOrNone(row[1])
        if entityId is None:
            entity.specimen = None
        elif entityId not in productionDatabase.specimenDict:
            raise StandardError, 'no Specimen entity with idSpecimen = %d' % entityId
        else:
            entity.specimen = productionDatabase.specimenDict[entityId]
            entity.specimen.sampleList.append(entity)
        entity.description = paftol.database.strOrNone(row[2])
        # many to one: action
        entityId = paftol.database.intOrNone(row[3])
        if entityId is None:
            entity.action = None
        elif entityId not in productionDatabase.actionDict:
            raise StandardError, 'no Action entity with idAction = %d' % entityId
        else:
            entity.action = productionDatabase.actionDict[entityId]
            entity.action.sampleList.append(entity)
        # many to one: extractionType
        entityId = paftol.database.intOrNone(row[4])
        if entityId is None:
            entity.extractionType = None
        elif entityId not in productionDatabase.extractionTypeDict:
            raise StandardError, 'no ExtractionType entity with idExtractionType = %d' % entityId
        else:
            entity.extractionType = productionDatabase.extractionTypeDict[entityId]
            entity.extractionType.sampleList.append(entity)
        # many to one: quality
        entityId = paftol.database.intOrNone(row[5])
        if entityId is None:
            entity.quality = None
        elif entityId not in productionDatabase.qualityDict:
            raise StandardError, 'no Quality entity with idQuality = %d' % entityId
        else:
            entity.quality = productionDatabase.qualityDict[entityId]
            entity.quality.sampleList.append(entity)
        entity.sampleConcentration = paftol.database.floatOrNone(row[6])
        entity.newSampleConcentration = paftol.database.floatOrNone(row[7])
        entity.gelImage = paftol.database.strOrNone(row[8])
        entity.sampleTapeStation = paftol.database.strOrNone(row[9])
        # many to one: dnaVolume
        entityId = paftol.database.intOrNone(row[10])
        if entityId is None:
            entity.dnaVolume = None
        elif entityId not in productionDatabase.dNAVolumeDict:
            raise StandardError, 'no DNAVolume entity with idDnaVolume = %d' % entityId
        else:
            entity.dnaVolume = productionDatabase.dNAVolumeDict[entityId]
            entity.dnaVolume.sampleList.append(entity)
        entity.compliant = paftol.database.strOrNone(row[11])
        entity.materialSource = paftol.database.strOrNone(row[12])
        entity.ageOfMaterial = paftol.database.intOrNone(row[13])
        entityDict[entity.idSample] = entity
    cursor.close()
    return entityDict


def loadSequenceDict(connection, productionDatabase):
    cursor = connection.cursor()
    entityDict = {}
    sqlStatement = 'SELECT `idSequencing`, `idLibrary`, `Platform`, `Location`, `SequencingRun`, `NumInferredCds`, `MedianHybpiperCdsLength`, `idStatus`, `HybridisationPool`, `R2FastqFile`, `R1FastqFile` FROM `Sequence`'
    cursor.execute(sqlStatement)
    for row in cursor:
        entity = Sequence()
        entity.idSequencing = paftol.database.intOrNone(row[0])
        # many to one: library
        entityId = paftol.database.intOrNone(row[1])
        if entityId is None:
            entity.library = None
        elif entityId not in productionDatabase.libraryDict:
            raise StandardError, 'no Library entity with idLibrary = %d' % entityId
        else:
            entity.library = productionDatabase.libraryDict[entityId]
            entity.library.sequenceList.append(entity)
        entity.platform = paftol.database.strOrNone(row[2])
        entity.location = paftol.database.strOrNone(row[3])
        entity.sequencingRun = paftol.database.strOrNone(row[4])
        entity.numInferredCds = paftol.database.intOrNone(row[5])
        entity.medianHybpiperCdsLength = paftol.database.floatOrNone(row[6])
        # many to one: status
        entityId = paftol.database.intOrNone(row[7])
        if entityId is None:
            entity.status = None
        elif entityId not in productionDatabase.statusDict:
            raise StandardError, 'no Status entity with idStatus = %d' % entityId
        else:
            entity.status = productionDatabase.statusDict[entityId]
            entity.status.sequenceList.append(entity)
        entity.hybridisationPool = paftol.database.strOrNone(row[8])
        entity.r2FastqFile = paftol.database.strOrNone(row[9])
        entity.r1FastqFile = paftol.database.strOrNone(row[10])
        entityDict[entity.idSequencing] = entity
    cursor.close()
    return entityDict


def loadSourceDict(connection, productionDatabase):
    cursor = connection.cursor()
    entityDict = {}
    sqlStatement = 'SELECT `idSource`, `Source` FROM `Source`'
    cursor.execute(sqlStatement)
    for row in cursor:
        entity = Source()
        entity.idSource = paftol.database.intOrNone(row[0])
        entity.source = paftol.database.strOrNone(row[1])
        entityDict[entity.idSource] = entity
    cursor.close()
    return entityDict


def loadSourceSpecimenDict(connection, productionDatabase):
    cursor = connection.cursor()
    entityDict = {}
    sqlStatement = 'SELECT `idSourceSpecimen`, `SourceSpecimen` FROM `SourceSpecimen`'
    cursor.execute(sqlStatement)
    for row in cursor:
        entity = SourceSpecimen()
        entity.idSourceSpecimen = paftol.database.intOrNone(row[0])
        entity.sourceSpecimen = paftol.database.strOrNone(row[1])
        entityDict[entity.idSourceSpecimen] = entity
    cursor.close()
    return entityDict


def loadSpeciesDict(connection, productionDatabase):
    cursor = connection.cursor()
    entityDict = {}
    sqlStatement = 'SELECT `idSpecies`, `Species`, `Source` FROM `Species`'
    cursor.execute(sqlStatement)
    for row in cursor:
        entity = Species()
        entity.idSpecies = paftol.database.intOrNone(row[0])
        entity.species = paftol.database.strOrNone(row[1])
        entity.source = paftol.database.strOrNone(row[2])
        entityDict[entity.idSpecies] = entity
    cursor.close()
    return entityDict


def loadSpecimenDict(connection, productionDatabase):
    cursor = connection.cursor()
    entityDict = {}
    sqlStatement = 'SELECT `idSpecimen`, `idGenus`, `idPaftol`, `idSpecies`, `Description`, `BankID`, `LCD`, `MSB`, `Collector`, `CollectorNo`, `VoucherNo`, `MuseumID`, `MuseumBarcode`, `OldSpeciesName`, `idSourceSpecimen`, `idProject` FROM `Specimen`'
    cursor.execute(sqlStatement)
    for row in cursor:
        entity = Specimen()
        entity.idSpecimen = paftol.database.intOrNone(row[0])
        # many to one: genus
        entityId = paftol.database.intOrNone(row[1])
        if entityId is None:
            entity.genus = None
        elif entityId not in productionDatabase.genusDict:
            raise StandardError, 'no Genus entity with idGenus = %d' % entityId
        else:
            entity.genus = productionDatabase.genusDict[entityId]
            entity.genus.specimenList.append(entity)
        entity.idPaftol = paftol.database.intOrNone(row[2])
        # many to one: species
        entityId = paftol.database.intOrNone(row[3])
        if entityId is None:
            entity.species = None
        elif entityId not in productionDatabase.speciesDict:
            raise StandardError, 'no Species entity with idSpecies = %d' % entityId
        else:
            entity.species = productionDatabase.speciesDict[entityId]
            entity.species.specimenList.append(entity)
        entity.description = paftol.database.strOrNone(row[4])
        entity.bankId = paftol.database.intOrNone(row[5])
        entity.lcd = paftol.database.strOrNone(row[6])
        entity.msb = paftol.database.intOrNone(row[7])
        entity.collector = paftol.database.strOrNone(row[8])
        entity.collectorNo = paftol.database.strOrNone(row[9])
        entity.voucherNo = paftol.database.strOrNone(row[10])
        entity.museumId = paftol.database.strOrNone(row[11])
        entity.museumBarcode = paftol.database.strOrNone(row[12])
        entity.oldSpeciesName = paftol.database.strOrNone(row[13])
        # many to one: sourceSpecimen
        entityId = paftol.database.intOrNone(row[14])
        if entityId is None:
            entity.sourceSpecimen = None
        elif entityId not in productionDatabase.sourceSpecimenDict:
            raise StandardError, 'no SourceSpecimen entity with idSourceSpecimen = %d' % entityId
        else:
            entity.sourceSpecimen = productionDatabase.sourceSpecimenDict[entityId]
            entity.sourceSpecimen.specimenList.append(entity)
        # many to one: project
        entityId = paftol.database.intOrNone(row[15])
        if entityId is None:
            entity.project = None
        elif entityId not in productionDatabase.projectDict:
            raise StandardError, 'no Project entity with idProject = %d' % entityId
        else:
            entity.project = productionDatabase.projectDict[entityId]
            entity.project.specimenList.append(entity)
        entityDict[entity.idSpecimen] = entity
    cursor.close()
    return entityDict


def loadStatusDict(connection, productionDatabase):
    cursor = connection.cursor()
    entityDict = {}
    sqlStatement = 'SELECT `idStatus`, `Status` FROM `Status`'
    cursor.execute(sqlStatement)
    for row in cursor:
        entity = Status()
        entity.idStatus = paftol.database.intOrNone(row[0])
        entity.status = paftol.database.strOrNone(row[1])
        entityDict[entity.idStatus] = entity
    cursor.close()
    return entityDict


class ProductionDatabase(object):

    def __init__(self, connection):
        self.actionDict = {}
        self.coordinatesDict = {}
        self.dBVersionDict = {}
        self.dNAVolumeDict = {}
        self.extractionTypeDict = {}
        self.familyDict = {}
        self.genusDict = {}
        self.indexesDict = {}
        self.libraryDict = {}
        self.orderDict = {}
        self.projectDict = {}
        self.projectLinkDict = {}
        self.projectUserDict = {}
        self.projectUserLinkDict = {}
        self.qualityDict = {}
        self.sampleDict = {}
        self.sequenceDict = {}
        self.sourceDict = {}
        self.sourceSpecimenDict = {}
        self.speciesDict = {}
        self.specimenDict = {}
        self.statusDict = {}
        self.actionDict = loadActionDict(connection, self)
        self.coordinatesDict = loadCoordinatesDict(connection, self)
        self.dBVersionDict = loadDBVersionDict(connection, self)
        self.dNAVolumeDict = loadDNAVolumeDict(connection, self)
        self.extractionTypeDict = loadExtractionTypeDict(connection, self)
        self.orderDict = loadOrderDict(connection, self)
        self.familyDict = loadFamilyDict(connection, self)
        self.sourceDict = loadSourceDict(connection, self)
        self.genusDict = loadGenusDict(connection, self)
        self.indexesDict = loadIndexesDict(connection, self)
        self.speciesDict = loadSpeciesDict(connection, self)
        self.sourceSpecimenDict = loadSourceSpecimenDict(connection, self)
        self.projectDict = loadProjectDict(connection, self)
        self.specimenDict = loadSpecimenDict(connection, self)
        self.qualityDict = loadQualityDict(connection, self)
        self.sampleDict = loadSampleDict(connection, self)
        self.statusDict = loadStatusDict(connection, self)
        self.libraryDict = loadLibraryDict(connection, self)
        self.projectLinkDict = loadProjectLinkDict(connection, self)
        self.projectUserDict = loadProjectUserDict(connection, self)
        self.projectUserLinkDict = loadProjectUserLinkDict(connection, self)
        self.sequenceDict = loadSequenceDict(connection, self)

    def __str__(self):
        s = ''
        s = s + 'action: %d\n' % len(self.actionDict)
        s = s + 'coordinates: %d\n' % len(self.coordinatesDict)
        s = s + 'dBVersion: %d\n' % len(self.dBVersionDict)
        s = s + 'dNAVolume: %d\n' % len(self.dNAVolumeDict)
        s = s + 'extractionType: %d\n' % len(self.extractionTypeDict)
        s = s + 'order: %d\n' % len(self.orderDict)
        s = s + 'family: %d\n' % len(self.familyDict)
        s = s + 'source: %d\n' % len(self.sourceDict)
        s = s + 'genus: %d\n' % len(self.genusDict)
        s = s + 'indexes: %d\n' % len(self.indexesDict)
        s = s + 'species: %d\n' % len(self.speciesDict)
        s = s + 'sourceSpecimen: %d\n' % len(self.sourceSpecimenDict)
        s = s + 'project: %d\n' % len(self.projectDict)
        s = s + 'specimen: %d\n' % len(self.specimenDict)
        s = s + 'quality: %d\n' % len(self.qualityDict)
        s = s + 'sample: %d\n' % len(self.sampleDict)
        s = s + 'status: %d\n' % len(self.statusDict)
        s = s + 'library: %d\n' % len(self.libraryDict)
        s = s + 'projectLink: %d\n' % len(self.projectLinkDict)
        s = s + 'projectUser: %d\n' % len(self.projectUserDict)
        s = s + 'projectUserLink: %d\n' % len(self.projectUserLinkDict)
        s = s + 'sequence: %d\n' % len(self.sequenceDict)
        return s

