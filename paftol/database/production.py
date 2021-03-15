#!/usr/bin/env python

# Copyright (c) 2020 The Board of Trustees of the Royal Botanic Gardens, Kew

import sys
import getopt
import re
import unicodedata

import mysql.connector

import paftol.database


class Action(object):

    def __init__(self, action=None):
        self.idAction = None
        self.action = action
        # one-to-many
        # ActionLookup: Sample.idAction REFERENCES Action(idAction)
        self.sampleActionList = []

    def insertIntoDatabase(self, cursor):
        sqlCmd = 'INSERT INTO `Action` (`Action`) VALUES (%s)'
        l = []
        l.append(self.action)
        cursor.execute(sqlCmd, tuple(l))


class BlacklistedReason(object):

    def __init__(self, blacklistedReason=None):
        self.idBlacklistedReason = None
        self.blacklistedReason = blacklistedReason
        # one-to-many
        # fk_Sequence_BlacklistedReason: Sequence.idBlacklistedReason REFERENCES BlacklistedReason(idBlacklistedReason)
        self.sequenceBlacklistedReasonList = []

    def insertIntoDatabase(self, cursor):
        sqlCmd = 'INSERT INTO `BlacklistedReason` (`BlacklistedReason`) VALUES (%s)'
        l = []
        l.append(self.blacklistedReason)
        cursor.execute(sqlCmd, tuple(l))


class Coordinates(object):

    def __init__(self, coordinate=None):
        self.idCoordinate = None
        self.coordinate = coordinate
        # one-to-many
        # fk_Library_Coordinates: Library.idCoordinate REFERENCES Coordinates(idCoordinate)
        self.libraryCoordinateList = []

    def insertIntoDatabase(self, cursor):
        sqlCmd = 'INSERT INTO `Coordinates` (`Coordinate`) VALUES (%s)'
        l = []
        l.append(self.coordinate)
        cursor.execute(sqlCmd, tuple(l))


class DBVersion(object):

    def __init__(self, dbName=None, dbDescription=None, dbVersion=None):
        self.id = None
        self.dbName = dbName
        self.dbDescription = dbDescription
        self.dbVersion = dbVersion
        # one-to-many

    def insertIntoDatabase(self, cursor):
        sqlCmd = 'INSERT INTO `DBVersion` (`DBName`, `DBDescription`, `DBVersion`) VALUES (%s, %s, %s)'
        l = []
        l.append(self.dbName)
        l.append(self.dbDescription)
        l.append(self.dbVersion)
        cursor.execute(sqlCmd, tuple(l))


class DNAVolume(object):

    def __init__(self, dnaVolume=None):
        self.idDnaVolume = None
        self.dnaVolume = dnaVolume
        # one-to-many
        # DNALookup: Sample.idDNAVolume REFERENCES DNAVolume(idDNAVolume)
        self.sampleDnaVolumeList = []

    def insertIntoDatabase(self, cursor):
        sqlCmd = 'INSERT INTO `DNAVolume` (`DNAVolume`) VALUES (%s)'
        l = []
        l.append(self.dnaVolume)
        cursor.execute(sqlCmd, tuple(l))


class DataRelease(object):

    def __init__(self, releaseNumber=None, dataRelease=None, taxonCount=None):
        self.idDataRelease = None
        self.releaseNumber = releaseNumber
        self.dataRelease = dataRelease
        self.taxonCount = taxonCount
        # one-to-many
        # fk_SequenceDataRelease_DataRelease: SequenceDataRelease.idDataRelease REFERENCES DataRelease(idDataRelease)
        self.sequenceDataReleaseDataReleaseList = []
        # fk_idDataRelease: SpecimenDataRelease.idDataRelease REFERENCES DataRelease(idDataRelease)
        self.specimenDataReleaseDataReleaseList = []

    def insertIntoDatabase(self, cursor):
        sqlCmd = 'INSERT INTO `DataRelease` (`ReleaseNumber`, `DataRelease`, `TaxonCount`) VALUES (%s, %s, %s)'
        l = []
        l.append(self.releaseNumber)
        l.append(self.dataRelease)
        l.append(self.taxonCount)
        cursor.execute(sqlCmd, tuple(l))


class DataSource(object):

    def __init__(self, dataSource=None):
        self.idDataSource = None
        self.dataSource = dataSource
        # one-to-many
        # fk_Project_DataSource: Project.idDataSource REFERENCES DataSource(idDataSource)
        self.projectDataSourceList = []

    def insertIntoDatabase(self, cursor):
        sqlCmd = 'INSERT INTO `DataSource` (`DataSource`) VALUES (%s)'
        l = []
        l.append(self.dataSource)
        cursor.execute(sqlCmd, tuple(l))


class ExtractionType(object):

    def __init__(self, extractionType=None):
        self.idExtractionType = None
        self.extractionType = extractionType
        # one-to-many
        # ExtractionLookup: Sample.idExtractionType REFERENCES ExtractionType(idExtractionType)
        self.sampleExtractionTypeList = []

    def insertIntoDatabase(self, cursor):
        sqlCmd = 'INSERT INTO `ExtractionType` (`ExtractionType`) VALUES (%s)'
        l = []
        l.append(self.extractionType)
        cursor.execute(sqlCmd, tuple(l))


class Family(object):

    def __init__(self, family=None, order=None):
        self.idFamily = None
        self.family = family
        self.order = order
        # one-to-many
        # FamilyLink: Genus.idFamily REFERENCES Family(idFamily)
        self.genusFamilyList = []

    def insertIntoDatabase(self, cursor):
        sqlCmd = 'INSERT INTO `Family` (`Family`, `idOrder`) VALUES (%s, %s)'
        l = []
        l.append(self.family)
        l.append(None if self.order is None else self.order.idOrder)
        cursor.execute(sqlCmd, tuple(l))


class GeneStats(object):

    def __init__(self, internalName=None, exemplarAccession=None, exemplarName=None, exemplarSpecies=None, exemplarHyperlink=None, newickFile=None, newickFilePathName=None, averageContigLength=None, depth=None, averageContigLengthPercentage=None, numSeq=None, numGenera=None, numSpecies=None):
        self.id = None
        self.internalName = internalName
        self.exemplarAccession = exemplarAccession
        self.exemplarName = exemplarName
        self.exemplarSpecies = exemplarSpecies
        self.exemplarHyperlink = exemplarHyperlink
        self.newickFile = newickFile
        self.newickFilePathName = newickFilePathName
        self.averageContigLength = averageContigLength
        self.depth = depth
        self.averageContigLengthPercentage = averageContigLengthPercentage
        self.numSeq = numSeq
        self.numGenera = numGenera
        self.numSpecies = numSpecies
        # one-to-many

    def insertIntoDatabase(self, cursor):
        sqlCmd = 'INSERT INTO `GeneStats` (`InternalName`, `ExemplarAccession`, `ExemplarName`, `ExemplarSpecies`, `ExemplarHyperlink`, `NewickFile`, `NewickFilePathName`, `AverageContigLength`, `Depth`, `AverageContigLengthPercentage`, `NumSeq`, `NumGenera`, `NumSpecies`) VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s)'
        l = []
        l.append(self.internalName)
        l.append(self.exemplarAccession)
        l.append(self.exemplarName)
        l.append(self.exemplarSpecies)
        l.append(self.exemplarHyperlink)
        l.append(self.newickFile)
        l.append(self.newickFilePathName)
        l.append(self.averageContigLength)
        l.append(self.depth)
        l.append(self.averageContigLengthPercentage)
        l.append(self.numSeq)
        l.append(self.numGenera)
        l.append(self.numSpecies)
        cursor.execute(sqlCmd, tuple(l))


class Genus(object):

    def __init__(self, family=None, genus=None, source=None, status=None, acceptedId=None, subfamily=None, tribe=None, subtribe=None, description=None, ipniId=None):
        self.idGenus = None
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
        # GenusLink: Specimen.idGenus REFERENCES Genus(idGenus)
        self.specimenGenusList = []

    def insertIntoDatabase(self, cursor):
        sqlCmd = 'INSERT INTO `Genus` (`idFamily`, `Genus`, `idSource`, `Status`, `AcceptedId`, `Subfamily`, `Tribe`, `Subtribe`, `Description`, `IPNIid`) VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s, %s)'
        l = []
        l.append(None if self.family is None else self.family.idFamily)
        l.append(self.genus)
        l.append(None if self.source is None else self.source.idSource)
        l.append(self.status)
        l.append(self.acceptedId)
        l.append(self.subfamily)
        l.append(self.tribe)
        l.append(self.subtribe)
        l.append(self.description)
        l.append(self.ipniId)
        cursor.execute(sqlCmd, tuple(l))


class Indexes(object):

    def __init__(self, indexes=None, indexNameFwd=None, seqFwdPlatform1=None, seqFwdPlatform2=None, indexNameRv=None, seqRv=None):
        self.idIndexes = None
        self.indexes = indexes
        self.indexNameFwd = indexNameFwd
        self.seqFwdPlatform1 = seqFwdPlatform1
        self.seqFwdPlatform2 = seqFwdPlatform2
        self.indexNameRv = indexNameRv
        self.seqRv = seqRv
        # one-to-many
        # fk_Library_Indexes: Library.idIndexes REFERENCES Indexes(idIndexes)
        self.libraryIndexesList = []

    def insertIntoDatabase(self, cursor):
        sqlCmd = 'INSERT INTO `Indexes` (`Indexes`, `IndexNameFwd`, `SeqFwdPlatform1`, `SeqFwdPlatform2`, `IndexNameRv`, `SeqRv`) VALUES (%s, %s, %s, %s, %s, %s)'
        l = []
        l.append(self.indexes)
        l.append(self.indexNameFwd)
        l.append(self.seqFwdPlatform1)
        l.append(self.seqFwdPlatform2)
        l.append(self.indexNameRv)
        l.append(self.seqRv)
        cursor.execute(sqlCmd, tuple(l))


class Library(object):

    def __init__(self, sample=None, libConcentration=None, libQuality=None, remainingVolume=None, libTapeStation=None, sonication=None, plate_No=None, coordinate=None, description=None, status=None, indexes=None, generateLibrary=None):
        self.idLibrary = None
        self.sample = sample
        self.libConcentration = libConcentration
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
        # LibraryLink: Sequence.idLibrary REFERENCES Library(idLibrary)
        self.sequenceLibraryList = []

    def insertIntoDatabase(self, cursor):
        sqlCmd = 'INSERT INTO `Library` (`idSample`, `LibConcentration`, `LibQuality`, `RemainingVolume`, `LibTapeStation`, `Sonication`, `Plate No`, `idCoordinate`, `Description`, `idStatus`, `idIndexes`, `GenerateLibrary`) VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s)'
        l = []
        l.append(None if self.sample is None else self.sample.idSample)
        l.append(self.libConcentration)
        l.append(self.libQuality)
        l.append(self.remainingVolume)
        l.append(self.libTapeStation)
        l.append(self.sonication)
        l.append(self.plate_No)
        l.append(None if self.coordinate is None else self.coordinate.idCoordinate)
        l.append(self.description)
        l.append(None if self.status is None else self.status.idStatus)
        l.append(None if self.indexes is None else self.indexes.idIndexes)
        l.append(self.generateLibrary)
        cursor.execute(sqlCmd, tuple(l))


class Location(object):

    def __init__(self, location=None):
        self.idLocation = None
        self.location = location
        # one-to-many
        # fk_Sequence_Location: Sequence.idLocation REFERENCES Location(idLocation)
        self.sequenceLocationList = []

    def insertIntoDatabase(self, cursor):
        sqlCmd = 'INSERT INTO `Location` (`Location`) VALUES (%s)'
        l = []
        l.append(self.location)
        cursor.execute(sqlCmd, tuple(l))


class Log(object):

    def __init__(self, transactionDate=None, user=None, file=None, transactionType=None):
        self.idLog = None
        self.transactionDate = transactionDate
        self.user = user
        self.file = file
        self.transactionType = transactionType
        # one-to-many

    def insertIntoDatabase(self, cursor):
        sqlCmd = 'INSERT INTO `Log` (`TransactionDate`, `User`, `File`, `TransactionType`) VALUES (%s, %s, %s, %s)'
        l = []
        l.append(self.transactionDate)
        l.append(self.user)
        l.append(self.file)
        l.append(self.transactionType)
        cursor.execute(sqlCmd, tuple(l))


class MaterialSource(object):

    def __init__(self, materialSource=None):
        self.idMaterialSource = None
        self.materialSource = materialSource
        # one-to-many
        # fk_Specimen_MaterialSource: Specimen.idMaterialSource REFERENCES MaterialSource(idMaterialSource)
        self.specimenMaterialSourceList = []

    def insertIntoDatabase(self, cursor):
        sqlCmd = 'INSERT INTO `MaterialSource` (`MaterialSource`) VALUES (%s)'
        l = []
        l.append(self.materialSource)
        cursor.execute(sqlCmd, tuple(l))


class MigrationsLog(object):

    def __init__(self, migrationName=None, migrationDate=None):
        self.idMigrationsLog = None
        self.migrationName = migrationName
        self.migrationDate = migrationDate
        # one-to-many

    def insertIntoDatabase(self, cursor):
        sqlCmd = 'INSERT INTO `MigrationsLog` (`MigrationName`, `MigrationDate`) VALUES (%s, %s)'
        l = []
        l.append(self.migrationName)
        l.append(self.migrationDate)
        cursor.execute(sqlCmd, tuple(l))


class Museum(object):

    def __init__(self, museumId=None, museumName=None):
        self.idMuseumId = None
        self.museumId = museumId
        self.museumName = museumName
        # one-to-many
        # fk_idMuseumID: Specimen.idMuseumID REFERENCES Museum(idMuseumID)
        self.specimenMuseumIdList = []

    def insertIntoDatabase(self, cursor):
        sqlCmd = 'INSERT INTO `Museum` (`MuseumID`, `MuseumName`) VALUES (%s, %s)'
        l = []
        l.append(self.museumId)
        l.append(self.museumName)
        cursor.execute(sqlCmd, tuple(l))


class Order(object):

    def __init__(self, order=None):
        self.idOrder = None
        self.order = order
        # one-to-many
        # OrderLink: Family.idOrder REFERENCES Order(idOrder)
        self.familyOrderList = []

    def insertIntoDatabase(self, cursor):
        sqlCmd = 'INSERT INTO `Order` (`Order`) VALUES (%s)'
        l = []
        l.append(self.order)
        cursor.execute(sqlCmd, tuple(l))


class Platform(object):

    def __init__(self, platform=None):
        self.idPlatform = None
        self.platform = platform
        # one-to-many
        # fk_Sequence_Platform: Sequence.idPlatform REFERENCES Platform(idPlatform)
        self.sequencePlatformList = []

    def insertIntoDatabase(self, cursor):
        sqlCmd = 'INSERT INTO `Platform` (`Platform`) VALUES (%s)'
        l = []
        l.append(self.platform)
        cursor.execute(sqlCmd, tuple(l))


class Project(object):

    def __init__(self, project=None, dataSource=None):
        self.idProject = None
        self.project = project
        self.dataSource = dataSource
        # one-to-many
        # fk_idProject: Specimen.idProject REFERENCES Project(idProject)
        self.specimenProjectList = []

    def insertIntoDatabase(self, cursor):
        sqlCmd = 'INSERT INTO `Project` (`Project`, `idDataSource`) VALUES (%s, %s)'
        l = []
        l.append(self.project)
        l.append(None if self.dataSource is None else self.dataSource.idDataSource)
        cursor.execute(sqlCmd, tuple(l))


class Quality(object):

    def __init__(self, quality=None):
        self.idQuality = None
        self.quality = quality
        # one-to-many
        # QualityLookup: Sample.idQuality REFERENCES Quality(idQuality)
        self.sampleQualityList = []

    def insertIntoDatabase(self, cursor):
        sqlCmd = 'INSERT INTO `Quality` (`Quality`) VALUES (%s)'
        l = []
        l.append(self.quality)
        cursor.execute(sqlCmd, tuple(l))


class Sample(object):

    def __init__(self, specimen=None, description=None, action=None, extractionType=None, quality=None, sampleConcentration=None, newSampleConcentration=None, gelImage=None, sampleTapeStation=None, dnaVolume=None, compliant=None, enaSampleNum=None, secEnaSampleNum=None):
        self.idSample = None
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
        self.enaSampleNum = enaSampleNum
        self.secEnaSampleNum = secEnaSampleNum
        # one-to-many
        # SampleLink: Library.idSample REFERENCES Sample(idSample)
        self.librarySampleList = []

    def insertIntoDatabase(self, cursor):
        sqlCmd = 'INSERT INTO `Sample` (`idSpecimen`, `Description`, `idAction`, `idExtractionType`, `idQuality`, `SampleConcentration`, `NewSampleConcentration`, `GelImage`, `SampleTapeStation`, `idDNAVolume`, `Compliant`, `ENASampleNum`, `SecENASampleNum`) VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s)'
        l = []
        l.append(None if self.specimen is None else self.specimen.idSpecimen)
        l.append(self.description)
        l.append(None if self.action is None else self.action.idAction)
        l.append(None if self.extractionType is None else self.extractionType.idExtractionType)
        l.append(None if self.quality is None else self.quality.idQuality)
        l.append(self.sampleConcentration)
        l.append(self.newSampleConcentration)
        l.append(self.gelImage)
        l.append(self.sampleTapeStation)
        l.append(None if self.dnaVolume is None else self.dnaVolume.idDnaVolume)
        l.append(self.compliant)
        l.append(self.enaSampleNum)
        l.append(self.secEnaSampleNum)
        cursor.execute(sqlCmd, tuple(l))


class Sequence(object):

    def __init__(self, sequenceId=None, library=None, platform=None, location=None, sequencingRun=None, numInferredCds=None, medianHybpiperCdsLength=None, status=None, hybridisationPool=None, r2FastqFile=None, r1FastqFile=None, blacklisted=None, blacklistedReason=None, sequencingStrategy=None, enaExpNumber=None, enaRunNumber=None, externalSequenceId=None):
        self.idSequencing = None
        self.sequenceId = sequenceId
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
        self.blacklisted = blacklisted
        self.blacklistedReason = blacklistedReason
        self.sequencingStrategy = sequencingStrategy
        self.enaExpNumber = enaExpNumber
        self.enaRunNumber = enaRunNumber
        self.externalSequenceId = externalSequenceId
        # one-to-many
        # fk_SequenceDataRelease_Sequence: SequenceDataRelease.idSequencing REFERENCES Sequence(idSequencing)
        self.sequenceDataReleaseSequencingList = []

    def insertIntoDatabase(self, cursor):
        sqlCmd = 'INSERT INTO `Sequence` (`SequenceID`, `idLibrary`, `idPlatform`, `idLocation`, `SequencingRun`, `NumInferredCds`, `MedianHybpiperCdsLength`, `idStatus`, `HybridisationPool`, `R2FastqFile`, `R1FastqFile`, `Blacklisted`, `idBlacklistedReason`, `idSequencingStrategy`, `ENAExpNumber`, `ENARunNumber`, `ExternalSequenceID`) VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s)'
        l = []
        l.append(self.sequenceId)
        l.append(None if self.library is None else self.library.idLibrary)
        l.append(None if self.platform is None else self.platform.idPlatform)
        l.append(None if self.location is None else self.location.idLocation)
        l.append(self.sequencingRun)
        l.append(self.numInferredCds)
        l.append(self.medianHybpiperCdsLength)
        l.append(None if self.status is None else self.status.idStatus)
        l.append(self.hybridisationPool)
        l.append(self.r2FastqFile)
        l.append(self.r1FastqFile)
        l.append(self.blacklisted)
        l.append(None if self.blacklistedReason is None else self.blacklistedReason.idBlacklistedReason)
        l.append(None if self.sequencingStrategy is None else self.sequencingStrategy.idSequencingStrategy)
        l.append(self.enaExpNumber)
        l.append(self.enaRunNumber)
        l.append(self.externalSequenceId)
        cursor.execute(sqlCmd, tuple(l))


class SequenceDataRelease(object):

    def __init__(self, sequencing=None, dataRelease=None):
        self.idSequenceDataRelease = None
        self.sequencing = sequencing
        self.dataRelease = dataRelease
        # one-to-many

    def insertIntoDatabase(self, cursor):
        sqlCmd = 'INSERT INTO `SequenceDataRelease` (`idSequencing`, `idDataRelease`) VALUES (%s, %s)'
        l = []
        l.append(None if self.sequencing is None else self.sequencing.idSequencing)
        l.append(None if self.dataRelease is None else self.dataRelease.idDataRelease)
        cursor.execute(sqlCmd, tuple(l))


class SequencingStrategy(object):

    def __init__(self, sequencingStrategy=None):
        self.idSequencingStrategy = None
        self.sequencingStrategy = sequencingStrategy
        # one-to-many
        # fk_Sequence_SequencingStrategy: Sequence.idSequencingStrategy REFERENCES SequencingStrategy(idSequencingStrategy)
        self.sequenceSequencingStrategyList = []

    def insertIntoDatabase(self, cursor):
        sqlCmd = 'INSERT INTO `SequencingStrategy` (`SequencingStrategy`) VALUES (%s)'
        l = []
        l.append(self.sequencingStrategy)
        cursor.execute(sqlCmd, tuple(l))


class Source(object):

    def __init__(self, source=None):
        self.idSource = None
        self.source = source
        # one-to-many
        # SourceLookup: Genus.idSource REFERENCES Source(idSource)
        self.genusSourceList = []

    def insertIntoDatabase(self, cursor):
        sqlCmd = 'INSERT INTO `Source` (`Source`) VALUES (%s)'
        l = []
        l.append(self.source)
        cursor.execute(sqlCmd, tuple(l))


class SourceSpecimen(object):

    def __init__(self, sourceSpecimen=None):
        self.idSourceSpecimen = None
        self.sourceSpecimen = sourceSpecimen
        # one-to-many
        # fk_SourceSpecimen: Specimen.idSourceSpecimen REFERENCES SourceSpecimen(idSourceSpecimen)
        self.specimenSourceSpecimenList = []

    def insertIntoDatabase(self, cursor):
        sqlCmd = 'INSERT INTO `SourceSpecimen` (`SourceSpecimen`) VALUES (%s)'
        l = []
        l.append(self.sourceSpecimen)
        cursor.execute(sqlCmd, tuple(l))


class Species(object):

    def __init__(self, species=None, source=None):
        self.idSpecies = None
        self.species = species
        self.source = source
        # one-to-many
        # SpeciesLookup: Specimen.idSpecies REFERENCES Species(idSpecies)
        self.specimenSpeciesList = []

    def insertIntoDatabase(self, cursor):
        sqlCmd = 'INSERT INTO `Species` (`Species`, `Source`) VALUES (%s, %s)'
        l = []
        l.append(self.species)
        l.append(self.source)
        cursor.execute(sqlCmd, tuple(l))


class Specimen(object):

    def __init__(self, genus=None, idPaftol=None, species=None, bankId=None, lcd=None, msb=None, collector=None, collectorNo=None, voucherNo=None, museumBarcode=None, oldSpeciesName=None, sourceSpecimen=None, project=None, idOriginCountry=None, materialSource=None, ageOfMaterial=None, museumId=None, specimenReference=None):
        self.idSpecimen = None
        self.genus = genus
        self.idPaftol = idPaftol
        self.species = species
        self.bankId = bankId
        self.lcd = lcd
        self.msb = msb
        self.collector = collector
        self.collectorNo = collectorNo
        self.voucherNo = voucherNo
        self.museumBarcode = museumBarcode
        self.oldSpeciesName = oldSpeciesName
        self.sourceSpecimen = sourceSpecimen
        self.project = project
        self.idOriginCountry = idOriginCountry
        self.materialSource = materialSource
        self.ageOfMaterial = ageOfMaterial
        self.museumId = museumId
        self.specimenReference = specimenReference
        # one-to-many
        # SpecimenLink: Sample.idSpecimen REFERENCES Specimen(idSpecimen)
        self.sampleSpecimenList = []
        # fk_idSpecimen: SpecimenDataRelease.idSpecimen REFERENCES Specimen(idSpecimen)
        self.specimenDataReleaseSpecimenList = []
        # SpecimenGeneStats_ibfk_1: SpecimenGeneStats.idSpecimen REFERENCES Specimen(idSpecimen)
        self.specimenGeneStatsSpecimenList = []
        # SpecimenRawReads_ibfk_1: SpecimenRawReads.idSpecimen REFERENCES Specimen(idSpecimen)
        self.specimenRawReadsSpecimenList = []
        # no python attribute: SourceSpecimen

    def insertIntoDatabase(self, cursor):
        sqlCmd = 'INSERT INTO `Specimen` (`idGenus`, `idPaftol`, `idSpecies`, `BankID`, `LCD`, `MSB`, `Collector`, `CollectorNo`, `VoucherNo`, `MuseumBarcode`, `OldSpeciesName`, `idSourceSpecimen`, `idProject`, `idOriginCountry`, `idMaterialSource`, `AgeOfMaterial`, `idMuseumID`, `SpecimenReference`) VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s)'
        l = []
        l.append(None if self.genus is None else self.genus.idGenus)
        l.append(self.idPaftol)
        l.append(None if self.species is None else self.species.idSpecies)
        l.append(self.bankId)
        l.append(self.lcd)
        l.append(self.msb)
        l.append(self.collector)
        l.append(self.collectorNo)
        l.append(self.voucherNo)
        l.append(self.museumBarcode)
        l.append(self.oldSpeciesName)
        l.append(None if self.sourceSpecimen is None else self.sourceSpecimen.idSourceSpecimen)
        l.append(None if self.project is None else self.project.idProject)
        l.append(self.idOriginCountry)
        l.append(None if self.materialSource is None else self.materialSource.idMaterialSource)
        l.append(self.ageOfMaterial)
        l.append(None if self.museumId is None else self.museumId.idMuseumId)
        l.append(self.specimenReference)
        cursor.execute(sqlCmd, tuple(l))


class SpecimenDataRelease(object):

    def __init__(self, specimen=None, dataRelease=None):
        self.idSpecimenDataRelease = None
        self.specimen = specimen
        self.dataRelease = dataRelease
        # one-to-many

    def insertIntoDatabase(self, cursor):
        sqlCmd = 'INSERT INTO `SpecimenDataRelease` (`idSpecimen`, `idDataRelease`) VALUES (%s, %s)'
        l = []
        l.append(None if self.specimen is None else self.specimen.idSpecimen)
        l.append(None if self.dataRelease is None else self.dataRelease.idDataRelease)
        cursor.execute(sqlCmd, tuple(l))


class SpecimenGeneStats(object):

    def __init__(self, numGene=None, recoveredLength=None):
        self.specimen = None
        self.numGene = numGene
        self.recoveredLength = recoveredLength
        # one-to-many

    def insertIntoDatabase(self, cursor):
        sqlCmd = 'INSERT INTO `SpecimenGeneStats` (`NumGene`, `RecoveredLength`) VALUES (%s, %s)'
        l = []
        l.append(None if self.specimen is None else self.specimen.idSpecimen)
        l.append(self.numGene)
        l.append(self.recoveredLength)
        cursor.execute(sqlCmd, tuple(l))


class SpecimenRawReads(object):

    def __init__(self, specimen=None, numReads=None, seqPlatform=None, enaExpNum=None, enaRunNum=None):
        self.id = None
        self.specimen = specimen
        self.numReads = numReads
        self.seqPlatform = seqPlatform
        self.enaExpNum = enaExpNum
        self.enaRunNum = enaRunNum
        # one-to-many

    def insertIntoDatabase(self, cursor):
        sqlCmd = 'INSERT INTO `SpecimenRawReads` (`idSpecimen`, `NumReads`, `SeqPlatform`, `ENAExpNum`, `ENARunNum`) VALUES (%s, %s, %s, %s, %s)'
        l = []
        l.append(None if self.specimen is None else self.specimen.idSpecimen)
        l.append(self.numReads)
        l.append(self.seqPlatform)
        l.append(self.enaExpNum)
        l.append(self.enaRunNum)
        cursor.execute(sqlCmd, tuple(l))


class Status(object):

    def __init__(self, status=None):
        self.idStatus = None
        self.status = status
        # one-to-many
        # StatusLookup: Library.idStatus REFERENCES Status(idStatus)
        self.libraryStatusList = []
        # fk_status: Sequence.idStatus REFERENCES Status(idStatus)
        self.sequenceStatusList = []

    def insertIntoDatabase(self, cursor):
        sqlCmd = 'INSERT INTO `Status` (`Status`) VALUES (%s)'
        l = []
        l.append(self.status)
        cursor.execute(sqlCmd, tuple(l))


class Iso_country_3166_1(object):

    def __init__(self, country=None, alpha_2_code=None, alpha_3_code=None, country_code=None, iso_3166_2_code=None, region=None, sub_region=None, intermediate_region=None, region_code=None, sub_region_code=None, intermediate_region_code=None):
        self.iso_country_id = None
        self.country = country
        self.alpha_2_code = alpha_2_code
        self.alpha_3_code = alpha_3_code
        self.country_code = country_code
        self.iso_3166_2_code = iso_3166_2_code
        self.region = region
        self.sub_region = sub_region
        self.intermediate_region = intermediate_region
        self.region_code = region_code
        self.sub_region_code = sub_region_code
        self.intermediate_region_code = intermediate_region_code
        # one-to-many

    def insertIntoDatabase(self, cursor):
        sqlCmd = 'INSERT INTO `iso_country_3166_1` (`country`, `alpha_2_code`, `alpha_3_code`, `country_code`, `iso_3166_2_code`, `region`, `sub_region`, `intermediate_region`, `region_code`, `sub_region_code`, `intermediate_region_code`) VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s)'
        l = []
        l.append(self.country)
        l.append(self.alpha_2_code)
        l.append(self.alpha_3_code)
        l.append(self.country_code)
        l.append(self.iso_3166_2_code)
        l.append(self.region)
        l.append(self.sub_region)
        l.append(self.intermediate_region)
        l.append(self.region_code)
        l.append(self.sub_region_code)
        l.append(self.intermediate_region_code)
        cursor.execute(sqlCmd, tuple(l))


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


def loadBlacklistedReasonDict(connection, productionDatabase):
    cursor = connection.cursor()
    entityDict = {}
    sqlStatement = 'SELECT `idBlacklistedReason`, `BlacklistedReason` FROM `BlacklistedReason`'
    cursor.execute(sqlStatement)
    for row in cursor:
        entity = BlacklistedReason()
        entity.idBlacklistedReason = paftol.database.intOrNone(row[0])
        entity.blacklistedReason = paftol.database.strOrNone(row[1])
        entityDict[entity.idBlacklistedReason] = entity
    cursor.close()
    return entityDict


def loadCoordinatesDict(connection, productionDatabase):
    cursor = connection.cursor()
    entityDict = {}
    sqlStatement = 'SELECT `idCoordinate`, `Coordinate` FROM `Coordinates`'
    cursor.execute(sqlStatement)
    for row in cursor:
        entity = Coordinates()
        entity.idCoordinate = paftol.database.intOrNone(row[0])
        entity.coordinate = paftol.database.strOrNone(row[1])
        entityDict[entity.idCoordinate] = entity
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


def loadDataReleaseDict(connection, productionDatabase):
    cursor = connection.cursor()
    entityDict = {}
    sqlStatement = 'SELECT `idDataRelease`, `ReleaseNumber`, `DataRelease`, `TaxonCount` FROM `DataRelease`'
    cursor.execute(sqlStatement)
    for row in cursor:
        entity = DataRelease()
        entity.idDataRelease = paftol.database.intOrNone(row[0])
        entity.releaseNumber = paftol.database.floatOrNone(row[1])
        entity.dataRelease = paftol.database.strOrNone(row[2])
        entity.taxonCount = paftol.database.strOrNone(row[3])
        entityDict[entity.idDataRelease] = entity
    cursor.close()
    return entityDict


def loadDataSourceDict(connection, productionDatabase):
    cursor = connection.cursor()
    entityDict = {}
    sqlStatement = 'SELECT `idDataSource`, `DataSource` FROM `DataSource`'
    cursor.execute(sqlStatement)
    for row in cursor:
        entity = DataSource()
        entity.idDataSource = paftol.database.intOrNone(row[0])
        entity.dataSource = paftol.database.strOrNone(row[1])
        entityDict[entity.idDataSource] = entity
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
            raise Exception('no Order entity with idOrder = %d' % entityId)
        else:
            entity.order = productionDatabase.orderDict[entityId]
            # type: int, name: idOrder, foreignTable: Order, foreignColumn: idOrder
            entity.order.familyOrderList.append(entity)
        entityDict[entity.idFamily] = entity
    cursor.close()
    return entityDict


def loadGeneStatsDict(connection, productionDatabase):
    cursor = connection.cursor()
    entityDict = {}
    sqlStatement = 'SELECT `id`, `InternalName`, `ExemplarAccession`, `ExemplarName`, `ExemplarSpecies`, `ExemplarHyperlink`, `NewickFile`, `NewickFilePathName`, `AverageContigLength`, `Depth`, `AverageContigLengthPercentage`, `NumSeq`, `NumGenera`, `NumSpecies` FROM `GeneStats`'
    cursor.execute(sqlStatement)
    for row in cursor:
        entity = GeneStats()
        entity.id = paftol.database.intOrNone(row[0])
        entity.internalName = paftol.database.strOrNone(row[1])
        entity.exemplarAccession = paftol.database.strOrNone(row[2])
        entity.exemplarName = paftol.database.strOrNone(row[3])
        entity.exemplarSpecies = paftol.database.strOrNone(row[4])
        entity.exemplarHyperlink = paftol.database.strOrNone(row[5])
        entity.newickFile = paftol.database.strOrNone(row[6])
        entity.newickFilePathName = paftol.database.strOrNone(row[7])
        entity.averageContigLength = paftol.database.floatOrNone(row[8])
        entity.depth = paftol.database.intOrNone(row[9])
        entity.averageContigLengthPercentage = paftol.database.floatOrNone(row[10])
        entity.numSeq = paftol.database.intOrNone(row[11])
        entity.numGenera = paftol.database.intOrNone(row[12])
        entity.numSpecies = paftol.database.intOrNone(row[13])
        entityDict[entity.id] = entity
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
            raise Exception('no Family entity with idFamily = %d' % entityId)
        else:
            entity.family = productionDatabase.familyDict[entityId]
            # type: int, name: idFamily, foreignTable: Family, foreignColumn: idFamily
            entity.family.genusFamilyList.append(entity)
        entity.genus = paftol.database.strOrNone(row[2])
        # many to one: source
        entityId = paftol.database.intOrNone(row[3])
        if entityId is None:
            entity.source = None
        elif entityId not in productionDatabase.sourceDict:
            raise Exception('no Source entity with idSource = %d' % entityId)
        else:
            entity.source = productionDatabase.sourceDict[entityId]
            # type: int, name: idSource, foreignTable: Source, foreignColumn: idSource
            entity.source.genusSourceList.append(entity)
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
    sqlStatement = 'SELECT `idIndexes`, `Indexes`, `IndexNameFwd`, `SeqFwdPlatform1`, `SeqFwdPlatform2`, `IndexNameRv`, `SeqRv` FROM `Indexes`'
    cursor.execute(sqlStatement)
    for row in cursor:
        entity = Indexes()
        entity.idIndexes = paftol.database.intOrNone(row[0])
        entity.indexes = paftol.database.strOrNone(row[1])
        entity.indexNameFwd = paftol.database.strOrNone(row[2])
        entity.seqFwdPlatform1 = paftol.database.strOrNone(row[3])
        entity.seqFwdPlatform2 = paftol.database.strOrNone(row[4])
        entity.indexNameRv = paftol.database.strOrNone(row[5])
        entity.seqRv = paftol.database.strOrNone(row[6])
        entityDict[entity.idIndexes] = entity
    cursor.close()
    return entityDict


def loadLibraryDict(connection, productionDatabase):
    cursor = connection.cursor()
    entityDict = {}
    sqlStatement = 'SELECT `idLibrary`, `idSample`, `LibConcentration`, `LibQuality`, `RemainingVolume`, `LibTapeStation`, `Sonication`, `Plate No`, `idCoordinate`, `Description`, `idStatus`, `idIndexes`, `GenerateLibrary` FROM `Library`'
    cursor.execute(sqlStatement)
    for row in cursor:
        entity = Library()
        entity.idLibrary = paftol.database.intOrNone(row[0])
        # many to one: sample
        entityId = paftol.database.intOrNone(row[1])
        if entityId is None:
            entity.sample = None
        elif entityId not in productionDatabase.sampleDict:
            raise Exception('no Sample entity with idSample = %d' % entityId)
        else:
            entity.sample = productionDatabase.sampleDict[entityId]
            # type: int, name: idSample, foreignTable: Sample, foreignColumn: idSample
            entity.sample.librarySampleList.append(entity)
        entity.libConcentration = paftol.database.floatOrNone(row[2])
        entity.libQuality = paftol.database.strOrNone(row[3])
        entity.remainingVolume = paftol.database.floatOrNone(row[4])
        entity.libTapeStation = paftol.database.strOrNone(row[5])
        entity.sonication = paftol.database.intOrNone(row[6])
        entity.plate_No = paftol.database.strOrNone(row[7])
        # many to one: coordinate
        entityId = paftol.database.intOrNone(row[8])
        if entityId is None:
            entity.coordinate = None
        elif entityId not in productionDatabase.coordinatesDict:
            raise Exception('no Coordinates entity with idCoordinate = %d' % entityId)
        else:
            entity.coordinate = productionDatabase.coordinatesDict[entityId]
            # type: int, name: idCoordinate, foreignTable: Coordinates, foreignColumn: idCoordinate
            entity.coordinate.libraryCoordinateList.append(entity)
        entity.description = paftol.database.strOrNone(row[9])
        # many to one: status
        entityId = paftol.database.intOrNone(row[10])
        if entityId is None:
            entity.status = None
        elif entityId not in productionDatabase.statusDict:
            raise Exception('no Status entity with idStatus = %d' % entityId)
        else:
            entity.status = productionDatabase.statusDict[entityId]
            # type: int, name: idStatus, foreignTable: Status, foreignColumn: idStatus
            entity.status.libraryStatusList.append(entity)
        # many to one: indexes
        entityId = paftol.database.intOrNone(row[11])
        if entityId is None:
            entity.indexes = None
        elif entityId not in productionDatabase.indexesDict:
            raise Exception('no Indexes entity with idIndexes = %d' % entityId)
        else:
            entity.indexes = productionDatabase.indexesDict[entityId]
            # type: int, name: idIndexes, foreignTable: Indexes, foreignColumn: idIndexes
            entity.indexes.libraryIndexesList.append(entity)
        entity.generateLibrary = paftol.database.intOrNone(row[12])
        entityDict[entity.idLibrary] = entity
    cursor.close()
    return entityDict


def loadLocationDict(connection, productionDatabase):
    cursor = connection.cursor()
    entityDict = {}
    sqlStatement = 'SELECT `idLocation`, `Location` FROM `Location`'
    cursor.execute(sqlStatement)
    for row in cursor:
        entity = Location()
        entity.idLocation = paftol.database.intOrNone(row[0])
        entity.location = paftol.database.strOrNone(row[1])
        entityDict[entity.idLocation] = entity
    cursor.close()
    return entityDict


def loadLogDict(connection, productionDatabase):
    cursor = connection.cursor()
    entityDict = {}
    sqlStatement = 'SELECT `idLog`, `TransactionDate`, `User`, `File`, `TransactionType` FROM `Log`'
    cursor.execute(sqlStatement)
    for row in cursor:
        entity = Log()
        entity.idLog = paftol.database.intOrNone(row[0])
        entity.transactionDate = paftol.database.strOrNone(row[1])
        entity.user = paftol.database.strOrNone(row[2])
        entity.file = paftol.database.strOrNone(row[3])
        entity.transactionType = paftol.database.strOrNone(row[4])
        entityDict[entity.idLog] = entity
    cursor.close()
    return entityDict


def loadMaterialSourceDict(connection, productionDatabase):
    cursor = connection.cursor()
    entityDict = {}
    sqlStatement = 'SELECT `idMaterialSource`, `MaterialSource` FROM `MaterialSource`'
    cursor.execute(sqlStatement)
    for row in cursor:
        entity = MaterialSource()
        entity.idMaterialSource = paftol.database.intOrNone(row[0])
        entity.materialSource = paftol.database.strOrNone(row[1])
        entityDict[entity.idMaterialSource] = entity
    cursor.close()
    return entityDict


def loadMigrationsLogDict(connection, productionDatabase):
    cursor = connection.cursor()
    entityDict = {}
    sqlStatement = 'SELECT `idMigrationsLog`, `MigrationName`, `MigrationDate` FROM `MigrationsLog`'
    cursor.execute(sqlStatement)
    for row in cursor:
        entity = MigrationsLog()
        entity.idMigrationsLog = paftol.database.intOrNone(row[0])
        entity.migrationName = paftol.database.strOrNone(row[1])
        entity.migrationDate = paftol.database.strOrNone(row[2])
        entityDict[entity.idMigrationsLog] = entity
    cursor.close()
    return entityDict


def loadMuseumDict(connection, productionDatabase):
    cursor = connection.cursor()
    entityDict = {}
    sqlStatement = 'SELECT `idMuseumID`, `MuseumID`, `MuseumName` FROM `Museum`'
    cursor.execute(sqlStatement)
    for row in cursor:
        entity = Museum()
        entity.idMuseumId = paftol.database.intOrNone(row[0])
        entity.museumId = paftol.database.strOrNone(row[1])
        entity.museumName = paftol.database.strOrNone(row[2])
        entityDict[entity.idMuseumId] = entity
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


def loadPlatformDict(connection, productionDatabase):
    cursor = connection.cursor()
    entityDict = {}
    sqlStatement = 'SELECT `idPlatform`, `Platform` FROM `Platform`'
    cursor.execute(sqlStatement)
    for row in cursor:
        entity = Platform()
        entity.idPlatform = paftol.database.intOrNone(row[0])
        entity.platform = paftol.database.strOrNone(row[1])
        entityDict[entity.idPlatform] = entity
    cursor.close()
    return entityDict


def loadProjectDict(connection, productionDatabase):
    cursor = connection.cursor()
    entityDict = {}
    sqlStatement = 'SELECT `idProject`, `Project`, `idDataSource` FROM `Project`'
    cursor.execute(sqlStatement)
    for row in cursor:
        entity = Project()
        entity.idProject = paftol.database.intOrNone(row[0])
        entity.project = paftol.database.strOrNone(row[1])
        # many to one: dataSource
        entityId = paftol.database.intOrNone(row[2])
        if entityId is None:
            entity.dataSource = None
        elif entityId not in productionDatabase.dataSourceDict:
            raise Exception('no DataSource entity with idDataSource = %d' % entityId)
        else:
            entity.dataSource = productionDatabase.dataSourceDict[entityId]
            # type: int, name: idDataSource, foreignTable: DataSource, foreignColumn: idDataSource
            entity.dataSource.projectDataSourceList.append(entity)
        entityDict[entity.idProject] = entity
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
    sqlStatement = 'SELECT `idSample`, `idSpecimen`, `Description`, `idAction`, `idExtractionType`, `idQuality`, `SampleConcentration`, `NewSampleConcentration`, `GelImage`, `SampleTapeStation`, `idDNAVolume`, `Compliant`, `ENASampleNum`, `SecENASampleNum` FROM `Sample`'
    cursor.execute(sqlStatement)
    for row in cursor:
        entity = Sample()
        entity.idSample = paftol.database.intOrNone(row[0])
        # many to one: specimen
        entityId = paftol.database.intOrNone(row[1])
        if entityId is None:
            entity.specimen = None
        elif entityId not in productionDatabase.specimenDict:
            raise Exception('no Specimen entity with idSpecimen = %d' % entityId)
        else:
            entity.specimen = productionDatabase.specimenDict[entityId]
            # type: int, name: idSpecimen, foreignTable: Specimen, foreignColumn: idSpecimen
            entity.specimen.sampleSpecimenList.append(entity)
        entity.description = paftol.database.strOrNone(row[2])
        # many to one: action
        entityId = paftol.database.intOrNone(row[3])
        if entityId is None:
            entity.action = None
        elif entityId not in productionDatabase.actionDict:
            raise Exception('no Action entity with idAction = %d' % entityId)
        else:
            entity.action = productionDatabase.actionDict[entityId]
            # type: int, name: idAction, foreignTable: Action, foreignColumn: idAction
            entity.action.sampleActionList.append(entity)
        # many to one: extractionType
        entityId = paftol.database.intOrNone(row[4])
        if entityId is None:
            entity.extractionType = None
        elif entityId not in productionDatabase.extractionTypeDict:
            raise Exception('no ExtractionType entity with idExtractionType = %d' % entityId)
        else:
            entity.extractionType = productionDatabase.extractionTypeDict[entityId]
            # type: int, name: idExtractionType, foreignTable: ExtractionType, foreignColumn: idExtractionType
            entity.extractionType.sampleExtractionTypeList.append(entity)
        # many to one: quality
        entityId = paftol.database.intOrNone(row[5])
        if entityId is None:
            entity.quality = None
        elif entityId not in productionDatabase.qualityDict:
            raise Exception('no Quality entity with idQuality = %d' % entityId)
        else:
            entity.quality = productionDatabase.qualityDict[entityId]
            # type: int, name: idQuality, foreignTable: Quality, foreignColumn: idQuality
            entity.quality.sampleQualityList.append(entity)
        entity.sampleConcentration = paftol.database.floatOrNone(row[6])
        entity.newSampleConcentration = paftol.database.floatOrNone(row[7])
        entity.gelImage = paftol.database.strOrNone(row[8])
        entity.sampleTapeStation = paftol.database.strOrNone(row[9])
        # many to one: dnaVolume
        entityId = paftol.database.intOrNone(row[10])
        if entityId is None:
            entity.dnaVolume = None
        elif entityId not in productionDatabase.dNAVolumeDict:
            raise Exception('no DNAVolume entity with idDnaVolume = %d' % entityId)
        else:
            entity.dnaVolume = productionDatabase.dNAVolumeDict[entityId]
            # type: int, name: idDNAVolume, foreignTable: DNAVolume, foreignColumn: idDNAVolume
            entity.dnaVolume.sampleDnaVolumeList.append(entity)
        entity.compliant = paftol.database.strOrNone(row[11])
        entity.enaSampleNum = paftol.database.strOrNone(row[12])
        entity.secEnaSampleNum = paftol.database.strOrNone(row[13])
        entityDict[entity.idSample] = entity
    cursor.close()
    return entityDict


def loadSequenceDict(connection, productionDatabase):
    cursor = connection.cursor()
    entityDict = {}
    sqlStatement = 'SELECT `idSequencing`, `SequenceID`, `idLibrary`, `idPlatform`, `idLocation`, `SequencingRun`, `NumInferredCds`, `MedianHybpiperCdsLength`, `idStatus`, `HybridisationPool`, `R2FastqFile`, `R1FastqFile`, `Blacklisted`, `idBlacklistedReason`, `idSequencingStrategy`, `ENAExpNumber`, `ENARunNumber`, `ExternalSequenceID` FROM `Sequence`'
    cursor.execute(sqlStatement)
    for row in cursor:
        entity = Sequence()
        entity.idSequencing = paftol.database.intOrNone(row[0])
        entity.sequenceId = paftol.database.intOrNone(row[1])
        # many to one: library
        entityId = paftol.database.intOrNone(row[2])
        if entityId is None:
            entity.library = None
        elif entityId not in productionDatabase.libraryDict:
            raise Exception('no Library entity with idLibrary = %d' % entityId)
        else:
            entity.library = productionDatabase.libraryDict[entityId]
            # type: int, name: idLibrary, foreignTable: Library, foreignColumn: idLibrary
            entity.library.sequenceLibraryList.append(entity)
        # many to one: platform
        entityId = paftol.database.intOrNone(row[3])
        if entityId is None:
            entity.platform = None
        elif entityId not in productionDatabase.platformDict:
            raise Exception('no Platform entity with idPlatform = %d' % entityId)
        else:
            entity.platform = productionDatabase.platformDict[entityId]
            # type: int, name: idPlatform, foreignTable: Platform, foreignColumn: idPlatform
            entity.platform.sequencePlatformList.append(entity)
        # many to one: location
        entityId = paftol.database.intOrNone(row[4])
        if entityId is None:
            entity.location = None
        elif entityId not in productionDatabase.locationDict:
            raise Exception('no Location entity with idLocation = %d' % entityId)
        else:
            entity.location = productionDatabase.locationDict[entityId]
            # type: int, name: idLocation, foreignTable: Location, foreignColumn: idLocation
            entity.location.sequenceLocationList.append(entity)
        entity.sequencingRun = paftol.database.strOrNone(row[5])
        entity.numInferredCds = paftol.database.intOrNone(row[6])
        entity.medianHybpiperCdsLength = paftol.database.floatOrNone(row[7])
        # many to one: status
        entityId = paftol.database.intOrNone(row[8])
        if entityId is None:
            entity.status = None
        elif entityId not in productionDatabase.statusDict:
            raise Exception('no Status entity with idStatus = %d' % entityId)
        else:
            entity.status = productionDatabase.statusDict[entityId]
            # type: int, name: idStatus, foreignTable: Status, foreignColumn: idStatus
            entity.status.sequenceStatusList.append(entity)
        entity.hybridisationPool = paftol.database.strOrNone(row[9])
        entity.r2FastqFile = paftol.database.strOrNone(row[10])
        entity.r1FastqFile = paftol.database.strOrNone(row[11])
        entity.blacklisted = paftol.database.intOrNone(row[12])
        # many to one: blacklistedReason
        entityId = paftol.database.intOrNone(row[13])
        if entityId is None:
            entity.blacklistedReason = None
        elif entityId not in productionDatabase.blacklistedReasonDict:
            raise Exception('no BlacklistedReason entity with idBlacklistedReason = %d' % entityId)
        else:
            entity.blacklistedReason = productionDatabase.blacklistedReasonDict[entityId]
            # type: int, name: idBlacklistedReason, foreignTable: BlacklistedReason, foreignColumn: idBlacklistedReason
            entity.blacklistedReason.sequenceBlacklistedReasonList.append(entity)
        # many to one: sequencingStrategy
        entityId = paftol.database.intOrNone(row[14])
        if entityId is None:
            entity.sequencingStrategy = None
        elif entityId not in productionDatabase.sequencingStrategyDict:
            raise Exception('no SequencingStrategy entity with idSequencingStrategy = %d' % entityId)
        else:
            entity.sequencingStrategy = productionDatabase.sequencingStrategyDict[entityId]
            # type: int, name: idSequencingStrategy, foreignTable: SequencingStrategy, foreignColumn: idSequencingStrategy
            entity.sequencingStrategy.sequenceSequencingStrategyList.append(entity)
        entity.enaExpNumber = paftol.database.strOrNone(row[15])
        entity.enaRunNumber = paftol.database.strOrNone(row[16])
        entity.externalSequenceId = paftol.database.strOrNone(row[17])
        entityDict[entity.idSequencing] = entity
    cursor.close()
    return entityDict


def loadSequenceDataReleaseDict(connection, productionDatabase):
    cursor = connection.cursor()
    entityDict = {}
    sqlStatement = 'SELECT `idSequenceDataRelease`, `idSequencing`, `idDataRelease` FROM `SequenceDataRelease`'
    cursor.execute(sqlStatement)
    for row in cursor:
        entity = SequenceDataRelease()
        entity.idSequenceDataRelease = paftol.database.intOrNone(row[0])
        # many to one: sequencing
        entityId = paftol.database.intOrNone(row[1])
        if entityId is None:
            entity.sequencing = None
        elif entityId not in productionDatabase.sequenceDict:
            raise Exception('no Sequence entity with idSequencing = %d' % entityId)
        else:
            entity.sequencing = productionDatabase.sequenceDict[entityId]
            # type: int, name: idSequencing, foreignTable: Sequence, foreignColumn: idSequencing
            entity.sequencing.sequenceDataReleaseSequencingList.append(entity)
        # many to one: dataRelease
        entityId = paftol.database.intOrNone(row[2])
        if entityId is None:
            entity.dataRelease = None
        elif entityId not in productionDatabase.dataReleaseDict:
            raise Exception('no DataRelease entity with idDataRelease = %d' % entityId)
        else:
            entity.dataRelease = productionDatabase.dataReleaseDict[entityId]
            # type: int, name: idDataRelease, foreignTable: DataRelease, foreignColumn: idDataRelease
            entity.dataRelease.sequenceDataReleaseDataReleaseList.append(entity)
        entityDict[entity.idSequenceDataRelease] = entity
    cursor.close()
    return entityDict


def loadSequencingStrategyDict(connection, productionDatabase):
    cursor = connection.cursor()
    entityDict = {}
    sqlStatement = 'SELECT `idSequencingStrategy`, `SequencingStrategy` FROM `SequencingStrategy`'
    cursor.execute(sqlStatement)
    for row in cursor:
        entity = SequencingStrategy()
        entity.idSequencingStrategy = paftol.database.intOrNone(row[0])
        entity.sequencingStrategy = paftol.database.strOrNone(row[1])
        entityDict[entity.idSequencingStrategy] = entity
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
    sqlStatement = 'SELECT `idSpecimen`, `idGenus`, `idPaftol`, `idSpecies`, `BankID`, `LCD`, `MSB`, `Collector`, `CollectorNo`, `VoucherNo`, `MuseumBarcode`, `OldSpeciesName`, `idSourceSpecimen`, `idProject`, `idOriginCountry`, `idMaterialSource`, `AgeOfMaterial`, `idMuseumID`, `SpecimenReference` FROM `Specimen`'
    cursor.execute(sqlStatement)
    for row in cursor:
        entity = Specimen()
        entity.idSpecimen = paftol.database.intOrNone(row[0])
        # many to one: genus
        entityId = paftol.database.intOrNone(row[1])
        if entityId is None:
            entity.genus = None
        elif entityId not in productionDatabase.genusDict:
            raise Exception('no Genus entity with idGenus = %d' % entityId)
        else:
            entity.genus = productionDatabase.genusDict[entityId]
            # type: int, name: idGenus, foreignTable: Genus, foreignColumn: idGenus
            entity.genus.specimenGenusList.append(entity)
        entity.idPaftol = paftol.database.intOrNone(row[2])
        # many to one: species
        entityId = paftol.database.intOrNone(row[3])
        if entityId is None:
            entity.species = None
        elif entityId not in productionDatabase.speciesDict:
            raise Exception('no Species entity with idSpecies = %d' % entityId)
        else:
            entity.species = productionDatabase.speciesDict[entityId]
            # type: int, name: idSpecies, foreignTable: Species, foreignColumn: idSpecies
            entity.species.specimenSpeciesList.append(entity)
        entity.bankId = paftol.database.intOrNone(row[4])
        entity.lcd = paftol.database.strOrNone(row[5])
        entity.msb = paftol.database.intOrNone(row[6])
        entity.collector = paftol.database.strOrNone(row[7])
        entity.collectorNo = paftol.database.strOrNone(row[8])
        entity.voucherNo = paftol.database.strOrNone(row[9])
        entity.museumBarcode = paftol.database.strOrNone(row[10])
        entity.oldSpeciesName = paftol.database.strOrNone(row[11])
        # many to one: sourceSpecimen
        entityId = paftol.database.intOrNone(row[12])
        if entityId is None:
            entity.sourceSpecimen = None
        elif entityId not in productionDatabase.sourceSpecimenDict:
            raise Exception('no SourceSpecimen entity with idSourceSpecimen = %d' % entityId)
        else:
            entity.sourceSpecimen = productionDatabase.sourceSpecimenDict[entityId]
            # type: int, name: idSourceSpecimen, foreignTable: SourceSpecimen, foreignColumn: idSourceSpecimen
            entity.sourceSpecimen.specimenSourceSpecimenList.append(entity)
        # many to one: project
        entityId = paftol.database.intOrNone(row[13])
        if entityId is None:
            entity.project = None
        elif entityId not in productionDatabase.projectDict:
            raise Exception('no Project entity with idProject = %d' % entityId)
        else:
            entity.project = productionDatabase.projectDict[entityId]
            # type: int, name: idProject, foreignTable: Project, foreignColumn: idProject
            entity.project.specimenProjectList.append(entity)
        entity.idOriginCountry = paftol.database.intOrNone(row[14])
        # many to one: materialSource
        entityId = paftol.database.intOrNone(row[15])
        if entityId is None:
            entity.materialSource = None
        elif entityId not in productionDatabase.materialSourceDict:
            raise Exception('no MaterialSource entity with idMaterialSource = %d' % entityId)
        else:
            entity.materialSource = productionDatabase.materialSourceDict[entityId]
            # type: int, name: idMaterialSource, foreignTable: MaterialSource, foreignColumn: idMaterialSource
            entity.materialSource.specimenMaterialSourceList.append(entity)
        entity.ageOfMaterial = paftol.database.intOrNone(row[16])
        # many to one: museumId
        entityId = paftol.database.intOrNone(row[17])
        if entityId is None:
            entity.museumId = None
        elif entityId not in productionDatabase.museumDict:
            raise Exception('no Museum entity with idMuseumId = %d' % entityId)
        else:
            entity.museumId = productionDatabase.museumDict[entityId]
            # type: int, name: idMuseumID, foreignTable: Museum, foreignColumn: idMuseumID
            entity.museumId.specimenMuseumIdList.append(entity)
        entity.specimenReference = paftol.database.strOrNone(row[18])
        entityDict[entity.idSpecimen] = entity
    cursor.close()
    return entityDict


def loadSpecimenDataReleaseDict(connection, productionDatabase):
    cursor = connection.cursor()
    entityDict = {}
    sqlStatement = 'SELECT `idSpecimenDataRelease`, `idSpecimen`, `idDataRelease` FROM `SpecimenDataRelease`'
    cursor.execute(sqlStatement)
    for row in cursor:
        entity = SpecimenDataRelease()
        entity.idSpecimenDataRelease = paftol.database.intOrNone(row[0])
        # many to one: specimen
        entityId = paftol.database.intOrNone(row[1])
        if entityId is None:
            entity.specimen = None
        elif entityId not in productionDatabase.specimenDict:
            raise Exception('no Specimen entity with idSpecimen = %d' % entityId)
        else:
            entity.specimen = productionDatabase.specimenDict[entityId]
            # type: int, name: idSpecimen, foreignTable: Specimen, foreignColumn: idSpecimen
            entity.specimen.specimenDataReleaseSpecimenList.append(entity)
        # many to one: dataRelease
        entityId = paftol.database.intOrNone(row[2])
        if entityId is None:
            entity.dataRelease = None
        elif entityId not in productionDatabase.dataReleaseDict:
            raise Exception('no DataRelease entity with idDataRelease = %d' % entityId)
        else:
            entity.dataRelease = productionDatabase.dataReleaseDict[entityId]
            # type: int, name: idDataRelease, foreignTable: DataRelease, foreignColumn: idDataRelease
            entity.dataRelease.specimenDataReleaseDataReleaseList.append(entity)
        entityDict[entity.idSpecimenDataRelease] = entity
    cursor.close()
    return entityDict


def loadSpecimenGeneStatsDict(connection, productionDatabase):
    cursor = connection.cursor()
    entityDict = {}
    sqlStatement = 'SELECT `idSpecimen`, `NumGene`, `RecoveredLength` FROM `SpecimenGeneStats`'
    cursor.execute(sqlStatement)
    for row in cursor:
        entity = SpecimenGeneStats()
        # many to one: specimen
        entityId = paftol.database.intOrNone(row[0])
        if entityId is None:
            entity.specimen = None
        elif entityId not in productionDatabase.specimenDict:
            raise Exception('no Specimen entity with idSpecimen = %d' % entityId)
        else:
            entity.specimen = productionDatabase.specimenDict[entityId]
            # type: int, name: idSpecimen, foreignTable: Specimen, foreignColumn: idSpecimen
            entity.specimen.specimenGeneStatsSpecimenList.append(entity)
        entity.numGene = paftol.database.intOrNone(row[1])
        entity.recoveredLength = paftol.database.floatOrNone(row[2])
        entityDict[entity.specimen] = entity
    cursor.close()
    return entityDict


def loadSpecimenRawReadsDict(connection, productionDatabase):
    cursor = connection.cursor()
    entityDict = {}
    sqlStatement = 'SELECT `id`, `idSpecimen`, `NumReads`, `SeqPlatform`, `ENAExpNum`, `ENARunNum` FROM `SpecimenRawReads`'
    cursor.execute(sqlStatement)
    for row in cursor:
        entity = SpecimenRawReads()
        entity.id = paftol.database.intOrNone(row[0])
        # many to one: specimen
        entityId = paftol.database.intOrNone(row[1])
        if entityId is None:
            entity.specimen = None
        elif entityId not in productionDatabase.specimenDict:
            raise Exception('no Specimen entity with idSpecimen = %d' % entityId)
        else:
            entity.specimen = productionDatabase.specimenDict[entityId]
            # type: int, name: idSpecimen, foreignTable: Specimen, foreignColumn: idSpecimen
            entity.specimen.specimenRawReadsSpecimenList.append(entity)
        entity.numReads = paftol.database.intOrNone(row[2])
        entity.seqPlatform = paftol.database.strOrNone(row[3])
        entity.enaExpNum = paftol.database.strOrNone(row[4])
        entity.enaRunNum = paftol.database.strOrNone(row[5])
        entityDict[entity.id] = entity
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


def loadIso_country_3166_1Dict(connection, productionDatabase):
    cursor = connection.cursor()
    entityDict = {}
    sqlStatement = 'SELECT `iso_country_id`, `country`, `alpha_2_code`, `alpha_3_code`, `country_code`, `iso_3166_2_code`, `region`, `sub_region`, `intermediate_region`, `region_code`, `sub_region_code`, `intermediate_region_code` FROM `iso_country_3166_1`'
    cursor.execute(sqlStatement)
    for row in cursor:
        entity = Iso_country_3166_1()
        entity.iso_country_id = paftol.database.intOrNone(row[0])
        entity.country = paftol.database.strOrNone(row[1])
        entity.alpha_2_code = paftol.database.strOrNone(row[2])
        entity.alpha_3_code = paftol.database.strOrNone(row[3])
        entity.country_code = paftol.database.intOrNone(row[4])
        entity.iso_3166_2_code = paftol.database.strOrNone(row[5])
        entity.region = paftol.database.strOrNone(row[6])
        entity.sub_region = paftol.database.strOrNone(row[7])
        entity.intermediate_region = paftol.database.strOrNone(row[8])
        entity.region_code = paftol.database.intOrNone(row[9])
        entity.sub_region_code = paftol.database.intOrNone(row[10])
        entity.intermediate_region_code = paftol.database.intOrNone(row[11])
        entityDict[entity.iso_country_id] = entity
    cursor.close()
    return entityDict


class ProductionDatabase(object):

    def __init__(self, connection):
        self.actionDict = {}
        self.blacklistedReasonDict = {}
        self.coordinatesDict = {}
        self.dBVersionDict = {}
        self.dNAVolumeDict = {}
        self.dataReleaseDict = {}
        self.dataSourceDict = {}
        self.extractionTypeDict = {}
        self.familyDict = {}
        self.geneStatsDict = {}
        self.genusDict = {}
        self.indexesDict = {}
        self.libraryDict = {}
        self.locationDict = {}
        self.logDict = {}
        self.materialSourceDict = {}
        self.migrationsLogDict = {}
        self.museumDict = {}
        self.orderDict = {}
        self.platformDict = {}
        self.projectDict = {}
        self.qualityDict = {}
        self.sampleDict = {}
        self.sequenceDict = {}
        self.sequenceDataReleaseDict = {}
        self.sequencingStrategyDict = {}
        self.sourceDict = {}
        self.sourceSpecimenDict = {}
        self.speciesDict = {}
        self.specimenDict = {}
        self.specimenDataReleaseDict = {}
        self.specimenGeneStatsDict = {}
        self.specimenRawReadsDict = {}
        self.statusDict = {}
        self.iso_country_3166_1Dict = {}
        self.actionDict = loadActionDict(connection, self)
        self.blacklistedReasonDict = loadBlacklistedReasonDict(connection, self)
        self.coordinatesDict = loadCoordinatesDict(connection, self)
        self.dBVersionDict = loadDBVersionDict(connection, self)
        self.dNAVolumeDict = loadDNAVolumeDict(connection, self)
        self.dataReleaseDict = loadDataReleaseDict(connection, self)
        self.dataSourceDict = loadDataSourceDict(connection, self)
        self.extractionTypeDict = loadExtractionTypeDict(connection, self)
        self.orderDict = loadOrderDict(connection, self)
        self.familyDict = loadFamilyDict(connection, self)
        self.geneStatsDict = loadGeneStatsDict(connection, self)
        self.sourceDict = loadSourceDict(connection, self)
        self.genusDict = loadGenusDict(connection, self)
        self.indexesDict = loadIndexesDict(connection, self)
        self.speciesDict = loadSpeciesDict(connection, self)
        self.sourceSpecimenDict = loadSourceSpecimenDict(connection, self)
        self.projectDict = loadProjectDict(connection, self)
        self.materialSourceDict = loadMaterialSourceDict(connection, self)
        self.museumDict = loadMuseumDict(connection, self)
        self.specimenDict = loadSpecimenDict(connection, self)
        self.qualityDict = loadQualityDict(connection, self)
        self.sampleDict = loadSampleDict(connection, self)
        self.statusDict = loadStatusDict(connection, self)
        self.libraryDict = loadLibraryDict(connection, self)
        self.locationDict = loadLocationDict(connection, self)
        self.logDict = loadLogDict(connection, self)
        self.migrationsLogDict = loadMigrationsLogDict(connection, self)
        self.platformDict = loadPlatformDict(connection, self)
        self.sequencingStrategyDict = loadSequencingStrategyDict(connection, self)
        self.sequenceDict = loadSequenceDict(connection, self)
        self.sequenceDataReleaseDict = loadSequenceDataReleaseDict(connection, self)
        self.specimenDataReleaseDict = loadSpecimenDataReleaseDict(connection, self)
        self.specimenGeneStatsDict = loadSpecimenGeneStatsDict(connection, self)
        self.specimenRawReadsDict = loadSpecimenRawReadsDict(connection, self)
        self.iso_country_3166_1Dict = loadIso_country_3166_1Dict(connection, self)

    def __str__(self):
        s = ''
        s = s + 'action: %d\n' % len(self.actionDict)
        s = s + 'blacklistedReason: %d\n' % len(self.blacklistedReasonDict)
        s = s + 'coordinates: %d\n' % len(self.coordinatesDict)
        s = s + 'dBVersion: %d\n' % len(self.dBVersionDict)
        s = s + 'dNAVolume: %d\n' % len(self.dNAVolumeDict)
        s = s + 'dataRelease: %d\n' % len(self.dataReleaseDict)
        s = s + 'dataSource: %d\n' % len(self.dataSourceDict)
        s = s + 'extractionType: %d\n' % len(self.extractionTypeDict)
        s = s + 'order: %d\n' % len(self.orderDict)
        s = s + 'family: %d\n' % len(self.familyDict)
        s = s + 'geneStats: %d\n' % len(self.geneStatsDict)
        s = s + 'source: %d\n' % len(self.sourceDict)
        s = s + 'genus: %d\n' % len(self.genusDict)
        s = s + 'indexes: %d\n' % len(self.indexesDict)
        s = s + 'species: %d\n' % len(self.speciesDict)
        s = s + 'sourceSpecimen: %d\n' % len(self.sourceSpecimenDict)
        s = s + 'project: %d\n' % len(self.projectDict)
        s = s + 'materialSource: %d\n' % len(self.materialSourceDict)
        s = s + 'museum: %d\n' % len(self.museumDict)
        s = s + 'specimen: %d\n' % len(self.specimenDict)
        s = s + 'quality: %d\n' % len(self.qualityDict)
        s = s + 'sample: %d\n' % len(self.sampleDict)
        s = s + 'status: %d\n' % len(self.statusDict)
        s = s + 'library: %d\n' % len(self.libraryDict)
        s = s + 'location: %d\n' % len(self.locationDict)
        s = s + 'log: %d\n' % len(self.logDict)
        s = s + 'migrationsLog: %d\n' % len(self.migrationsLogDict)
        s = s + 'platform: %d\n' % len(self.platformDict)
        s = s + 'sequencingStrategy: %d\n' % len(self.sequencingStrategyDict)
        s = s + 'sequence: %d\n' % len(self.sequenceDict)
        s = s + 'sequenceDataRelease: %d\n' % len(self.sequenceDataReleaseDict)
        s = s + 'specimenDataRelease: %d\n' % len(self.specimenDataReleaseDict)
        s = s + 'specimenGeneStats: %d\n' % len(self.specimenGeneStatsDict)
        s = s + 'specimenRawReads: %d\n' % len(self.specimenRawReadsDict)
        s = s + 'iso_country_3166_1: %d\n' % len(self.iso_country_3166_1Dict)
        return s

