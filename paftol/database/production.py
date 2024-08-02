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


class AdditionalBaitKit(object):

    def __init__(self, additionalBaitKit=None):
        self.idAdditionalBaitKit = None
        self.additionalBaitKit = additionalBaitKit
        # one-to-many
        # fk_AdditionalBaitKit_Sequence: Sequence.idAdditionalBaitKit REFERENCES AdditionalBaitKit(idAdditionalBaitKit)
        self.sequenceAdditionalBaitKitList = []

    def insertIntoDatabase(self, cursor):
        sqlCmd = 'INSERT INTO `AdditionalBaitKit` (`AdditionalBaitKit`) VALUES (%s)'
        l = []
        l.append(self.additionalBaitKit)
        cursor.execute(sqlCmd, tuple(l))


class AlternativeName(object):

    def __init__(self, specimen=None, order=None, family=None, taxonName=None, nameNote=None, dateAdded=None, nameSource=None):
        self.idAlternativeName = None
        self.specimen = specimen
        self.order = order
        self.family = family
        self.taxonName = taxonName
        self.nameNote = nameNote
        self.dateAdded = dateAdded
        self.nameSource = nameSource
        # one-to-many

    def insertIntoDatabase(self, cursor):
        sqlCmd = 'INSERT INTO `AlternativeName` (`idSpecimen`, `Order`, `Family`, `TaxonName`, `NameNote`, `DateAdded`, `idNameSource`) VALUES (%s, %s, %s, %s, %s, %s, %s)'
        l = []
        l.append(None if self.specimen is None else self.specimen.idSpecimen)
        l.append(self.order)
        l.append(self.family)
        l.append(self.taxonName)
        l.append(self.nameNote)
        l.append(self.dateAdded)
        l.append(None if self.nameSource is None else self.nameSource.idNameSource)
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
        # fk_DataSourceAndMethodInRelease_idDataRelease: DataSourceAndMethodInRelease.idDataRelease REFERENCES DataRelease(idDataRelease)
        self.dataSourceAndMethodInReleaseDataReleaseList = []
        # fk_SequenceDataRelease_DataRelease: SequenceDataRelease.idDataRelease REFERENCES DataRelease(idDataRelease)
        self.sequenceDataReleaseDataReleaseList = []
        # fk_SpeciesTree_idDataRelease: SpeciesTree.idDataRelease REFERENCES DataRelease(idDataRelease)
        self.speciesTreeDataReleaseList = []

    def insertIntoDatabase(self, cursor):
        sqlCmd = 'INSERT INTO `DataRelease` (`ReleaseNumber`, `DataRelease`, `TaxonCount`) VALUES (%s, %s, %s)'
        l = []
        l.append(self.releaseNumber)
        l.append(self.dataRelease)
        l.append(self.taxonCount)
        cursor.execute(sqlCmd, tuple(l))


class DataSource(object):

    def __init__(self, dataSource=None, sequenceType=None, isAssembled=None):
        self.idDataSource = None
        self.dataSource = dataSource
        self.sequenceType = sequenceType
        self.isAssembled = isAssembled
        # one-to-many
        # fk_DataSourceAndMethodInRelease_idDataSource: DataSourceAndMethodInRelease.idDataSource REFERENCES DataSource(idDataSource)
        self.dataSourceAndMethodInReleaseDataSourceList = []
        # fk_Specimen_DataSource: Specimen.idDataSource REFERENCES DataSource(idDataSource)
        self.specimenDataSourceList = []

    def insertIntoDatabase(self, cursor):
        sqlCmd = 'INSERT INTO `DataSource` (`DataSource`, `SequenceType`, `IsAssembled`) VALUES (%s, %s, %s)'
        l = []
        l.append(self.dataSource)
        l.append(self.sequenceType)
        l.append(self.isAssembled)
        cursor.execute(sqlCmd, tuple(l))


class DataSourceAndMethodInRelease(object):

    def __init__(self, dataRelease=None, recoveryMethod=None, dataSource=None):
        self.idDataSourceAndMethodInRelease = None
        self.dataRelease = dataRelease
        self.recoveryMethod = recoveryMethod
        self.dataSource = dataSource
        # one-to-many

    def insertIntoDatabase(self, cursor):
        sqlCmd = 'INSERT INTO `DataSourceAndMethodInRelease` (`idDataRelease`, `idRecoveryMethod`, `idDataSource`) VALUES (%s, %s, %s)'
        l = []
        l.append(None if self.dataRelease is None else self.dataRelease.idDataRelease)
        l.append(None if self.recoveryMethod is None else self.recoveryMethod.idRecoveryMethod)
        l.append(None if self.dataSource is None else self.dataSource.idDataSource)
        cursor.execute(sqlCmd, tuple(l))


class DecisionReason(object):

    def __init__(self, decisionReason=None):
        self.idDecisionReason = None
        self.decisionReason = decisionReason
        # one-to-many
        # fk_SequenceDataRelease_DecisionReason: SequenceDataRelease.idDecisionReason REFERENCES DecisionReason(idDecisionReason)
        self.sequenceDataReleaseDecisionReasonList = []

    def insertIntoDatabase(self, cursor):
        sqlCmd = 'INSERT INTO `DecisionReason` (`DecisionReason`) VALUES (%s)'
        l = []
        l.append(self.decisionReason)
        cursor.execute(sqlCmd, tuple(l))


class ExemplarGene(object):

    def __init__(self, ac=None, gn=None, de=None, os=None, url=None):
        self.idExemplarGene = None
        self.ac = ac
        self.gn = gn
        self.de = de
        self.os = os
        self.url = url
        # one-to-many
        # fk_Gene_idExemplarGene: Gene.idExemplarGene REFERENCES ExemplarGene(idExemplarGene)
        self.geneExemplarGeneList = []

    def insertIntoDatabase(self, cursor):
        sqlCmd = 'INSERT INTO `ExemplarGene` (`AC`, `GN`, `DE`, `OS`, `URL`) VALUES (%s, %s, %s, %s, %s)'
        l = []
        l.append(self.ac)
        l.append(self.gn)
        l.append(self.de)
        l.append(self.os)
        l.append(self.url)
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


class FastqStats(object):

    def __init__(self, isTrimmed=None, numReads=None, qual28=None, meanA=None, meanC=None, meanG=None, meanT=None, stddevA=None, stddevC=None, stddevG=None, stddevT=None, meanN=None, stddevN=None, meanAdapterContent=None, maxAdapterContent=None):
        self.idFastqStats = None
        self.isTrimmed = isTrimmed
        self.numReads = numReads
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
        # one-to-many
        # fk_RawFastqFile_idFastqStats: RawFastqFile.idFastqStats REFERENCES FastqStats(idFastqStats)
        self.rawFastqFileFastqStatsList = []
        # fk_TrimmedRawFastqFile_idFastqStats: TrimmedRawFastqFile.idFastqStats REFERENCES FastqStats(idFastqStats)
        self.trimmedRawFastqFileFastqStatsList = []

    def insertIntoDatabase(self, cursor):
        sqlCmd = 'INSERT INTO `FastqStats` (`IsTrimmed`, `NumReads`, `Qual28`, `MeanA`, `MeanC`, `MeanG`, `MeanT`, `StddevA`, `StddevC`, `StddevG`, `StddevT`, `MeanN`, `StddevN`, `MeanAdapterContent`, `MaxAdapterContent`) VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s)'
        l = []
        l.append(self.isTrimmed)
        l.append(self.numReads)
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
        cursor.execute(sqlCmd, tuple(l))


class Gene(object):

    def __init__(self, geneName=None, geneType=None, exemplarGene=None):
        self.idGene = None
        self.geneName = geneName
        self.geneType = geneType
        self.exemplarGene = exemplarGene
        # one-to-many
        # fk_GeneTree_idGene: GeneTree.idGene REFERENCES Gene(idGene)
        self.geneTreeGeneList = []
        # fk_ReferenceTarget_idGene: ReferenceTarget.idGene REFERENCES Gene(idGene)
        self.referenceTargetGeneList = []

    def insertIntoDatabase(self, cursor):
        sqlCmd = 'INSERT INTO `Gene` (`GeneName`, `idGeneType`, `idExemplarGene`) VALUES (%s, %s, %s)'
        l = []
        l.append(self.geneName)
        l.append(None if self.geneType is None else self.geneType.idGeneType)
        l.append(None if self.exemplarGene is None else self.exemplarGene.idExemplarGene)
        cursor.execute(sqlCmd, tuple(l))


class GeneStats(object):

    def __init__(self, internalName=None, exemplarAccession=None, exemplarName=None, exemplarSpecies=None, exemplarHyperlink=None, newickFile=None, newickFilePathName=None, averageContigLength=None, depth=None, averageContigLengthPercentage=None, numSeq=None, numGenera=None, numSpecies=None):
        self.idGeneStats = None
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


class GeneTree(object):

    def __init__(self, gene=None, unAlnFastaFilePathName=None, alnFastaFilePathName=None, newickTree=None, newickFilePathName=None, speciesTree=None):
        self.idGeneTree = None
        self.gene = gene
        self.unAlnFastaFilePathName = unAlnFastaFilePathName
        self.alnFastaFilePathName = alnFastaFilePathName
        self.newickTree = newickTree
        self.newickFilePathName = newickFilePathName
        self.speciesTree = speciesTree
        # one-to-many

    def insertIntoDatabase(self, cursor):
        sqlCmd = 'INSERT INTO `GeneTree` (`idGene`, `UnAlnFastaFilePathName`, `AlnFastaFilePathName`, `NewickTree`, `NewickFilePathName`, `idSpeciesTree`) VALUES (%s, %s, %s, %s, %s, %s)'
        l = []
        l.append(None if self.gene is None else self.gene.idGene)
        l.append(self.unAlnFastaFilePathName)
        l.append(self.alnFastaFilePathName)
        l.append(self.newickTree)
        l.append(self.newickFilePathName)
        l.append(None if self.speciesTree is None else self.speciesTree.idSpeciesTree)
        cursor.execute(sqlCmd, tuple(l))


class GeneType(object):

    def __init__(self, geneTypeName=None):
        self.idGeneType = None
        self.geneTypeName = geneTypeName
        # one-to-many
        # fk_Gene_idGeneType: Gene.idGeneType REFERENCES GeneType(idGeneType)
        self.geneGeneTypeList = []

    def insertIntoDatabase(self, cursor):
        sqlCmd = 'INSERT INTO `GeneType` (`GeneTypeName`) VALUES (%s)'
        l = []
        l.append(self.geneTypeName)
        cursor.execute(sqlCmd, tuple(l))


class Herbarium(object):

    def __init__(self, herbariumCode=None, herbariumName=None, herbariumUrL=None):
        self.idHerbarium = None
        self.herbariumCode = herbariumCode
        self.herbariumName = herbariumName
        self.herbariumUrL = herbariumUrL
        # one-to-many
        # fk_Specimen_Herbarium: Specimen.idHerbarium REFERENCES Herbarium(idHerbarium)
        self.specimenHerbariumList = []

    def insertIntoDatabase(self, cursor):
        sqlCmd = 'INSERT INTO `Herbarium` (`HerbariumCode`, `HerbariumName`, `HerbariumURL`) VALUES (%s, %s, %s)'
        l = []
        l.append(self.herbariumCode)
        l.append(self.herbariumName)
        l.append(self.herbariumUrL)
        cursor.execute(sqlCmd, tuple(l))


class ISOCountry(object):

    def __init__(self, country=None, alpha2Code=None, alpha3Code=None, countryCode=None, isoCode=None, region=None, subRegion=None, intermediateRegion=None, regionCode=None, subRegionCode=None, intermediateRegionCode=None):
        self.idIsoCountry = None
        self.country = country
        self.alpha2Code = alpha2Code
        self.alpha3Code = alpha3Code
        self.countryCode = countryCode
        self.isoCode = isoCode
        self.region = region
        self.subRegion = subRegion
        self.intermediateRegion = intermediateRegion
        self.regionCode = regionCode
        self.subRegionCode = subRegionCode
        self.intermediateRegionCode = intermediateRegionCode
        # one-to-many
        # fk_Specimen_ISOCountry: Specimen.idISOCountry REFERENCES ISOCountry(idISOCountry)
        self.specimenIsoCountryList = []

    def insertIntoDatabase(self, cursor):
        sqlCmd = 'INSERT INTO `ISOCountry` (`Country`, `Alpha2Code`, `Alpha3Code`, `CountryCode`, `ISOCode`, `Region`, `SubRegion`, `IntermediateRegion`, `RegionCode`, `SubRegionCode`, `IntermediateRegionCode`) VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s)'
        l = []
        l.append(self.country)
        l.append(self.alpha2Code)
        l.append(self.alpha3Code)
        l.append(self.countryCode)
        l.append(self.isoCode)
        l.append(self.region)
        l.append(self.subRegion)
        l.append(self.intermediateRegion)
        l.append(self.regionCode)
        l.append(self.subRegionCode)
        l.append(self.intermediateRegionCode)
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

    def __init__(self, sample=None, libConcentration=None, libQuality=None, libTapeStation=None, sonication=None, plateNumber=None, coordinate=None, description=None, status=None, indexes=None, generateLibrary=None, pcrCycles=None):
        self.idLibrary = None
        self.sample = sample
        self.libConcentration = libConcentration
        self.libQuality = libQuality
        self.libTapeStation = libTapeStation
        self.sonication = sonication
        self.plateNumber = plateNumber
        self.coordinate = coordinate
        self.description = description
        self.status = status
        self.indexes = indexes
        self.generateLibrary = generateLibrary
        self.pcrCycles = pcrCycles
        # one-to-many
        # LibraryLink: Sequence.idLibrary REFERENCES Library(idLibrary)
        self.sequenceLibraryList = []

    def insertIntoDatabase(self, cursor):
        sqlCmd = 'INSERT INTO `Library` (`idSample`, `LibConcentration`, `LibQuality`, `LibTapeStation`, `Sonication`, `PlateNumber`, `idCoordinate`, `Description`, `idStatus`, `idIndexes`, `GenerateLibrary`, `PCRCycles`) VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s)'
        l = []
        l.append(None if self.sample is None else self.sample.idSample)
        l.append(self.libConcentration)
        l.append(self.libQuality)
        l.append(self.libTapeStation)
        l.append(self.sonication)
        l.append(self.plateNumber)
        l.append(None if self.coordinate is None else self.coordinate.idCoordinate)
        l.append(self.description)
        l.append(None if self.status is None else self.status.idStatus)
        l.append(None if self.indexes is None else self.indexes.idIndexes)
        l.append(self.generateLibrary)
        l.append(self.pcrCycles)
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


class MergedSequence(object):

    def __init__(self, sequence=None, sequenceOrigin=None):
        self.idMergedSequence = None
        self.sequence = sequence
        self.sequenceOrigin = sequenceOrigin
        # one-to-many

    def insertIntoDatabase(self, cursor):
        sqlCmd = 'INSERT INTO `MergedSequence` (`idSequence`, `idSequenceOrigin`) VALUES (%s, %s)'
        l = []
        l.append(None if self.sequence is None else self.sequence.idSequence)
        l.append(None if self.sequenceOrigin is None else self.sequenceOrigin.idSequence)
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


class NameSource(object):

    def __init__(self, nameSource=None):
        self.idNameSource = None
        self.nameSource = nameSource
        # one-to-many
        # AlternativeName_ibfk_2: AlternativeName.idNameSource REFERENCES NameSource(idNameSource)
        self.alternativeNameNameSourceList = []

    def insertIntoDatabase(self, cursor):
        sqlCmd = 'INSERT INTO `NameSource` (`NameSource`) VALUES (%s)'
        l = []
        l.append(self.nameSource)
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

    def __init__(self, project=None):
        self.idProject = None
        self.project = project
        # one-to-many
        # fk_idProject: Specimen.idProject REFERENCES Project(idProject)
        self.specimenProjectList = []

    def insertIntoDatabase(self, cursor):
        sqlCmd = 'INSERT INTO `Project` (`Project`) VALUES (%s)'
        l = []
        l.append(self.project)
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


class RawFastaFile(object):

    def __init__(self, filePathName=None, uncompressedFileMd5sum=None, numRecords=None, sequence=None):
        self.idRawFastaFile = None
        self.filePathName = filePathName
        self.uncompressedFileMd5sum = uncompressedFileMd5sum
        self.numRecords = numRecords
        self.sequence = sequence
        # one-to-many

    def insertIntoDatabase(self, cursor):
        sqlCmd = 'INSERT INTO `RawFastaFile` (`FilePathName`, `UncompressedFileMd5sum`, `NumRecords`, `idSequence`) VALUES (%s, %s, %s, %s)'
        l = []
        l.append(self.filePathName)
        l.append(self.uncompressedFileMd5sum)
        l.append(self.numRecords)
        l.append(None if self.sequence is None else self.sequence.idSequence)
        cursor.execute(sqlCmd, tuple(l))


class RawFastqFile(object):

    def __init__(self, filePathName=None, uncompressedFileMd5sum=None, orientation=None, sequence=None, fastqStats=None):
        self.idRawFastqFile = None
        self.filePathName = filePathName
        self.uncompressedFileMd5sum = uncompressedFileMd5sum
        self.orientation = orientation
        self.sequence = sequence
        self.fastqStats = fastqStats
        # one-to-many
        # fk_TrimmedRawFastqFile_idRawFastqFile: TrimmedRawFastqFile.idRawFastqFile REFERENCES RawFastqFile(idRawFastqFile)
        self.trimmedRawFastqFileRawFastqFileList = []

    def insertIntoDatabase(self, cursor):
        sqlCmd = 'INSERT INTO `RawFastqFile` (`FilePathName`, `UncompressedFileMd5sum`, `Orientation`, `idSequence`, `idFastqStats`) VALUES (%s, %s, %s, %s, %s)'
        l = []
        l.append(self.filePathName)
        l.append(self.uncompressedFileMd5sum)
        l.append(self.orientation)
        l.append(None if self.sequence is None else self.sequence.idSequence)
        l.append(None if self.fastqStats is None else self.fastqStats.idFastqStats)
        cursor.execute(sqlCmd, tuple(l))


class RecoveredContig(object):

    def __init__(self, sequenceRecovery=None, seqlength=None, contigMd5sum=None, representativeReferenceTarget=None):
        self.idRecoveredContig = None
        self.sequenceRecovery = sequenceRecovery
        self.seqlength = seqlength
        self.contigMd5sum = contigMd5sum
        self.representativeReferenceTarget = representativeReferenceTarget
        # one-to-many

    def insertIntoDatabase(self, cursor):
        sqlCmd = 'INSERT INTO `RecoveredContig` (`idSequenceRecovery`, `Seqlength`, `ContigMd5sum`, `idRepresentativeReferenceTarget`) VALUES (%s, %s, %s, %s)'
        l = []
        l.append(None if self.sequenceRecovery is None else self.sequenceRecovery.idSequenceRecovery)
        l.append(self.seqlength)
        l.append(self.contigMd5sum)
        l.append(None if self.representativeReferenceTarget is None else self.representativeReferenceTarget.idReferenceTarget)
        cursor.execute(sqlCmd, tuple(l))


class RecoveryMethod(object):

    def __init__(self, software=None, softwareVersion=None, softwareParameters=None, targetSet=None):
        self.idRecoveryMethod = None
        self.software = software
        self.softwareVersion = softwareVersion
        self.softwareParameters = softwareParameters
        self.targetSet = targetSet
        # one-to-many
        # fk_DataSourceAndMethodInRelease_idRecoveryMethod: DataSourceAndMethodInRelease.idRecoveryMethod REFERENCES RecoveryMethod(idRecoveryMethod)
        self.dataSourceAndMethodInReleaseRecoveryMethodList = []
        # fk_SequenceRecovery_idRecoveryMethod: SequenceRecovery.idRecoveryMethod REFERENCES RecoveryMethod(idRecoveryMethod)
        self.sequenceRecoveryRecoveryMethodList = []

    def insertIntoDatabase(self, cursor):
        sqlCmd = 'INSERT INTO `RecoveryMethod` (`Software`, `SoftwareVersion`, `SoftwareParameters`, `idTargetSet`) VALUES (%s, %s, %s, %s)'
        l = []
        l.append(self.software)
        l.append(self.softwareVersion)
        l.append(self.softwareParameters)
        l.append(None if self.targetSet is None else self.targetSet.idTargetSet)
        cursor.execute(sqlCmd, tuple(l))


class ReferenceTarget(object):

    def __init__(self, gene=None, organism=None, targetLength=None, targetSet=None):
        self.idReferenceTarget = None
        self.gene = gene
        self.organism = organism
        self.targetLength = targetLength
        self.targetSet = targetSet
        # one-to-many
        # fk_RecoveredContig_idRepresentativeReferenceTarget: RecoveredContig.idRepresentativeReferenceTarget REFERENCES ReferenceTarget(idRepresentativeReferenceTarget)
        self.recoveredContigRepresentativeReferenceTargetList = []

    def insertIntoDatabase(self, cursor):
        sqlCmd = 'INSERT INTO `ReferenceTarget` (`idGene`, `Organism`, `TargetLength`, `idTargetSet`) VALUES (%s, %s, %s, %s)'
        l = []
        l.append(None if self.gene is None else self.gene.idGene)
        l.append(self.organism)
        l.append(self.targetLength)
        l.append(None if self.targetSet is None else self.targetSet.idTargetSet)
        cursor.execute(sqlCmd, tuple(l))


class Sample(object):

    def __init__(self, description=None, action=None, extractionType=None, quality=None, sampleConcentration=None, gelImage=None, sampleTapeStation=None, dnaVolume=None, enaSampleNum=None, secEnaSampleNum=None, specimen=None):
        self.idSample = None
        self.description = description
        self.action = action
        self.extractionType = extractionType
        self.quality = quality
        self.sampleConcentration = sampleConcentration
        self.gelImage = gelImage
        self.sampleTapeStation = sampleTapeStation
        self.dnaVolume = dnaVolume
        self.enaSampleNum = enaSampleNum
        self.secEnaSampleNum = secEnaSampleNum
        self.specimen = specimen
        # one-to-many
        # SampleLink: Library.idSample REFERENCES Sample(idSample)
        self.librarySampleList = []

    def insertIntoDatabase(self, cursor):
        sqlCmd = 'INSERT INTO `Sample` (`Description`, `idAction`, `idExtractionType`, `idQuality`, `SampleConcentration`, `GelImage`, `SampleTapeStation`, `idDNAVolume`, `ENASampleNum`, `SecENASampleNum`, `idSpecimen`) VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s)'
        l = []
        l.append(self.description)
        l.append(None if self.action is None else self.action.idAction)
        l.append(None if self.extractionType is None else self.extractionType.idExtractionType)
        l.append(None if self.quality is None else self.quality.idQuality)
        l.append(self.sampleConcentration)
        l.append(self.gelImage)
        l.append(self.sampleTapeStation)
        l.append(None if self.dnaVolume is None else self.dnaVolume.idDnaVolume)
        l.append(self.enaSampleNum)
        l.append(self.secEnaSampleNum)
        l.append(None if self.specimen is None else self.specimen.idSpecimen)
        cursor.execute(sqlCmd, tuple(l))


class Sequence(object):

    def __init__(self, externalSequenceId=None, hasDuplicate=None, isMerged=None, library=None, platform=None, location=None, sequencingRun=None, numInferredCds=None, medianHybpiperCdsLength=None, status=None, hybridisationPool=None, r2FastqFile=None, r1FastqFile=None, blacklisted=None, blacklistedReason=None, blacklistedPerson=None, sequencingStrategy=None, enaExpNumber=None, enaRunNumber=None, suspiciousPlacement=None, additionalBaitKit=None):
        self.idSequence = None
        self.externalSequenceId = externalSequenceId
        self.hasDuplicate = hasDuplicate
        self.isMerged = isMerged
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
        self.blacklistedPerson = blacklistedPerson
        self.sequencingStrategy = sequencingStrategy
        self.enaExpNumber = enaExpNumber
        self.enaRunNumber = enaRunNumber
        self.suspiciousPlacement = suspiciousPlacement
        self.additionalBaitKit = additionalBaitKit
        # one-to-many
        # fk_MergedSequence_idSequence: MergedSequence.idSequence REFERENCES Sequence(idSequence)
        self.mergedSequenceSequenceList = []
        # fk_MergedSequence_idSequenceOrigin: MergedSequence.idSequenceOrigin REFERENCES Sequence(idSequenceOrigin)
        self.mergedSequenceSequenceOriginList = []
        # fk_RawFastaFile_idSequence: RawFastaFile.idSequence REFERENCES Sequence(idSequence)
        self.rawFastaFileSequenceList = []
        # fk_RawFastqFile_idSequence: RawFastqFile.idSequence REFERENCES Sequence(idSequence)
        self.rawFastqFileSequenceList = []
        # fk_SequenceDataRelease_Sequence: SequenceDataRelease.idSequence REFERENCES Sequence(idSequence)
        self.sequenceDataReleaseSequenceList = []
        # fk_SequenceGeneStats_Sequence: SequenceGeneStats.idSequence REFERENCES Sequence(idSequence)
        self.sequenceGeneStatsSequenceList = []
        # fk_SequenceRawReads_Sequence: SequenceRawReads.idSequence REFERENCES Sequence(idSequence)
        self.sequenceRawReadsSequenceList = []
        # fk_SequenceRecovery_idSequence: SequenceRecovery.idSequence REFERENCES Sequence(idSequence)
        self.sequenceRecoverySequenceList = []

    def insertIntoDatabase(self, cursor):
        sqlCmd = 'INSERT INTO `Sequence` (`ExternalSequenceID`, `HasDuplicate`, `IsMerged`, `idLibrary`, `idPlatform`, `idLocation`, `SequencingRun`, `NumInferredCds`, `MedianHybpiperCdsLength`, `idStatus`, `HybridisationPool`, `R2FastqFile`, `R1FastqFile`, `Blacklisted`, `idBlacklistedReason`, `BlacklistedPerson`, `idSequencingStrategy`, `ENAExpNumber`, `ENARunNumber`, `SuspiciousPlacement`, `idAdditionalBaitKit`) VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s)'
        l = []
        l.append(self.externalSequenceId)
        l.append(self.hasDuplicate)
        l.append(self.isMerged)
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
        l.append(self.blacklistedPerson)
        l.append(None if self.sequencingStrategy is None else self.sequencingStrategy.idSequencingStrategy)
        l.append(self.enaExpNumber)
        l.append(self.enaRunNumber)
        l.append(self.suspiciousPlacement)
        l.append(None if self.additionalBaitKit is None else self.additionalBaitKit.idAdditionalBaitKit)
        cursor.execute(sqlCmd, tuple(l))


class SequenceDataRelease(object):

    def __init__(self, sequence=None, dataRelease=None, barcodeValidation=None, phylogeneticValidation=None, validationResult=None, decisionReason=None, validationComments=None):
        self.idSequenceDataRelease = None
        self.sequence = sequence
        self.dataRelease = dataRelease
        self.barcodeValidation = barcodeValidation
        self.phylogeneticValidation = phylogeneticValidation
        self.validationResult = validationResult
        self.decisionReason = decisionReason
        self.validationComments = validationComments
        # one-to-many

    def insertIntoDatabase(self, cursor):
        sqlCmd = 'INSERT INTO `SequenceDataRelease` (`idSequence`, `idDataRelease`, `idBarcodeValidation`, `idPhylogeneticValidation`, `idValidationResult`, `idDecisionReason`, `ValidationComments`) VALUES (%s, %s, %s, %s, %s, %s, %s)'
        l = []
        l.append(None if self.sequence is None else self.sequence.idSequence)
        l.append(None if self.dataRelease is None else self.dataRelease.idDataRelease)
        l.append(None if self.barcodeValidation is None else self.barcodeValidation.idTestResult)
        l.append(None if self.phylogeneticValidation is None else self.phylogeneticValidation.idTestResult)
        l.append(None if self.validationResult is None else self.validationResult.idValidationResult)
        l.append(None if self.decisionReason is None else self.decisionReason.idDecisionReason)
        l.append(self.validationComments)
        cursor.execute(sqlCmd, tuple(l))


class SequenceGeneStats(object):

    def __init__(self, sequence=None, numRecoveredGenes=None, sumContigLength=None):
        self.idSequenceGeneStats = None
        self.sequence = sequence
        self.numRecoveredGenes = numRecoveredGenes
        self.sumContigLength = sumContigLength
        # one-to-many

    def insertIntoDatabase(self, cursor):
        sqlCmd = 'INSERT INTO `SequenceGeneStats` (`idSequence`, `NumRecoveredGenes`, `SumContigLength`) VALUES (%s, %s, %s)'
        l = []
        l.append(None if self.sequence is None else self.sequence.idSequence)
        l.append(self.numRecoveredGenes)
        l.append(self.sumContigLength)
        cursor.execute(sqlCmd, tuple(l))


class SequenceRawReads(object):

    def __init__(self, sequence=None, numReads=None):
        self.idSequenceRawReads = None
        self.sequence = sequence
        self.numReads = numReads
        # one-to-many

    def insertIntoDatabase(self, cursor):
        sqlCmd = 'INSERT INTO `SequenceRawReads` (`idSequence`, `NumReads`) VALUES (%s, %s)'
        l = []
        l.append(None if self.sequence is None else self.sequence.idSequence)
        l.append(self.numReads)
        cursor.execute(sqlCmd, tuple(l))


class SequenceRecovery(object):

    def __init__(self, sequence=None, contigFastaFilePathName=None, contigFastaFileMd5sum=None, numMappedReads=None, numUnmappedReads=None, softwareVersion=None, cmdLine=None, numRecoveredContigsCheck=None, recoveryMethod=None):
        self.idSequenceRecovery = None
        self.sequence = sequence
        self.contigFastaFilePathName = contigFastaFilePathName
        self.contigFastaFileMd5sum = contigFastaFileMd5sum
        self.numMappedReads = numMappedReads
        self.numUnmappedReads = numUnmappedReads
        self.softwareVersion = softwareVersion
        self.cmdLine = cmdLine
        self.numRecoveredContigsCheck = numRecoveredContigsCheck
        self.recoveryMethod = recoveryMethod
        # one-to-many
        # fk_RecoveredContig_idSequenceRecovery: RecoveredContig.idSequenceRecovery REFERENCES SequenceRecovery(idSequenceRecovery)
        self.recoveredContigSequenceRecoveryList = []
        # fk_TrimmedRawFastqFile_idSequenceRecovery: TrimmedRawFastqFile.idSequenceRecovery REFERENCES SequenceRecovery(idSequenceRecovery)
        self.trimmedRawFastqFileSequenceRecoveryList = []

    def insertIntoDatabase(self, cursor):
        sqlCmd = 'INSERT INTO `SequenceRecovery` (`idSequence`, `ContigFastaFilePathName`, `ContigFastaFileMd5sum`, `NumMappedReads`, `NumUnmappedReads`, `SoftwareVersion`, `CmdLine`, `NumRecoveredContigsCheck`, `idRecoveryMethod`) VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s)'
        l = []
        l.append(None if self.sequence is None else self.sequence.idSequence)
        l.append(self.contigFastaFilePathName)
        l.append(self.contigFastaFileMd5sum)
        l.append(self.numMappedReads)
        l.append(self.numUnmappedReads)
        l.append(self.softwareVersion)
        l.append(self.cmdLine)
        l.append(self.numRecoveredContigsCheck)
        l.append(None if self.recoveryMethod is None else self.recoveryMethod.idRecoveryMethod)
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


class SpeciesTree(object):

    def __init__(self, newickTree=None, newickFilePathName=None, pipelineVersion=None, pipelineCmdLine=None, dataRelease=None):
        self.idSpeciesTree = None
        self.newickTree = newickTree
        self.newickFilePathName = newickFilePathName
        self.pipelineVersion = pipelineVersion
        self.pipelineCmdLine = pipelineCmdLine
        self.dataRelease = dataRelease
        # one-to-many
        # fk_GeneTree_idSpeciesTree: GeneTree.idSpeciesTree REFERENCES SpeciesTree(idSpeciesTree)
        self.geneTreeSpeciesTreeList = []

    def insertIntoDatabase(self, cursor):
        sqlCmd = 'INSERT INTO `SpeciesTree` (`NewickTree`, `NewickFilePathName`, `PipelineVersion`, `PipelineCmdLine`, `idDataRelease`) VALUES (%s, %s, %s, %s, %s)'
        l = []
        l.append(self.newickTree)
        l.append(self.newickFilePathName)
        l.append(self.pipelineVersion)
        l.append(self.pipelineCmdLine)
        l.append(None if self.dataRelease is None else self.dataRelease.idDataRelease)
        cursor.execute(sqlCmd, tuple(l))


class Specimen(object):

    def __init__(self, POWO=None, bankId=None, lcd=None, msb=None, collector=None, collectorNo=None, voucherNo=None, museumBarcode=None, oldSpeciesName=None, sourceSpecimen=None, project=None, isoCountry=None, materialSource=None, ageOfMaterial=None, herbarium=None, specimenReference=None, herbcatUrL=None, blacklisted=None, holdDate=None, dataSource=None):
        self.idSpecimen = None
        self.POWO = POWO
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
        self.isoCountry = isoCountry
        self.materialSource = materialSource
        self.ageOfMaterial = ageOfMaterial
        self.herbarium = herbarium
        self.specimenReference = specimenReference
        self.herbcatUrL = herbcatUrL
        self.blacklisted = blacklisted
        self.holdDate = holdDate
        self.dataSource = dataSource
        # one-to-many
        # AlternativeName_ibfk_1: AlternativeName.idSpecimen REFERENCES Specimen(idSpecimen)
        self.alternativeNameSpecimenList = []
        # fk_Sample_Specimen: Sample.idSpecimen REFERENCES Specimen(idSpecimen)
        self.sampleSpecimenList = []
        # no python attribute: SourceSpecimen

    def insertIntoDatabase(self, cursor):
        sqlCmd = 'INSERT INTO `Specimen` (`POWOId`, `BankID`, `LCD`, `MSB`, `Collector`, `CollectorNo`, `VoucherNo`, `MuseumBarcode`, `OldSpeciesName`, `idSourceSpecimen`, `idProject`, `idISOCountry`, `idMaterialSource`, `AgeOfMaterial`, `idHerbarium`, `SpecimenReference`, `HerbcatURL`, `Blacklisted`, `HoldDate`, `idDataSource`) VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s)'
        l = []
        l.append(None if self.POWO is None else self.POWO.powoId)
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
        l.append(None if self.isoCountry is None else self.isoCountry.idIsoCountry)
        l.append(None if self.materialSource is None else self.materialSource.idMaterialSource)
        l.append(self.ageOfMaterial)
        l.append(None if self.herbarium is None else self.herbarium.idHerbarium)
        l.append(self.specimenReference)
        l.append(self.herbcatUrL)
        l.append(self.blacklisted)
        l.append(self.holdDate)
        l.append(None if self.dataSource is None else self.dataSource.idDataSource)
        cursor.execute(sqlCmd, tuple(l))


class SpecimenGeneStats(object):

    def __init__(self, numGene=None, recoveredLength=None):
        self.idSpecimen = None
        self.numGene = numGene
        self.recoveredLength = recoveredLength
        # one-to-many

    def insertIntoDatabase(self, cursor):
        sqlCmd = 'INSERT INTO `SpecimenGeneStats` (`NumGene`, `RecoveredLength`) VALUES (%s, %s)'
        l = []
        l.append(self.numGene)
        l.append(self.recoveredLength)
        cursor.execute(sqlCmd, tuple(l))


class SpecimenRawReads(object):

    def __init__(self, idSpecimen=None, numReads=None, seqPlatform=None, enaExpNum=None, enaRunNum=None):
        self.idSpecimenRawReads = None
        self.idSpecimen = idSpecimen
        self.numReads = numReads
        self.seqPlatform = seqPlatform
        self.enaExpNum = enaExpNum
        self.enaRunNum = enaRunNum
        # one-to-many

    def insertIntoDatabase(self, cursor):
        sqlCmd = 'INSERT INTO `SpecimenRawReads` (`idSpecimen`, `NumReads`, `SeqPlatform`, `ENAExpNum`, `ENARunNum`) VALUES (%s, %s, %s, %s, %s)'
        l = []
        l.append(self.idSpecimen)
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


class TargetSet(object):

    def __init__(self, targetSetName=None, targetsFastaFile=None, targetsFastaFileMd5sum=None, numTargetSequences=None):
        self.idTargetSet = None
        self.targetSetName = targetSetName
        self.targetsFastaFile = targetsFastaFile
        self.targetsFastaFileMd5sum = targetsFastaFileMd5sum
        self.numTargetSequences = numTargetSequences
        # one-to-many
        # fk_RecoveryMethod_idTargetSet: RecoveryMethod.idTargetSet REFERENCES TargetSet(idTargetSet)
        self.recoveryMethodTargetSetList = []
        # fk_ReferenceTarget_idTargetSet: ReferenceTarget.idTargetSet REFERENCES TargetSet(idTargetSet)
        self.referenceTargetTargetSetList = []

    def insertIntoDatabase(self, cursor):
        sqlCmd = 'INSERT INTO `TargetSet` (`TargetSetName`, `TargetsFastaFile`, `TargetsFastaFileMd5sum`, `NumTargetSequences`) VALUES (%s, %s, %s, %s)'
        l = []
        l.append(self.targetSetName)
        l.append(self.targetsFastaFile)
        l.append(self.targetsFastaFileMd5sum)
        l.append(self.numTargetSequences)
        cursor.execute(sqlCmd, tuple(l))


class TestResult(object):

    def __init__(self, testResult=None):
        self.idTestResult = None
        self.testResult = testResult
        # one-to-many
        # fk_SequenceDataRelease_TestResult_barcode: SequenceDataRelease.idBarcodeValidation REFERENCES TestResult(idBarcodeValidation)
        self.sequenceDataReleaseBarcodeValidationList = []
        # fk_SequenceDataRelease_TestResult_phylogenetic: SequenceDataRelease.idPhylogeneticValidation REFERENCES TestResult(idPhylogeneticValidation)
        self.sequenceDataReleasePhylogeneticValidationList = []

    def insertIntoDatabase(self, cursor):
        sqlCmd = 'INSERT INTO `TestResult` (`TestResult`) VALUES (%s)'
        l = []
        l.append(self.testResult)
        cursor.execute(sqlCmd, tuple(l))


class TrimmedRawFastqFile(object):

    def __init__(self, rawFastqFile=None, sequenceRecovery=None, fastqStats=None):
        self.idTrimmedRawFastqFile = None
        self.rawFastqFile = rawFastqFile
        self.sequenceRecovery = sequenceRecovery
        self.fastqStats = fastqStats
        # one-to-many

    def insertIntoDatabase(self, cursor):
        sqlCmd = 'INSERT INTO `TrimmedRawFastqFile` (`idRawFastqFile`, `idSequenceRecovery`, `idFastqStats`) VALUES (%s, %s, %s)'
        l = []
        l.append(None if self.rawFastqFile is None else self.rawFastqFile.idRawFastqFile)
        l.append(None if self.sequenceRecovery is None else self.sequenceRecovery.idSequenceRecovery)
        l.append(None if self.fastqStats is None else self.fastqStats.idFastqStats)
        cursor.execute(sqlCmd, tuple(l))


class ValidationResult(object):

    def __init__(self, validationResult=None):
        self.idValidationResult = None
        self.validationResult = validationResult
        # one-to-many
        # fk_SequenceDataRelease_ValidationResult: SequenceDataRelease.idValidationResult REFERENCES ValidationResult(idValidationResult)
        self.sequenceDataReleaseValidationResultList = []

    def insertIntoDatabase(self, cursor):
        sqlCmd = 'INSERT INTO `ValidationResult` (`ValidationResult`) VALUES (%s)'
        l = []
        l.append(self.validationResult)
        cursor.execute(sqlCmd, tuple(l))


class WCVPFamily(object):

    def __init__(self, family=None, ipniId=None, acceptedIpniId=None, apgivId=None, checklistDb=None, peerReviewed=None, taxonStatusId=None, apgTaxonRemarks=None, includedInWcvP=None, order=None):
        self.idFamily = None
        self.family = family
        self.ipniId = ipniId
        self.acceptedIpniId = acceptedIpniId
        self.apgivId = apgivId
        self.checklistDb = checklistDb
        self.peerReviewed = peerReviewed
        self.taxonStatusId = taxonStatusId
        self.apgTaxonRemarks = apgTaxonRemarks
        self.includedInWcvP = includedInWcvP
        self.order = order
        # one-to-many
        # fk_Name_Family: WCVPName.idFamily REFERENCES WCVPFamily(idFamily)
        self.wcvpNameFamilyList = []

    def insertIntoDatabase(self, cursor):
        sqlCmd = 'INSERT INTO `WCVPFamily` (`Family`, `IPNIId`, `AcceptedIPNIId`, `APGIVId`, `ChecklistDB`, `PeerReviewed`, `TaxonStatusId`, `APGTaxonRemarks`, `IncludedInWCVP`, `idOrder`) VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s, %s)'
        l = []
        l.append(self.family)
        l.append(self.ipniId)
        l.append(self.acceptedIpniId)
        l.append(self.apgivId)
        l.append(self.checklistDb)
        l.append(self.peerReviewed)
        l.append(self.taxonStatusId)
        l.append(self.apgTaxonRemarks)
        l.append(self.includedInWcvP)
        l.append(None if self.order is None else self.order.idOrder)
        cursor.execute(sqlCmd, tuple(l))


class WCVPHigherRank(object):

    def __init__(self, kingdom=None, phylum=None, className=None, order=None, checklistDb=None, peerReviewed=None, ipniId=None, taxonStatusId=None, acceptedIpniId=None, apgivId=None, apgivOrderId=None, apgTaxonRemarks=None, includedInWcvP=None, subclass=None, isAngiosperm=None):
        self.family = None
        self.kingdom = kingdom
        self.phylum = phylum
        self.className = className
        self.order = order
        self.checklistDb = checklistDb
        self.peerReviewed = peerReviewed
        self.ipniId = ipniId
        self.taxonStatusId = taxonStatusId
        self.acceptedIpniId = acceptedIpniId
        self.apgivId = apgivId
        self.apgivOrderId = apgivOrderId
        self.apgTaxonRemarks = apgTaxonRemarks
        self.includedInWcvP = includedInWcvP
        self.subclass = subclass
        self.isAngiosperm = isAngiosperm
        # one-to-many

    def insertIntoDatabase(self, cursor):
        sqlCmd = 'INSERT INTO `WCVPHigherRank` (`Kingdom`, `Phylum`, `Class`, `Order`, `ChecklistDB`, `PeerReviewed`, `IPNIId`, `TaxonStatusId`, `AcceptedIPNIId`, `APGIVId`, `APGIVOrderId`, `APGTaxonRemarks`, `IncludedInWCVP`, `Subclass`, `isAngiosperm`) VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s)'
        l = []
        l.append(self.kingdom)
        l.append(self.phylum)
        l.append(self.className)
        l.append(self.order)
        l.append(self.checklistDb)
        l.append(self.peerReviewed)
        l.append(self.ipniId)
        l.append(self.taxonStatusId)
        l.append(self.acceptedIpniId)
        l.append(self.apgivId)
        l.append(self.apgivOrderId)
        l.append(self.apgTaxonRemarks)
        l.append(self.includedInWcvP)
        l.append(self.subclass)
        l.append(self.isAngiosperm)
        cursor.execute(sqlCmd, tuple(l))


class WCVPLoad(object):

    def __init__(self, ipniId=None, taxonRank=None, taxonStatus=None, family=None, genusHybrid=None, genus=None, speciesHybrid=None, species=None, infraspecificRank=None, infraspecies=None, parentheticalAuthor=None, primaryAuthor=None, publicationAuthor=None, placeOfPublication=None, volumeAndPage=None, firstPublished=None, nomenclaturalRemarks=None, geographicArea=None, lifeformDescription=None, climateDescription=None, taxonName=None, taxonAuthors=None, acceptedPlantName=None, basionymPlantNameId=None, replacedSynonymAuthor=None, homotypicSynonym=None, parentPlantName=None, powoId=None, hybridFormula=None, reviewed=None):
        self.plantNameId = None
        self.ipniId = ipniId
        self.taxonRank = taxonRank
        self.taxonStatus = taxonStatus
        self.family = family
        self.genusHybrid = genusHybrid
        self.genus = genus
        self.speciesHybrid = speciesHybrid
        self.species = species
        self.infraspecificRank = infraspecificRank
        self.infraspecies = infraspecies
        self.parentheticalAuthor = parentheticalAuthor
        self.primaryAuthor = primaryAuthor
        self.publicationAuthor = publicationAuthor
        self.placeOfPublication = placeOfPublication
        self.volumeAndPage = volumeAndPage
        self.firstPublished = firstPublished
        self.nomenclaturalRemarks = nomenclaturalRemarks
        self.geographicArea = geographicArea
        self.lifeformDescription = lifeformDescription
        self.climateDescription = climateDescription
        self.taxonName = taxonName
        self.taxonAuthors = taxonAuthors
        self.acceptedPlantName = acceptedPlantName
        self.basionymPlantNameId = basionymPlantNameId
        self.replacedSynonymAuthor = replacedSynonymAuthor
        self.homotypicSynonym = homotypicSynonym
        self.parentPlantName = parentPlantName
        self.powoId = powoId
        self.hybridFormula = hybridFormula
        self.reviewed = reviewed
        # one-to-many

    def insertIntoDatabase(self, cursor):
        sqlCmd = 'INSERT INTO `WCVPLoad` (`IPNIId`, `TaxonRank`, `TaxonStatus`, `Family`, `GenusHybrid`, `Genus`, `SpeciesHybrid`, `Species`, `InfraspecificRank`, `Infraspecies`, `ParentheticalAuthor`, `PrimaryAuthor`, `PublicationAuthor`, `PlaceOfPublication`, `VolumeAndPage`, `FirstPublished`, `NomenclaturalRemarks`, `GeographicArea`, `LifeformDescription`, `ClimateDescription`, `TaxonName`, `TaxonAuthors`, `AcceptedPlantName`, `BasionymPlantNameId`, `ReplacedSynonymAuthor`, `HomotypicSynonym`, `ParentPlantName`, `POWOId`, `HybridFormula`, `Reviewed`) VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s)'
        l = []
        l.append(self.ipniId)
        l.append(self.taxonRank)
        l.append(self.taxonStatus)
        l.append(self.family)
        l.append(self.genusHybrid)
        l.append(self.genus)
        l.append(self.speciesHybrid)
        l.append(self.species)
        l.append(self.infraspecificRank)
        l.append(self.infraspecies)
        l.append(self.parentheticalAuthor)
        l.append(self.primaryAuthor)
        l.append(self.publicationAuthor)
        l.append(self.placeOfPublication)
        l.append(self.volumeAndPage)
        l.append(self.firstPublished)
        l.append(self.nomenclaturalRemarks)
        l.append(self.geographicArea)
        l.append(self.lifeformDescription)
        l.append(self.climateDescription)
        l.append(self.taxonName)
        l.append(self.taxonAuthors)
        l.append(self.acceptedPlantName)
        l.append(self.basionymPlantNameId)
        l.append(self.replacedSynonymAuthor)
        l.append(self.homotypicSynonym)
        l.append(self.parentPlantName)
        l.append(self.powoId)
        l.append(self.hybridFormula)
        l.append(self.reviewed)
        cursor.execute(sqlCmd, tuple(l))


class WCVPName(object):

    def __init__(self, plantNameId=None, ipniId=None, taxonRank=None, taxonStatus=None, genusHybrid=None, genus=None, speciesHybrid=None, species=None, infraspecificRank=None, infraspecies=None, parentheticalAuthor=None, primaryAuthor=None, publicationAuthor=None, placeOfPublication=None, volumeAndPage=None, firstPublished=None, nomenclaturalRemarks=None, geographicArea=None, lifeformDescription=None, climateDescription=None, taxonName=None, taxonAuthors=None, acceptedPlantName=None, basionymPlantNameId=None, replacedSynonymAuthor=None, homotypicSynonym=None, parentPlantName=None, hybridFormula=None, reviewed=None, family=None):
        self.plantNameId = plantNameId
        self.ipniId = ipniId
        self.taxonRank = taxonRank
        self.taxonStatus = taxonStatus
        self.genusHybrid = genusHybrid
        self.genus = genus
        self.speciesHybrid = speciesHybrid
        self.species = species
        self.infraspecificRank = infraspecificRank
        self.infraspecies = infraspecies
        self.parentheticalAuthor = parentheticalAuthor
        self.primaryAuthor = primaryAuthor
        self.publicationAuthor = publicationAuthor
        self.placeOfPublication = placeOfPublication
        self.volumeAndPage = volumeAndPage
        self.firstPublished = firstPublished
        self.nomenclaturalRemarks = nomenclaturalRemarks
        self.geographicArea = geographicArea
        self.lifeformDescription = lifeformDescription
        self.climateDescription = climateDescription
        self.taxonName = taxonName
        self.taxonAuthors = taxonAuthors
        self.acceptedPlantName = acceptedPlantName
        self.basionymPlantNameId = basionymPlantNameId
        self.replacedSynonymAuthor = replacedSynonymAuthor
        self.homotypicSynonym = homotypicSynonym
        self.parentPlantName = parentPlantName
        self.powoId = None
        self.hybridFormula = hybridFormula
        self.reviewed = reviewed
        self.family = family
        # one-to-many
        # fk_Specimen_Name: Specimen.POWOId REFERENCES WCVPName(POWOId)
        self.specimenPOWOList = []

    def insertIntoDatabase(self, cursor):
        sqlCmd = 'INSERT INTO `WCVPName` (`PlantNameId`, `IPNIId`, `TaxonRank`, `TaxonStatus`, `GenusHybrid`, `Genus`, `SpeciesHybrid`, `Species`, `InfraspecificRank`, `Infraspecies`, `ParentheticalAuthor`, `PrimaryAuthor`, `PublicationAuthor`, `PlaceOfPublication`, `VolumeAndPage`, `FirstPublished`, `NomenclaturalRemarks`, `GeographicArea`, `LifeformDescription`, `ClimateDescription`, `TaxonName`, `TaxonAuthors`, `AcceptedPlantName`, `BasionymPlantNameId`, `ReplacedSynonymAuthor`, `HomotypicSynonym`, `ParentPlantName`, `HybridFormula`, `Reviewed`, `idFamily`) VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s)'
        l = []
        l.append(self.plantNameId)
        l.append(self.ipniId)
        l.append(self.taxonRank)
        l.append(self.taxonStatus)
        l.append(self.genusHybrid)
        l.append(self.genus)
        l.append(self.speciesHybrid)
        l.append(self.species)
        l.append(self.infraspecificRank)
        l.append(self.infraspecies)
        l.append(self.parentheticalAuthor)
        l.append(self.primaryAuthor)
        l.append(self.publicationAuthor)
        l.append(self.placeOfPublication)
        l.append(self.volumeAndPage)
        l.append(self.firstPublished)
        l.append(self.nomenclaturalRemarks)
        l.append(self.geographicArea)
        l.append(self.lifeformDescription)
        l.append(self.climateDescription)
        l.append(self.taxonName)
        l.append(self.taxonAuthors)
        l.append(self.acceptedPlantName)
        l.append(self.basionymPlantNameId)
        l.append(self.replacedSynonymAuthor)
        l.append(self.homotypicSynonym)
        l.append(self.parentPlantName)
        l.append(self.hybridFormula)
        l.append(self.reviewed)
        l.append(None if self.family is None else self.family.idFamily)
        cursor.execute(sqlCmd, tuple(l))


class WCVPOrder(object):

    def __init__(self, order=None, apgivOrderId=None, isAngiosperm=None, subclass=None, className=None, phylum=None, kingdom=None):
        self.idOrder = None
        self.order = order
        self.apgivOrderId = apgivOrderId
        self.isAngiosperm = isAngiosperm
        self.subclass = subclass
        self.className = className
        self.phylum = phylum
        self.kingdom = kingdom
        # one-to-many
        # fk_Family_Order: WCVPFamily.idOrder REFERENCES WCVPOrder(idOrder)
        self.wcvpFamilyOrderList = []

    def insertIntoDatabase(self, cursor):
        sqlCmd = 'INSERT INTO `WCVPOrder` (`Order`, `APGIVOrderId`, `isAngiosperm`, `Subclass`, `Class`, `Phylum`, `Kingdom`) VALUES (%s, %s, %s, %s, %s, %s, %s)'
        l = []
        l.append(self.order)
        l.append(self.apgivOrderId)
        l.append(self.isAngiosperm)
        l.append(self.subclass)
        l.append(self.className)
        l.append(self.phylum)
        l.append(self.kingdom)
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


def loadAdditionalBaitKitDict(connection, productionDatabase):
    cursor = connection.cursor()
    entityDict = {}
    sqlStatement = 'SELECT `idAdditionalBaitKit`, `AdditionalBaitKit` FROM `AdditionalBaitKit`'
    cursor.execute(sqlStatement)
    for row in cursor:
        entity = AdditionalBaitKit()
        entity.idAdditionalBaitKit = paftol.database.intOrNone(row[0])
        entity.additionalBaitKit = paftol.database.strOrNone(row[1])
        entityDict[entity.idAdditionalBaitKit] = entity
    cursor.close()
    return entityDict


def loadAlternativeNameDict(connection, productionDatabase):
    cursor = connection.cursor()
    entityDict = {}
    sqlStatement = 'SELECT `idAlternativeName`, `idSpecimen`, `Order`, `Family`, `TaxonName`, `NameNote`, `DateAdded`, `idNameSource` FROM `AlternativeName`'
    cursor.execute(sqlStatement)
    for row in cursor:
        entity = AlternativeName()
        entity.idAlternativeName = paftol.database.intOrNone(row[0])
        # many to one: specimen
        entityId = paftol.database.intOrNone(row[1])
        if entityId is None:
            entity.specimen = None
        elif entityId not in productionDatabase.specimenDict:
            raise StandardError, 'no Specimen entity with idSpecimen = %d' % entityId
        else:
            entity.specimen = productionDatabase.specimenDict[entityId]
            # type: int, name: idSpecimen, foreignTable: Specimen, foreignColumn: idSpecimen
            entity.specimen.alternativeNameSpecimenList.append(entity)
        entity.order = paftol.database.strOrNone(row[2])
        entity.family = paftol.database.strOrNone(row[3])
        entity.taxonName = paftol.database.strOrNone(row[4])
        entity.nameNote = paftol.database.strOrNone(row[5])
        entity.dateAdded = paftol.database.strOrNone(row[6])
        # many to one: nameSource
        entityId = paftol.database.intOrNone(row[7])
        if entityId is None:
            entity.nameSource = None
        elif entityId not in productionDatabase.nameSourceDict:
            raise StandardError, 'no NameSource entity with idNameSource = %d' % entityId
        else:
            entity.nameSource = productionDatabase.nameSourceDict[entityId]
            # type: int, name: idNameSource, foreignTable: NameSource, foreignColumn: idNameSource
            entity.nameSource.alternativeNameNameSourceList.append(entity)
        entityDict[entity.idAlternativeName] = entity
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
    sqlStatement = 'SELECT `idDataSource`, `DataSource`, `SequenceType`, `IsAssembled` FROM `DataSource`'
    cursor.execute(sqlStatement)
    for row in cursor:
        entity = DataSource()
        entity.idDataSource = paftol.database.intOrNone(row[0])
        entity.dataSource = paftol.database.strOrNone(row[1])
        entity.sequenceType = paftol.database.strOrNone(row[2])
        entity.isAssembled = paftol.database.intOrNone(row[3])
        entityDict[entity.idDataSource] = entity
    cursor.close()
    return entityDict


def loadDataSourceAndMethodInReleaseDict(connection, productionDatabase):
    cursor = connection.cursor()
    entityDict = {}
    sqlStatement = 'SELECT `idDataSourceAndMethodInRelease`, `idDataRelease`, `idRecoveryMethod`, `idDataSource` FROM `DataSourceAndMethodInRelease`'
    cursor.execute(sqlStatement)
    for row in cursor:
        entity = DataSourceAndMethodInRelease()
        entity.idDataSourceAndMethodInRelease = paftol.database.intOrNone(row[0])
        # many to one: dataRelease
        entityId = paftol.database.intOrNone(row[1])
        if entityId is None:
            entity.dataRelease = None
        elif entityId not in productionDatabase.dataReleaseDict:
            raise StandardError, 'no DataRelease entity with idDataRelease = %d' % entityId
        else:
            entity.dataRelease = productionDatabase.dataReleaseDict[entityId]
            # type: int, name: idDataRelease, foreignTable: DataRelease, foreignColumn: idDataRelease
            entity.dataRelease.dataSourceAndMethodInReleaseDataReleaseList.append(entity)
        # many to one: recoveryMethod
        entityId = paftol.database.intOrNone(row[2])
        if entityId is None:
            entity.recoveryMethod = None
        elif entityId not in productionDatabase.recoveryMethodDict:
            raise StandardError, 'no RecoveryMethod entity with idRecoveryMethod = %d' % entityId
        else:
            entity.recoveryMethod = productionDatabase.recoveryMethodDict[entityId]
            # type: int, name: idRecoveryMethod, foreignTable: RecoveryMethod, foreignColumn: idRecoveryMethod
            entity.recoveryMethod.dataSourceAndMethodInReleaseRecoveryMethodList.append(entity)
        # many to one: dataSource
        entityId = paftol.database.intOrNone(row[3])
        if entityId is None:
            entity.dataSource = None
        elif entityId not in productionDatabase.dataSourceDict:
            raise StandardError, 'no DataSource entity with idDataSource = %d' % entityId
        else:
            entity.dataSource = productionDatabase.dataSourceDict[entityId]
            # type: int, name: idDataSource, foreignTable: DataSource, foreignColumn: idDataSource
            entity.dataSource.dataSourceAndMethodInReleaseDataSourceList.append(entity)
        entityDict[entity.idDataSourceAndMethodInRelease] = entity
    cursor.close()
    return entityDict


def loadDecisionReasonDict(connection, productionDatabase):
    cursor = connection.cursor()
    entityDict = {}
    sqlStatement = 'SELECT `idDecisionReason`, `DecisionReason` FROM `DecisionReason`'
    cursor.execute(sqlStatement)
    for row in cursor:
        entity = DecisionReason()
        entity.idDecisionReason = paftol.database.intOrNone(row[0])
        entity.decisionReason = paftol.database.strOrNone(row[1])
        entityDict[entity.idDecisionReason] = entity
    cursor.close()
    return entityDict


def loadExemplarGeneDict(connection, productionDatabase):
    cursor = connection.cursor()
    entityDict = {}
    sqlStatement = 'SELECT `idExemplarGene`, `AC`, `GN`, `DE`, `OS`, `URL` FROM `ExemplarGene`'
    cursor.execute(sqlStatement)
    for row in cursor:
        entity = ExemplarGene()
        entity.idExemplarGene = paftol.database.intOrNone(row[0])
        entity.ac = paftol.database.strOrNone(row[1])
        entity.gn = paftol.database.strOrNone(row[2])
        entity.de = paftol.database.strOrNone(row[3])
        entity.os = paftol.database.strOrNone(row[4])
        entity.url = paftol.database.strOrNone(row[5])
        entityDict[entity.idExemplarGene] = entity
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


def loadFastqStatsDict(connection, productionDatabase):
    cursor = connection.cursor()
    entityDict = {}
    sqlStatement = 'SELECT `idFastqStats`, `IsTrimmed`, `NumReads`, `Qual28`, `MeanA`, `MeanC`, `MeanG`, `MeanT`, `StddevA`, `StddevC`, `StddevG`, `StddevT`, `MeanN`, `StddevN`, `MeanAdapterContent`, `MaxAdapterContent` FROM `FastqStats`'
    cursor.execute(sqlStatement)
    for row in cursor:
        entity = FastqStats()
        entity.idFastqStats = paftol.database.intOrNone(row[0])
        entity.isTrimmed = paftol.database.intOrNone(row[1])
        entity.numReads = paftol.database.intOrNone(row[2])
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
        entityDict[entity.idFastqStats] = entity
    cursor.close()
    return entityDict


def loadGeneDict(connection, productionDatabase):
    cursor = connection.cursor()
    entityDict = {}
    sqlStatement = 'SELECT `idGene`, `GeneName`, `idGeneType`, `idExemplarGene` FROM `Gene`'
    cursor.execute(sqlStatement)
    for row in cursor:
        entity = Gene()
        entity.idGene = paftol.database.intOrNone(row[0])
        entity.geneName = paftol.database.strOrNone(row[1])
        # many to one: geneType
        entityId = paftol.database.intOrNone(row[2])
        if entityId is None:
            entity.geneType = None
        elif entityId not in productionDatabase.geneTypeDict:
            raise StandardError, 'no GeneType entity with idGeneType = %d' % entityId
        else:
            entity.geneType = productionDatabase.geneTypeDict[entityId]
            # type: int, name: idGeneType, foreignTable: GeneType, foreignColumn: idGeneType
            entity.geneType.geneGeneTypeList.append(entity)
        # many to one: exemplarGene
        entityId = paftol.database.intOrNone(row[3])
        if entityId is None:
            entity.exemplarGene = None
        elif entityId not in productionDatabase.exemplarGeneDict:
            raise StandardError, 'no ExemplarGene entity with idExemplarGene = %d' % entityId
        else:
            entity.exemplarGene = productionDatabase.exemplarGeneDict[entityId]
            # type: int, name: idExemplarGene, foreignTable: ExemplarGene, foreignColumn: idExemplarGene
            entity.exemplarGene.geneExemplarGeneList.append(entity)
        entityDict[entity.idGene] = entity
    cursor.close()
    return entityDict


def loadGeneStatsDict(connection, productionDatabase):
    cursor = connection.cursor()
    entityDict = {}
    sqlStatement = 'SELECT `idGeneStats`, `InternalName`, `ExemplarAccession`, `ExemplarName`, `ExemplarSpecies`, `ExemplarHyperlink`, `NewickFile`, `NewickFilePathName`, `AverageContigLength`, `Depth`, `AverageContigLengthPercentage`, `NumSeq`, `NumGenera`, `NumSpecies` FROM `GeneStats`'
    cursor.execute(sqlStatement)
    for row in cursor:
        entity = GeneStats()
        entity.idGeneStats = paftol.database.intOrNone(row[0])
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
        entityDict[entity.idGeneStats] = entity
    cursor.close()
    return entityDict


def loadGeneTreeDict(connection, productionDatabase):
    cursor = connection.cursor()
    entityDict = {}
    sqlStatement = 'SELECT `idGeneTree`, `idGene`, `UnAlnFastaFilePathName`, `AlnFastaFilePathName`, `NewickTree`, `NewickFilePathName`, `idSpeciesTree` FROM `GeneTree`'
    cursor.execute(sqlStatement)
    for row in cursor:
        entity = GeneTree()
        entity.idGeneTree = paftol.database.intOrNone(row[0])
        # many to one: gene
        entityId = paftol.database.intOrNone(row[1])
        if entityId is None:
            entity.gene = None
        elif entityId not in productionDatabase.geneDict:
            raise StandardError, 'no Gene entity with idGene = %d' % entityId
        else:
            entity.gene = productionDatabase.geneDict[entityId]
            # type: int, name: idGene, foreignTable: Gene, foreignColumn: idGene
            entity.gene.geneTreeGeneList.append(entity)
        entity.unAlnFastaFilePathName = paftol.database.strOrNone(row[2])
        entity.alnFastaFilePathName = paftol.database.strOrNone(row[3])
        entity.newickTree = paftol.database.strOrNone(row[4])
        entity.newickFilePathName = paftol.database.strOrNone(row[5])
        # many to one: speciesTree
        entityId = paftol.database.intOrNone(row[6])
        if entityId is None:
            entity.speciesTree = None
        elif entityId not in productionDatabase.speciesTreeDict:
            raise StandardError, 'no SpeciesTree entity with idSpeciesTree = %d' % entityId
        else:
            entity.speciesTree = productionDatabase.speciesTreeDict[entityId]
            # type: int, name: idSpeciesTree, foreignTable: SpeciesTree, foreignColumn: idSpeciesTree
            entity.speciesTree.geneTreeSpeciesTreeList.append(entity)
        entityDict[entity.idGeneTree] = entity
    cursor.close()
    return entityDict


def loadGeneTypeDict(connection, productionDatabase):
    cursor = connection.cursor()
    entityDict = {}
    sqlStatement = 'SELECT `idGeneType`, `GeneTypeName` FROM `GeneType`'
    cursor.execute(sqlStatement)
    for row in cursor:
        entity = GeneType()
        entity.idGeneType = paftol.database.intOrNone(row[0])
        entity.geneTypeName = paftol.database.strOrNone(row[1])
        entityDict[entity.idGeneType] = entity
    cursor.close()
    return entityDict


def loadHerbariumDict(connection, productionDatabase):
    cursor = connection.cursor()
    entityDict = {}
    sqlStatement = 'SELECT `idHerbarium`, `HerbariumCode`, `HerbariumName`, `HerbariumURL` FROM `Herbarium`'
    cursor.execute(sqlStatement)
    for row in cursor:
        entity = Herbarium()
        entity.idHerbarium = paftol.database.intOrNone(row[0])
        entity.herbariumCode = paftol.database.strOrNone(row[1])
        entity.herbariumName = paftol.database.strOrNone(row[2])
        entity.herbariumUrL = paftol.database.strOrNone(row[3])
        entityDict[entity.idHerbarium] = entity
    cursor.close()
    return entityDict


def loadISOCountryDict(connection, productionDatabase):
    cursor = connection.cursor()
    entityDict = {}
    sqlStatement = 'SELECT `idISOCountry`, `Country`, `Alpha2Code`, `Alpha3Code`, `CountryCode`, `ISOCode`, `Region`, `SubRegion`, `IntermediateRegion`, `RegionCode`, `SubRegionCode`, `IntermediateRegionCode` FROM `ISOCountry`'
    cursor.execute(sqlStatement)
    for row in cursor:
        entity = ISOCountry()
        entity.idIsoCountry = paftol.database.intOrNone(row[0])
        entity.country = paftol.database.strOrNone(row[1])
        entity.alpha2Code = paftol.database.strOrNone(row[2])
        entity.alpha3Code = paftol.database.strOrNone(row[3])
        entity.countryCode = paftol.database.intOrNone(row[4])
        entity.isoCode = paftol.database.strOrNone(row[5])
        entity.region = paftol.database.strOrNone(row[6])
        entity.subRegion = paftol.database.strOrNone(row[7])
        entity.intermediateRegion = paftol.database.strOrNone(row[8])
        entity.regionCode = paftol.database.intOrNone(row[9])
        entity.subRegionCode = paftol.database.intOrNone(row[10])
        entity.intermediateRegionCode = paftol.database.intOrNone(row[11])
        entityDict[entity.idIsoCountry] = entity
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
    sqlStatement = 'SELECT `idLibrary`, `idSample`, `LibConcentration`, `LibQuality`, `LibTapeStation`, `Sonication`, `PlateNumber`, `idCoordinate`, `Description`, `idStatus`, `idIndexes`, `GenerateLibrary`, `PCRCycles` FROM `Library`'
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
            # type: int, name: idSample, foreignTable: Sample, foreignColumn: idSample
            entity.sample.librarySampleList.append(entity)
        entity.libConcentration = paftol.database.floatOrNone(row[2])
        entity.libQuality = paftol.database.strOrNone(row[3])
        entity.libTapeStation = paftol.database.strOrNone(row[4])
        entity.sonication = paftol.database.intOrNone(row[5])
        entity.plateNumber = paftol.database.strOrNone(row[6])
        # many to one: coordinate
        entityId = paftol.database.intOrNone(row[7])
        if entityId is None:
            entity.coordinate = None
        elif entityId not in productionDatabase.coordinatesDict:
            raise StandardError, 'no Coordinates entity with idCoordinate = %d' % entityId
        else:
            entity.coordinate = productionDatabase.coordinatesDict[entityId]
            # type: int, name: idCoordinate, foreignTable: Coordinates, foreignColumn: idCoordinate
            entity.coordinate.libraryCoordinateList.append(entity)
        entity.description = paftol.database.strOrNone(row[8])
        # many to one: status
        entityId = paftol.database.intOrNone(row[9])
        if entityId is None:
            entity.status = None
        elif entityId not in productionDatabase.statusDict:
            raise StandardError, 'no Status entity with idStatus = %d' % entityId
        else:
            entity.status = productionDatabase.statusDict[entityId]
            # type: int, name: idStatus, foreignTable: Status, foreignColumn: idStatus
            entity.status.libraryStatusList.append(entity)
        # many to one: indexes
        entityId = paftol.database.intOrNone(row[10])
        if entityId is None:
            entity.indexes = None
        elif entityId not in productionDatabase.indexesDict:
            raise StandardError, 'no Indexes entity with idIndexes = %d' % entityId
        else:
            entity.indexes = productionDatabase.indexesDict[entityId]
            # type: int, name: idIndexes, foreignTable: Indexes, foreignColumn: idIndexes
            entity.indexes.libraryIndexesList.append(entity)
        entity.generateLibrary = paftol.database.intOrNone(row[11])
        entity.pcrCycles = paftol.database.intOrNone(row[12])
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


def loadMergedSequenceDict(connection, productionDatabase):
    cursor = connection.cursor()
    entityDict = {}
    sqlStatement = 'SELECT `idMergedSequence`, `idSequence`, `idSequenceOrigin` FROM `MergedSequence`'
    cursor.execute(sqlStatement)
    for row in cursor:
        entity = MergedSequence()
        entity.idMergedSequence = paftol.database.intOrNone(row[0])
        # many to one: sequence
        entityId = paftol.database.intOrNone(row[1])
        if entityId is None:
            entity.sequence = None
        elif entityId not in productionDatabase.sequenceDict:
            raise StandardError, 'no Sequence entity with idSequence = %d' % entityId
        else:
            entity.sequence = productionDatabase.sequenceDict[entityId]
            # type: int, name: idSequence, foreignTable: Sequence, foreignColumn: idSequence
            entity.sequence.mergedSequenceSequenceList.append(entity)
        # many to one: sequenceOrigin
        entityId = paftol.database.intOrNone(row[2])
        if entityId is None:
            entity.sequenceOrigin = None
        elif entityId not in productionDatabase.sequenceDict:
            raise StandardError, 'no Sequence entity with idSequence = %d' % entityId
        else:
            entity.sequenceOrigin = productionDatabase.sequenceDict[entityId]
            # type: int, name: idSequenceOrigin, foreignTable: Sequence, foreignColumn: idSequence
            entity.sequenceOrigin.mergedSequenceSequenceOriginList.append(entity)
        entityDict[entity.idMergedSequence] = entity
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


def loadNameSourceDict(connection, productionDatabase):
    cursor = connection.cursor()
    entityDict = {}
    sqlStatement = 'SELECT `idNameSource`, `NameSource` FROM `NameSource`'
    cursor.execute(sqlStatement)
    for row in cursor:
        entity = NameSource()
        entity.idNameSource = paftol.database.intOrNone(row[0])
        entity.nameSource = paftol.database.strOrNone(row[1])
        entityDict[entity.idNameSource] = entity
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
    sqlStatement = 'SELECT `idProject`, `Project` FROM `Project`'
    cursor.execute(sqlStatement)
    for row in cursor:
        entity = Project()
        entity.idProject = paftol.database.intOrNone(row[0])
        entity.project = paftol.database.strOrNone(row[1])
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


def loadRawFastaFileDict(connection, productionDatabase):
    cursor = connection.cursor()
    entityDict = {}
    sqlStatement = 'SELECT `idRawFastaFile`, `FilePathName`, `UncompressedFileMd5sum`, `NumRecords`, `idSequence` FROM `RawFastaFile`'
    cursor.execute(sqlStatement)
    for row in cursor:
        entity = RawFastaFile()
        entity.idRawFastaFile = paftol.database.intOrNone(row[0])
        entity.filePathName = paftol.database.strOrNone(row[1])
        entity.uncompressedFileMd5sum = paftol.database.strOrNone(row[2])
        entity.numRecords = paftol.database.intOrNone(row[3])
        # many to one: sequence
        entityId = paftol.database.intOrNone(row[4])
        if entityId is None:
            entity.sequence = None
        elif entityId not in productionDatabase.sequenceDict:
            raise StandardError, 'no Sequence entity with idSequence = %d' % entityId
        else:
            entity.sequence = productionDatabase.sequenceDict[entityId]
            # type: int, name: idSequence, foreignTable: Sequence, foreignColumn: idSequence
            entity.sequence.rawFastaFileSequenceList.append(entity)
        entityDict[entity.idRawFastaFile] = entity
    cursor.close()
    return entityDict


def loadRawFastqFileDict(connection, productionDatabase):
    cursor = connection.cursor()
    entityDict = {}
    sqlStatement = 'SELECT `idRawFastqFile`, `FilePathName`, `UncompressedFileMd5sum`, `Orientation`, `idSequence`, `idFastqStats` FROM `RawFastqFile`'
    cursor.execute(sqlStatement)
    for row in cursor:
        entity = RawFastqFile()
        entity.idRawFastqFile = paftol.database.intOrNone(row[0])
        entity.filePathName = paftol.database.strOrNone(row[1])
        entity.uncompressedFileMd5sum = paftol.database.strOrNone(row[2])
        entity.orientation = paftol.database.strOrNone(row[3])
        # many to one: sequence
        entityId = paftol.database.intOrNone(row[4])
        if entityId is None:
            entity.sequence = None
        elif entityId not in productionDatabase.sequenceDict:
            raise StandardError, 'no Sequence entity with idSequence = %d' % entityId
        else:
            entity.sequence = productionDatabase.sequenceDict[entityId]
            # type: int, name: idSequence, foreignTable: Sequence, foreignColumn: idSequence
            entity.sequence.rawFastqFileSequenceList.append(entity)
        # many to one: fastqStats
        entityId = paftol.database.intOrNone(row[5])
        if entityId is None:
            entity.fastqStats = None
        elif entityId not in productionDatabase.fastqStatsDict:
            raise StandardError, 'no FastqStats entity with idFastqStats = %d' % entityId
        else:
            entity.fastqStats = productionDatabase.fastqStatsDict[entityId]
            # type: int, name: idFastqStats, foreignTable: FastqStats, foreignColumn: idFastqStats
            entity.fastqStats.rawFastqFileFastqStatsList.append(entity)
        entityDict[entity.idRawFastqFile] = entity
    cursor.close()
    return entityDict


def loadRecoveredContigDict(connection, productionDatabase):
    cursor = connection.cursor()
    entityDict = {}
    sqlStatement = 'SELECT `idRecoveredContig`, `idSequenceRecovery`, `Seqlength`, `ContigMd5sum`, `idRepresentativeReferenceTarget` FROM `RecoveredContig`'
    cursor.execute(sqlStatement)
    for row in cursor:
        entity = RecoveredContig()
        entity.idRecoveredContig = paftol.database.intOrNone(row[0])
        # many to one: sequenceRecovery
        entityId = paftol.database.intOrNone(row[1])
        if entityId is None:
            entity.sequenceRecovery = None
        elif entityId not in productionDatabase.sequenceRecoveryDict:
            raise StandardError, 'no SequenceRecovery entity with idSequenceRecovery = %d' % entityId
        else:
            entity.sequenceRecovery = productionDatabase.sequenceRecoveryDict[entityId]
            # type: int, name: idSequenceRecovery, foreignTable: SequenceRecovery, foreignColumn: idSequenceRecovery
            entity.sequenceRecovery.recoveredContigSequenceRecoveryList.append(entity)
        entity.seqlength = paftol.database.intOrNone(row[2])
        entity.contigMd5sum = paftol.database.strOrNone(row[3])
        # many to one: representativeReferenceTarget
        entityId = paftol.database.intOrNone(row[4])
        if entityId is None:
            entity.representativeReferenceTarget = None
        elif entityId not in productionDatabase.referenceTargetDict:
            raise StandardError, 'no ReferenceTarget entity with idReferenceTarget = %d' % entityId
        else:
            entity.representativeReferenceTarget = productionDatabase.referenceTargetDict[entityId]
            # type: int, name: idRepresentativeReferenceTarget, foreignTable: ReferenceTarget, foreignColumn: idReferenceTarget
            entity.representativeReferenceTarget.recoveredContigRepresentativeReferenceTargetList.append(entity)
        entityDict[entity.idRecoveredContig] = entity
    cursor.close()
    return entityDict


def loadRecoveryMethodDict(connection, productionDatabase):
    cursor = connection.cursor()
    entityDict = {}
    sqlStatement = 'SELECT `idRecoveryMethod`, `Software`, `SoftwareVersion`, `SoftwareParameters`, `idTargetSet` FROM `RecoveryMethod`'
    cursor.execute(sqlStatement)
    for row in cursor:
        entity = RecoveryMethod()
        entity.idRecoveryMethod = paftol.database.intOrNone(row[0])
        entity.software = paftol.database.strOrNone(row[1])
        entity.softwareVersion = paftol.database.strOrNone(row[2])
        entity.softwareParameters = paftol.database.strOrNone(row[3])
        # many to one: targetSet
        entityId = paftol.database.intOrNone(row[4])
        if entityId is None:
            entity.targetSet = None
        elif entityId not in productionDatabase.targetSetDict:
            raise StandardError, 'no TargetSet entity with idTargetSet = %d' % entityId
        else:
            entity.targetSet = productionDatabase.targetSetDict[entityId]
            # type: int, name: idTargetSet, foreignTable: TargetSet, foreignColumn: idTargetSet
            entity.targetSet.recoveryMethodTargetSetList.append(entity)
        entityDict[entity.idRecoveryMethod] = entity
    cursor.close()
    return entityDict


def loadReferenceTargetDict(connection, productionDatabase):
    cursor = connection.cursor()
    entityDict = {}
    sqlStatement = 'SELECT `idReferenceTarget`, `idGene`, `Organism`, `TargetLength`, `idTargetSet` FROM `ReferenceTarget`'
    cursor.execute(sqlStatement)
    for row in cursor:
        entity = ReferenceTarget()
        entity.idReferenceTarget = paftol.database.intOrNone(row[0])
        # many to one: gene
        entityId = paftol.database.intOrNone(row[1])
        if entityId is None:
            entity.gene = None
        elif entityId not in productionDatabase.geneDict:
            raise StandardError, 'no Gene entity with idGene = %d' % entityId
        else:
            entity.gene = productionDatabase.geneDict[entityId]
            # type: int, name: idGene, foreignTable: Gene, foreignColumn: idGene
            entity.gene.referenceTargetGeneList.append(entity)
        entity.organism = paftol.database.strOrNone(row[2])
        entity.targetLength = paftol.database.intOrNone(row[3])
        # many to one: targetSet
        entityId = paftol.database.intOrNone(row[4])
        if entityId is None:
            entity.targetSet = None
        elif entityId not in productionDatabase.targetSetDict:
            raise StandardError, 'no TargetSet entity with idTargetSet = %d' % entityId
        else:
            entity.targetSet = productionDatabase.targetSetDict[entityId]
            # type: int, name: idTargetSet, foreignTable: TargetSet, foreignColumn: idTargetSet
            entity.targetSet.referenceTargetTargetSetList.append(entity)
        entityDict[entity.idReferenceTarget] = entity
    cursor.close()
    return entityDict


def loadSampleDict(connection, productionDatabase):
    cursor = connection.cursor()
    entityDict = {}
    sqlStatement = 'SELECT `idSample`, `Description`, `idAction`, `idExtractionType`, `idQuality`, `SampleConcentration`, `GelImage`, `SampleTapeStation`, `idDNAVolume`, `ENASampleNum`, `SecENASampleNum`, `idSpecimen` FROM `Sample`'
    cursor.execute(sqlStatement)
    for row in cursor:
        entity = Sample()
        entity.idSample = paftol.database.intOrNone(row[0])
        entity.description = paftol.database.strOrNone(row[1])
        # many to one: action
        entityId = paftol.database.intOrNone(row[2])
        if entityId is None:
            entity.action = None
        elif entityId not in productionDatabase.actionDict:
            raise StandardError, 'no Action entity with idAction = %d' % entityId
        else:
            entity.action = productionDatabase.actionDict[entityId]
            # type: int, name: idAction, foreignTable: Action, foreignColumn: idAction
            entity.action.sampleActionList.append(entity)
        # many to one: extractionType
        entityId = paftol.database.intOrNone(row[3])
        if entityId is None:
            entity.extractionType = None
        elif entityId not in productionDatabase.extractionTypeDict:
            raise StandardError, 'no ExtractionType entity with idExtractionType = %d' % entityId
        else:
            entity.extractionType = productionDatabase.extractionTypeDict[entityId]
            # type: int, name: idExtractionType, foreignTable: ExtractionType, foreignColumn: idExtractionType
            entity.extractionType.sampleExtractionTypeList.append(entity)
        # many to one: quality
        entityId = paftol.database.intOrNone(row[4])
        if entityId is None:
            entity.quality = None
        elif entityId not in productionDatabase.qualityDict:
            raise StandardError, 'no Quality entity with idQuality = %d' % entityId
        else:
            entity.quality = productionDatabase.qualityDict[entityId]
            # type: int, name: idQuality, foreignTable: Quality, foreignColumn: idQuality
            entity.quality.sampleQualityList.append(entity)
        entity.sampleConcentration = paftol.database.floatOrNone(row[5])
        entity.gelImage = paftol.database.strOrNone(row[6])
        entity.sampleTapeStation = paftol.database.strOrNone(row[7])
        # many to one: dnaVolume
        entityId = paftol.database.intOrNone(row[8])
        if entityId is None:
            entity.dnaVolume = None
        elif entityId not in productionDatabase.dNAVolumeDict:
            raise StandardError, 'no DNAVolume entity with idDnaVolume = %d' % entityId
        else:
            entity.dnaVolume = productionDatabase.dNAVolumeDict[entityId]
            # type: int, name: idDNAVolume, foreignTable: DNAVolume, foreignColumn: idDNAVolume
            entity.dnaVolume.sampleDnaVolumeList.append(entity)
        entity.enaSampleNum = paftol.database.strOrNone(row[9])
        entity.secEnaSampleNum = paftol.database.strOrNone(row[10])
        # many to one: specimen
        entityId = paftol.database.intOrNone(row[11])
        if entityId is None:
            entity.specimen = None
        elif entityId not in productionDatabase.specimenDict:
            raise StandardError, 'no Specimen entity with idSpecimen = %d' % entityId
        else:
            entity.specimen = productionDatabase.specimenDict[entityId]
            # type: int, name: idSpecimen, foreignTable: Specimen, foreignColumn: idSpecimen
            entity.specimen.sampleSpecimenList.append(entity)
        entityDict[entity.idSample] = entity
    cursor.close()
    return entityDict


def loadSequenceDict(connection, productionDatabase):
    cursor = connection.cursor()
    entityDict = {}
    sqlStatement = 'SELECT `idSequence`, `ExternalSequenceID`, `HasDuplicate`, `IsMerged`, `idLibrary`, `idPlatform`, `idLocation`, `SequencingRun`, `NumInferredCds`, `MedianHybpiperCdsLength`, `idStatus`, `HybridisationPool`, `R2FastqFile`, `R1FastqFile`, `Blacklisted`, `idBlacklistedReason`, `BlacklistedPerson`, `idSequencingStrategy`, `ENAExpNumber`, `ENARunNumber`, `SuspiciousPlacement`, `idAdditionalBaitKit` FROM `Sequence`'
    cursor.execute(sqlStatement)
    for row in cursor:
        entity = Sequence()
        entity.idSequence = paftol.database.intOrNone(row[0])
        entity.externalSequenceId = paftol.database.strOrNone(row[1])
        entity.hasDuplicate = paftol.database.intOrNone(row[2])
        entity.isMerged = paftol.database.intOrNone(row[3])
        # many to one: library
        entityId = paftol.database.intOrNone(row[4])
        if entityId is None:
            entity.library = None
        elif entityId not in productionDatabase.libraryDict:
            raise StandardError, 'no Library entity with idLibrary = %d' % entityId
        else:
            entity.library = productionDatabase.libraryDict[entityId]
            # type: int, name: idLibrary, foreignTable: Library, foreignColumn: idLibrary
            entity.library.sequenceLibraryList.append(entity)
        # many to one: platform
        entityId = paftol.database.intOrNone(row[5])
        if entityId is None:
            entity.platform = None
        elif entityId not in productionDatabase.platformDict:
            raise StandardError, 'no Platform entity with idPlatform = %d' % entityId
        else:
            entity.platform = productionDatabase.platformDict[entityId]
            # type: int, name: idPlatform, foreignTable: Platform, foreignColumn: idPlatform
            entity.platform.sequencePlatformList.append(entity)
        # many to one: location
        entityId = paftol.database.intOrNone(row[6])
        if entityId is None:
            entity.location = None
        elif entityId not in productionDatabase.locationDict:
            raise StandardError, 'no Location entity with idLocation = %d' % entityId
        else:
            entity.location = productionDatabase.locationDict[entityId]
            # type: int, name: idLocation, foreignTable: Location, foreignColumn: idLocation
            entity.location.sequenceLocationList.append(entity)
        entity.sequencingRun = paftol.database.strOrNone(row[7])
        entity.numInferredCds = paftol.database.intOrNone(row[8])
        entity.medianHybpiperCdsLength = paftol.database.floatOrNone(row[9])
        # many to one: status
        entityId = paftol.database.intOrNone(row[10])
        if entityId is None:
            entity.status = None
        elif entityId not in productionDatabase.statusDict:
            raise StandardError, 'no Status entity with idStatus = %d' % entityId
        else:
            entity.status = productionDatabase.statusDict[entityId]
            # type: int, name: idStatus, foreignTable: Status, foreignColumn: idStatus
            entity.status.sequenceStatusList.append(entity)
        entity.hybridisationPool = paftol.database.strOrNone(row[11])
        entity.r2FastqFile = paftol.database.strOrNone(row[12])
        entity.r1FastqFile = paftol.database.strOrNone(row[13])
        entity.blacklisted = paftol.database.intOrNone(row[14])
        # many to one: blacklistedReason
        entityId = paftol.database.intOrNone(row[15])
        if entityId is None:
            entity.blacklistedReason = None
        elif entityId not in productionDatabase.blacklistedReasonDict:
            raise StandardError, 'no BlacklistedReason entity with idBlacklistedReason = %d' % entityId
        else:
            entity.blacklistedReason = productionDatabase.blacklistedReasonDict[entityId]
            # type: int, name: idBlacklistedReason, foreignTable: BlacklistedReason, foreignColumn: idBlacklistedReason
            entity.blacklistedReason.sequenceBlacklistedReasonList.append(entity)
        entity.blacklistedPerson = paftol.database.strOrNone(row[16])
        # many to one: sequencingStrategy
        entityId = paftol.database.intOrNone(row[17])
        if entityId is None:
            entity.sequencingStrategy = None
        elif entityId not in productionDatabase.sequencingStrategyDict:
            raise StandardError, 'no SequencingStrategy entity with idSequencingStrategy = %d' % entityId
        else:
            entity.sequencingStrategy = productionDatabase.sequencingStrategyDict[entityId]
            # type: int, name: idSequencingStrategy, foreignTable: SequencingStrategy, foreignColumn: idSequencingStrategy
            entity.sequencingStrategy.sequenceSequencingStrategyList.append(entity)
        entity.enaExpNumber = paftol.database.strOrNone(row[18])
        entity.enaRunNumber = paftol.database.strOrNone(row[19])
        entity.suspiciousPlacement = paftol.database.intOrNone(row[20])
        # many to one: additionalBaitKit
        entityId = paftol.database.intOrNone(row[21])
        if entityId is None:
            entity.additionalBaitKit = None
        elif entityId not in productionDatabase.additionalBaitKitDict:
            raise StandardError, 'no AdditionalBaitKit entity with idAdditionalBaitKit = %d' % entityId
        else:
            entity.additionalBaitKit = productionDatabase.additionalBaitKitDict[entityId]
            # type: int, name: idAdditionalBaitKit, foreignTable: AdditionalBaitKit, foreignColumn: idAdditionalBaitKit
            entity.additionalBaitKit.sequenceAdditionalBaitKitList.append(entity)
        entityDict[entity.idSequence] = entity
    cursor.close()
    return entityDict


def loadSequenceDataReleaseDict(connection, productionDatabase):
    cursor = connection.cursor()
    entityDict = {}
    sqlStatement = 'SELECT `idSequenceDataRelease`, `idSequence`, `idDataRelease`, `idBarcodeValidation`, `idPhylogeneticValidation`, `idValidationResult`, `idDecisionReason`, `ValidationComments` FROM `SequenceDataRelease`'
    cursor.execute(sqlStatement)
    for row in cursor:
        entity = SequenceDataRelease()
        entity.idSequenceDataRelease = paftol.database.intOrNone(row[0])
        # many to one: sequence
        entityId = paftol.database.intOrNone(row[1])
        if entityId is None:
            entity.sequence = None
        elif entityId not in productionDatabase.sequenceDict:
            raise StandardError, 'no Sequence entity with idSequence = %d' % entityId
        else:
            entity.sequence = productionDatabase.sequenceDict[entityId]
            # type: int, name: idSequence, foreignTable: Sequence, foreignColumn: idSequence
            entity.sequence.sequenceDataReleaseSequenceList.append(entity)
        # many to one: dataRelease
        entityId = paftol.database.intOrNone(row[2])
        if entityId is None:
            entity.dataRelease = None
        elif entityId not in productionDatabase.dataReleaseDict:
            raise StandardError, 'no DataRelease entity with idDataRelease = %d' % entityId
        else:
            entity.dataRelease = productionDatabase.dataReleaseDict[entityId]
            # type: int, name: idDataRelease, foreignTable: DataRelease, foreignColumn: idDataRelease
            entity.dataRelease.sequenceDataReleaseDataReleaseList.append(entity)
        # many to one: barcodeValidation
        entityId = paftol.database.intOrNone(row[3])
        if entityId is None:
            entity.barcodeValidation = None
        elif entityId not in productionDatabase.testResultDict:
            raise StandardError, 'no TestResult entity with idTestResult = %d' % entityId
        else:
            entity.barcodeValidation = productionDatabase.testResultDict[entityId]
            # type: int, name: idBarcodeValidation, foreignTable: TestResult, foreignColumn: idTestResult
            entity.barcodeValidation.sequenceDataReleaseBarcodeValidationList.append(entity)
        # many to one: phylogeneticValidation
        entityId = paftol.database.intOrNone(row[4])
        if entityId is None:
            entity.phylogeneticValidation = None
        elif entityId not in productionDatabase.testResultDict:
            raise StandardError, 'no TestResult entity with idTestResult = %d' % entityId
        else:
            entity.phylogeneticValidation = productionDatabase.testResultDict[entityId]
            # type: int, name: idPhylogeneticValidation, foreignTable: TestResult, foreignColumn: idTestResult
            entity.phylogeneticValidation.sequenceDataReleasePhylogeneticValidationList.append(entity)
        # many to one: validationResult
        entityId = paftol.database.intOrNone(row[5])
        if entityId is None:
            entity.validationResult = None
        elif entityId not in productionDatabase.validationResultDict:
            raise StandardError, 'no ValidationResult entity with idValidationResult = %d' % entityId
        else:
            entity.validationResult = productionDatabase.validationResultDict[entityId]
            # type: int, name: idValidationResult, foreignTable: ValidationResult, foreignColumn: idValidationResult
            entity.validationResult.sequenceDataReleaseValidationResultList.append(entity)
        # many to one: decisionReason
        entityId = paftol.database.intOrNone(row[6])
        if entityId is None:
            entity.decisionReason = None
        elif entityId not in productionDatabase.decisionReasonDict:
            raise StandardError, 'no DecisionReason entity with idDecisionReason = %d' % entityId
        else:
            entity.decisionReason = productionDatabase.decisionReasonDict[entityId]
            # type: int, name: idDecisionReason, foreignTable: DecisionReason, foreignColumn: idDecisionReason
            entity.decisionReason.sequenceDataReleaseDecisionReasonList.append(entity)
        entity.validationComments = paftol.database.strOrNone(row[7])
        entityDict[entity.idSequenceDataRelease] = entity
    cursor.close()
    return entityDict


def loadSequenceGeneStatsDict(connection, productionDatabase):
    cursor = connection.cursor()
    entityDict = {}
    sqlStatement = 'SELECT `idSequenceGeneStats`, `idSequence`, `NumRecoveredGenes`, `SumContigLength` FROM `SequenceGeneStats`'
    cursor.execute(sqlStatement)
    for row in cursor:
        entity = SequenceGeneStats()
        entity.idSequenceGeneStats = paftol.database.intOrNone(row[0])
        # many to one: sequence
        entityId = paftol.database.intOrNone(row[1])
        if entityId is None:
            entity.sequence = None
        elif entityId not in productionDatabase.sequenceDict:
            raise StandardError, 'no Sequence entity with idSequence = %d' % entityId
        else:
            entity.sequence = productionDatabase.sequenceDict[entityId]
            # type: int, name: idSequence, foreignTable: Sequence, foreignColumn: idSequence
            entity.sequence.sequenceGeneStatsSequenceList.append(entity)
        entity.numRecoveredGenes = paftol.database.intOrNone(row[2])
        entity.sumContigLength = paftol.database.floatOrNone(row[3])
        entityDict[entity.idSequenceGeneStats] = entity
    cursor.close()
    return entityDict


def loadSequenceRawReadsDict(connection, productionDatabase):
    cursor = connection.cursor()
    entityDict = {}
    sqlStatement = 'SELECT `idSequenceRawReads`, `idSequence`, `NumReads` FROM `SequenceRawReads`'
    cursor.execute(sqlStatement)
    for row in cursor:
        entity = SequenceRawReads()
        entity.idSequenceRawReads = paftol.database.intOrNone(row[0])
        # many to one: sequence
        entityId = paftol.database.intOrNone(row[1])
        if entityId is None:
            entity.sequence = None
        elif entityId not in productionDatabase.sequenceDict:
            raise StandardError, 'no Sequence entity with idSequence = %d' % entityId
        else:
            entity.sequence = productionDatabase.sequenceDict[entityId]
            # type: int, name: idSequence, foreignTable: Sequence, foreignColumn: idSequence
            entity.sequence.sequenceRawReadsSequenceList.append(entity)
        entity.numReads = paftol.database.intOrNone(row[2])
        entityDict[entity.idSequenceRawReads] = entity
    cursor.close()
    return entityDict


def loadSequenceRecoveryDict(connection, productionDatabase):
    cursor = connection.cursor()
    entityDict = {}
    sqlStatement = 'SELECT `idSequenceRecovery`, `idSequence`, `ContigFastaFilePathName`, `ContigFastaFileMd5sum`, `NumMappedReads`, `NumUnmappedReads`, `SoftwareVersion`, `CmdLine`, `NumRecoveredContigsCheck`, `idRecoveryMethod` FROM `SequenceRecovery`'
    cursor.execute(sqlStatement)
    for row in cursor:
        entity = SequenceRecovery()
        entity.idSequenceRecovery = paftol.database.intOrNone(row[0])
        # many to one: sequence
        entityId = paftol.database.intOrNone(row[1])
        if entityId is None:
            entity.sequence = None
        elif entityId not in productionDatabase.sequenceDict:
            raise StandardError, 'no Sequence entity with idSequence = %d' % entityId
        else:
            entity.sequence = productionDatabase.sequenceDict[entityId]
            # type: int, name: idSequence, foreignTable: Sequence, foreignColumn: idSequence
            entity.sequence.sequenceRecoverySequenceList.append(entity)
        entity.contigFastaFilePathName = paftol.database.strOrNone(row[2])
        entity.contigFastaFileMd5sum = paftol.database.strOrNone(row[3])
        entity.numMappedReads = paftol.database.intOrNone(row[4])
        entity.numUnmappedReads = paftol.database.intOrNone(row[5])
        entity.softwareVersion = paftol.database.strOrNone(row[6])
        entity.cmdLine = paftol.database.strOrNone(row[7])
        entity.numRecoveredContigsCheck = paftol.database.intOrNone(row[8])
        # many to one: recoveryMethod
        entityId = paftol.database.intOrNone(row[9])
        if entityId is None:
            entity.recoveryMethod = None
        elif entityId not in productionDatabase.recoveryMethodDict:
            raise StandardError, 'no RecoveryMethod entity with idRecoveryMethod = %d' % entityId
        else:
            entity.recoveryMethod = productionDatabase.recoveryMethodDict[entityId]
            # type: int, name: idRecoveryMethod, foreignTable: RecoveryMethod, foreignColumn: idRecoveryMethod
            entity.recoveryMethod.sequenceRecoveryRecoveryMethodList.append(entity)
        entityDict[entity.idSequenceRecovery] = entity
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


def loadSpeciesTreeDict(connection, productionDatabase):
    cursor = connection.cursor()
    entityDict = {}
    sqlStatement = 'SELECT `idSpeciesTree`, `NewickTree`, `NewickFilePathName`, `PipelineVersion`, `PipelineCmdLine`, `idDataRelease` FROM `SpeciesTree`'
    cursor.execute(sqlStatement)
    for row in cursor:
        entity = SpeciesTree()
        entity.idSpeciesTree = paftol.database.intOrNone(row[0])
        entity.newickTree = paftol.database.strOrNone(row[1])
        entity.newickFilePathName = paftol.database.strOrNone(row[2])
        entity.pipelineVersion = paftol.database.strOrNone(row[3])
        entity.pipelineCmdLine = paftol.database.strOrNone(row[4])
        # many to one: dataRelease
        entityId = paftol.database.intOrNone(row[5])
        if entityId is None:
            entity.dataRelease = None
        elif entityId not in productionDatabase.dataReleaseDict:
            raise StandardError, 'no DataRelease entity with idDataRelease = %d' % entityId
        else:
            entity.dataRelease = productionDatabase.dataReleaseDict[entityId]
            # type: int, name: idDataRelease, foreignTable: DataRelease, foreignColumn: idDataRelease
            entity.dataRelease.speciesTreeDataReleaseList.append(entity)
        entityDict[entity.idSpeciesTree] = entity
    cursor.close()
    return entityDict


def loadSpecimenDict(connection, productionDatabase):
    cursor = connection.cursor()
    entityDict = {}
    sqlStatement = 'SELECT `idSpecimen`, `POWOId`, `BankID`, `LCD`, `MSB`, `Collector`, `CollectorNo`, `VoucherNo`, `MuseumBarcode`, `OldSpeciesName`, `idSourceSpecimen`, `idProject`, `idISOCountry`, `idMaterialSource`, `AgeOfMaterial`, `idHerbarium`, `SpecimenReference`, `HerbcatURL`, `Blacklisted`, `HoldDate`, `idDataSource` FROM `Specimen`'
    cursor.execute(sqlStatement)
    for row in cursor:
        entity = Specimen()
        entity.idSpecimen = paftol.database.intOrNone(row[0])
        # many to one: POWO
        entityId = paftol.database.strOrNone(row[1])
        if entityId is None:
            entity.POWO = None
        elif entityId not in productionDatabase.wCVPNameDict:
            raise StandardError, 'no WCVPName entity with powoId = %d' % entityId
        else:
            entity.POWO = productionDatabase.wCVPNameDict[entityId]
            # type: varchar, name: POWOId, foreignTable: WCVPName, foreignColumn: POWOId
            entity.POWO.specimenPOWOList.append(entity)
        entity.bankId = paftol.database.intOrNone(row[2])
        entity.lcd = paftol.database.strOrNone(row[3])
        entity.msb = paftol.database.intOrNone(row[4])
        entity.collector = paftol.database.strOrNone(row[5])
        entity.collectorNo = paftol.database.strOrNone(row[6])
        entity.voucherNo = paftol.database.strOrNone(row[7])
        entity.museumBarcode = paftol.database.strOrNone(row[8])
        entity.oldSpeciesName = paftol.database.strOrNone(row[9])
        # many to one: sourceSpecimen
        entityId = paftol.database.intOrNone(row[10])
        if entityId is None:
            entity.sourceSpecimen = None
        elif entityId not in productionDatabase.sourceSpecimenDict:
            raise StandardError, 'no SourceSpecimen entity with idSourceSpecimen = %d' % entityId
        else:
            entity.sourceSpecimen = productionDatabase.sourceSpecimenDict[entityId]
            # type: int, name: idSourceSpecimen, foreignTable: SourceSpecimen, foreignColumn: idSourceSpecimen
            entity.sourceSpecimen.specimenSourceSpecimenList.append(entity)
        # many to one: project
        entityId = paftol.database.intOrNone(row[11])
        if entityId is None:
            entity.project = None
        elif entityId not in productionDatabase.projectDict:
            raise StandardError, 'no Project entity with idProject = %d' % entityId
        else:
            entity.project = productionDatabase.projectDict[entityId]
            # type: int, name: idProject, foreignTable: Project, foreignColumn: idProject
            entity.project.specimenProjectList.append(entity)
        # many to one: isoCountry
        entityId = paftol.database.intOrNone(row[12])
        if entityId is None:
            entity.isoCountry = None
        elif entityId not in productionDatabase.iSOCountryDict:
            raise StandardError, 'no ISOCountry entity with idIsoCountry = %d' % entityId
        else:
            entity.isoCountry = productionDatabase.iSOCountryDict[entityId]
            # type: int, name: idISOCountry, foreignTable: ISOCountry, foreignColumn: idISOCountry
            entity.isoCountry.specimenIsoCountryList.append(entity)
        # many to one: materialSource
        entityId = paftol.database.intOrNone(row[13])
        if entityId is None:
            entity.materialSource = None
        elif entityId not in productionDatabase.materialSourceDict:
            raise StandardError, 'no MaterialSource entity with idMaterialSource = %d' % entityId
        else:
            entity.materialSource = productionDatabase.materialSourceDict[entityId]
            # type: int, name: idMaterialSource, foreignTable: MaterialSource, foreignColumn: idMaterialSource
            entity.materialSource.specimenMaterialSourceList.append(entity)
        entity.ageOfMaterial = paftol.database.intOrNone(row[14])
        # many to one: herbarium
        entityId = paftol.database.intOrNone(row[15])
        if entityId is None:
            entity.herbarium = None
        elif entityId not in productionDatabase.herbariumDict:
            raise StandardError, 'no Herbarium entity with idHerbarium = %d' % entityId
        else:
            entity.herbarium = productionDatabase.herbariumDict[entityId]
            # type: int, name: idHerbarium, foreignTable: Herbarium, foreignColumn: idHerbarium
            entity.herbarium.specimenHerbariumList.append(entity)
        entity.specimenReference = paftol.database.strOrNone(row[16])
        entity.herbcatUrL = paftol.database.intOrNone(row[17])
        entity.blacklisted = paftol.database.intOrNone(row[18])
        entity.holdDate = paftol.database.strOrNone(row[19])
        # many to one: dataSource
        entityId = paftol.database.intOrNone(row[20])
        if entityId is None:
            entity.dataSource = None
        elif entityId not in productionDatabase.dataSourceDict:
            raise StandardError, 'no DataSource entity with idDataSource = %d' % entityId
        else:
            entity.dataSource = productionDatabase.dataSourceDict[entityId]
            # type: int, name: idDataSource, foreignTable: DataSource, foreignColumn: idDataSource
            entity.dataSource.specimenDataSourceList.append(entity)
        entityDict[entity.idSpecimen] = entity
    cursor.close()
    return entityDict


def loadSpecimenGeneStatsDict(connection, productionDatabase):
    cursor = connection.cursor()
    entityDict = {}
    sqlStatement = 'SELECT `idSpecimen`, `NumGene`, `RecoveredLength` FROM `SpecimenGeneStats`'
    cursor.execute(sqlStatement)
    for row in cursor:
        entity = SpecimenGeneStats()
        entity.idSpecimen = paftol.database.intOrNone(row[0])
        entity.numGene = paftol.database.intOrNone(row[1])
        entity.recoveredLength = paftol.database.floatOrNone(row[2])
        entityDict[entity.idSpecimen] = entity
    cursor.close()
    return entityDict


def loadSpecimenRawReadsDict(connection, productionDatabase):
    cursor = connection.cursor()
    entityDict = {}
    sqlStatement = 'SELECT `idSpecimenRawReads`, `idSpecimen`, `NumReads`, `SeqPlatform`, `ENAExpNum`, `ENARunNum` FROM `SpecimenRawReads`'
    cursor.execute(sqlStatement)
    for row in cursor:
        entity = SpecimenRawReads()
        entity.idSpecimenRawReads = paftol.database.intOrNone(row[0])
        entity.idSpecimen = paftol.database.intOrNone(row[1])
        entity.numReads = paftol.database.intOrNone(row[2])
        entity.seqPlatform = paftol.database.strOrNone(row[3])
        entity.enaExpNum = paftol.database.strOrNone(row[4])
        entity.enaRunNum = paftol.database.strOrNone(row[5])
        entityDict[entity.idSpecimenRawReads] = entity
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


def loadTargetSetDict(connection, productionDatabase):
    cursor = connection.cursor()
    entityDict = {}
    sqlStatement = 'SELECT `idTargetSet`, `TargetSetName`, `TargetsFastaFile`, `TargetsFastaFileMd5sum`, `NumTargetSequences` FROM `TargetSet`'
    cursor.execute(sqlStatement)
    for row in cursor:
        entity = TargetSet()
        entity.idTargetSet = paftol.database.intOrNone(row[0])
        entity.targetSetName = paftol.database.strOrNone(row[1])
        entity.targetsFastaFile = paftol.database.strOrNone(row[2])
        entity.targetsFastaFileMd5sum = paftol.database.strOrNone(row[3])
        entity.numTargetSequences = paftol.database.intOrNone(row[4])
        entityDict[entity.idTargetSet] = entity
    cursor.close()
    return entityDict


def loadTestResultDict(connection, productionDatabase):
    cursor = connection.cursor()
    entityDict = {}
    sqlStatement = 'SELECT `idTestResult`, `TestResult` FROM `TestResult`'
    cursor.execute(sqlStatement)
    for row in cursor:
        entity = TestResult()
        entity.idTestResult = paftol.database.intOrNone(row[0])
        entity.testResult = paftol.database.strOrNone(row[1])
        entityDict[entity.idTestResult] = entity
    cursor.close()
    return entityDict


def loadTrimmedRawFastqFileDict(connection, productionDatabase):
    cursor = connection.cursor()
    entityDict = {}
    sqlStatement = 'SELECT `idTrimmedRawFastqFile`, `idRawFastqFile`, `idSequenceRecovery`, `idFastqStats` FROM `TrimmedRawFastqFile`'
    cursor.execute(sqlStatement)
    for row in cursor:
        entity = TrimmedRawFastqFile()
        entity.idTrimmedRawFastqFile = paftol.database.intOrNone(row[0])
        # many to one: rawFastqFile
        entityId = paftol.database.intOrNone(row[1])
        if entityId is None:
            entity.rawFastqFile = None
        elif entityId not in productionDatabase.rawFastqFileDict:
            raise StandardError, 'no RawFastqFile entity with idRawFastqFile = %d' % entityId
        else:
            entity.rawFastqFile = productionDatabase.rawFastqFileDict[entityId]
            # type: int, name: idRawFastqFile, foreignTable: RawFastqFile, foreignColumn: idRawFastqFile
            entity.rawFastqFile.trimmedRawFastqFileRawFastqFileList.append(entity)
        # many to one: sequenceRecovery
        entityId = paftol.database.intOrNone(row[2])
        if entityId is None:
            entity.sequenceRecovery = None
        elif entityId not in productionDatabase.sequenceRecoveryDict:
            raise StandardError, 'no SequenceRecovery entity with idSequenceRecovery = %d' % entityId
        else:
            entity.sequenceRecovery = productionDatabase.sequenceRecoveryDict[entityId]
            # type: int, name: idSequenceRecovery, foreignTable: SequenceRecovery, foreignColumn: idSequenceRecovery
            entity.sequenceRecovery.trimmedRawFastqFileSequenceRecoveryList.append(entity)
        # many to one: fastqStats
        entityId = paftol.database.intOrNone(row[3])
        if entityId is None:
            entity.fastqStats = None
        elif entityId not in productionDatabase.fastqStatsDict:
            raise StandardError, 'no FastqStats entity with idFastqStats = %d' % entityId
        else:
            entity.fastqStats = productionDatabase.fastqStatsDict[entityId]
            # type: int, name: idFastqStats, foreignTable: FastqStats, foreignColumn: idFastqStats
            entity.fastqStats.trimmedRawFastqFileFastqStatsList.append(entity)
        entityDict[entity.idTrimmedRawFastqFile] = entity
    cursor.close()
    return entityDict


def loadValidationResultDict(connection, productionDatabase):
    cursor = connection.cursor()
    entityDict = {}
    sqlStatement = 'SELECT `idValidationResult`, `ValidationResult` FROM `ValidationResult`'
    cursor.execute(sqlStatement)
    for row in cursor:
        entity = ValidationResult()
        entity.idValidationResult = paftol.database.intOrNone(row[0])
        entity.validationResult = paftol.database.strOrNone(row[1])
        entityDict[entity.idValidationResult] = entity
    cursor.close()
    return entityDict


def loadWCVPFamilyDict(connection, productionDatabase):
    cursor = connection.cursor()
    entityDict = {}
    sqlStatement = 'SELECT `idFamily`, `Family`, `IPNIId`, `AcceptedIPNIId`, `APGIVId`, `ChecklistDB`, `PeerReviewed`, `TaxonStatusId`, `APGTaxonRemarks`, `IncludedInWCVP`, `idOrder` FROM `WCVPFamily`'
    cursor.execute(sqlStatement)
    for row in cursor:
        entity = WCVPFamily()
        entity.idFamily = paftol.database.intOrNone(row[0])
        entity.family = paftol.database.strOrNone(row[1])
        entity.ipniId = paftol.database.strOrNone(row[2])
        entity.acceptedIpniId = paftol.database.strOrNone(row[3])
        entity.apgivId = paftol.database.intOrNone(row[4])
        entity.checklistDb = paftol.database.strOrNone(row[5])
        entity.peerReviewed = paftol.database.strOrNone(row[6])
        entity.taxonStatusId = paftol.database.strOrNone(row[7])
        entity.apgTaxonRemarks = paftol.database.strOrNone(row[8])
        entity.includedInWcvP = paftol.database.strOrNone(row[9])
        # many to one: order
        entityId = paftol.database.intOrNone(row[10])
        if entityId is None:
            entity.order = None
        elif entityId not in productionDatabase.wCVPOrderDict:
            raise StandardError, 'no WCVPOrder entity with idOrder = %d' % entityId
        else:
            entity.order = productionDatabase.wCVPOrderDict[entityId]
            # type: int, name: idOrder, foreignTable: WCVPOrder, foreignColumn: idOrder
            entity.order.wcvpFamilyOrderList.append(entity)
        entityDict[entity.idFamily] = entity
    cursor.close()
    return entityDict


def loadWCVPHigherRankDict(connection, productionDatabase):
    cursor = connection.cursor()
    entityDict = {}
    sqlStatement = 'SELECT `Family`, `Kingdom`, `Phylum`, `Class`, `Order`, `ChecklistDB`, `PeerReviewed`, `IPNIId`, `TaxonStatusId`, `AcceptedIPNIId`, `APGIVId`, `APGIVOrderId`, `APGTaxonRemarks`, `IncludedInWCVP`, `Subclass`, `isAngiosperm` FROM `WCVPHigherRank`'
    cursor.execute(sqlStatement)
    for row in cursor:
        entity = WCVPHigherRank()
        entity.family = paftol.database.strOrNone(row[0])
        entity.kingdom = paftol.database.strOrNone(row[1])
        entity.phylum = paftol.database.strOrNone(row[2])
        entity.className = paftol.database.strOrNone(row[3])
        entity.order = paftol.database.strOrNone(row[4])
        entity.checklistDb = paftol.database.strOrNone(row[5])
        entity.peerReviewed = paftol.database.strOrNone(row[6])
        entity.ipniId = paftol.database.strOrNone(row[7])
        entity.taxonStatusId = paftol.database.strOrNone(row[8])
        entity.acceptedIpniId = paftol.database.strOrNone(row[9])
        entity.apgivId = paftol.database.intOrNone(row[10])
        entity.apgivOrderId = paftol.database.intOrNone(row[11])
        entity.apgTaxonRemarks = paftol.database.strOrNone(row[12])
        entity.includedInWcvP = paftol.database.strOrNone(row[13])
        entity.subclass = paftol.database.strOrNone(row[14])
        entity.isAngiosperm = paftol.database.intOrNone(row[15])
        entityDict[entity.family] = entity
    cursor.close()
    return entityDict


def loadWCVPLoadDict(connection, productionDatabase):
    cursor = connection.cursor()
    entityDict = {}
    sqlStatement = 'SELECT `PlantNameId`, `IPNIId`, `TaxonRank`, `TaxonStatus`, `Family`, `GenusHybrid`, `Genus`, `SpeciesHybrid`, `Species`, `InfraspecificRank`, `Infraspecies`, `ParentheticalAuthor`, `PrimaryAuthor`, `PublicationAuthor`, `PlaceOfPublication`, `VolumeAndPage`, `FirstPublished`, `NomenclaturalRemarks`, `GeographicArea`, `LifeformDescription`, `ClimateDescription`, `TaxonName`, `TaxonAuthors`, `AcceptedPlantName`, `BasionymPlantNameId`, `ReplacedSynonymAuthor`, `HomotypicSynonym`, `ParentPlantName`, `POWOId`, `HybridFormula`, `Reviewed` FROM `WCVPLoad`'
    cursor.execute(sqlStatement)
    for row in cursor:
        entity = WCVPLoad()
        entity.plantNameId = paftol.database.strOrNone(row[0])
        entity.ipniId = paftol.database.strOrNone(row[1])
        entity.taxonRank = paftol.database.strOrNone(row[2])
        entity.taxonStatus = paftol.database.strOrNone(row[3])
        entity.family = paftol.database.strOrNone(row[4])
        entity.genusHybrid = paftol.database.strOrNone(row[5])
        entity.genus = paftol.database.strOrNone(row[6])
        entity.speciesHybrid = paftol.database.strOrNone(row[7])
        entity.species = paftol.database.strOrNone(row[8])
        entity.infraspecificRank = paftol.database.strOrNone(row[9])
        entity.infraspecies = paftol.database.strOrNone(row[10])
        entity.parentheticalAuthor = paftol.database.strOrNone(row[11])
        entity.primaryAuthor = paftol.database.strOrNone(row[12])
        entity.publicationAuthor = paftol.database.strOrNone(row[13])
        entity.placeOfPublication = paftol.database.strOrNone(row[14])
        entity.volumeAndPage = paftol.database.strOrNone(row[15])
        entity.firstPublished = paftol.database.strOrNone(row[16])
        entity.nomenclaturalRemarks = paftol.database.strOrNone(row[17])
        entity.geographicArea = paftol.database.strOrNone(row[18])
        entity.lifeformDescription = paftol.database.strOrNone(row[19])
        entity.climateDescription = paftol.database.strOrNone(row[20])
        entity.taxonName = paftol.database.strOrNone(row[21])
        entity.taxonAuthors = paftol.database.strOrNone(row[22])
        entity.acceptedPlantName = paftol.database.strOrNone(row[23])
        entity.basionymPlantNameId = paftol.database.strOrNone(row[24])
        entity.replacedSynonymAuthor = paftol.database.strOrNone(row[25])
        entity.homotypicSynonym = paftol.database.strOrNone(row[26])
        entity.parentPlantName = paftol.database.strOrNone(row[27])
        entity.powoId = paftol.database.strOrNone(row[28])
        entity.hybridFormula = paftol.database.strOrNone(row[29])
        entity.reviewed = paftol.database.strOrNone(row[30])
        entityDict[entity.plantNameId] = entity
    cursor.close()
    return entityDict


def loadWCVPNameDict(connection, productionDatabase):
    cursor = connection.cursor()
    entityDict = {}
    sqlStatement = 'SELECT `PlantNameId`, `IPNIId`, `TaxonRank`, `TaxonStatus`, `GenusHybrid`, `Genus`, `SpeciesHybrid`, `Species`, `InfraspecificRank`, `Infraspecies`, `ParentheticalAuthor`, `PrimaryAuthor`, `PublicationAuthor`, `PlaceOfPublication`, `VolumeAndPage`, `FirstPublished`, `NomenclaturalRemarks`, `GeographicArea`, `LifeformDescription`, `ClimateDescription`, `TaxonName`, `TaxonAuthors`, `AcceptedPlantName`, `BasionymPlantNameId`, `ReplacedSynonymAuthor`, `HomotypicSynonym`, `ParentPlantName`, `POWOId`, `HybridFormula`, `Reviewed`, `idFamily` FROM `WCVPName`'
    cursor.execute(sqlStatement)
    for row in cursor:
        entity = WCVPName()
        entity.plantNameId = paftol.database.strOrNone(row[0])
        entity.ipniId = paftol.database.strOrNone(row[1])
        entity.taxonRank = paftol.database.strOrNone(row[2])
        entity.taxonStatus = paftol.database.strOrNone(row[3])
        entity.genusHybrid = paftol.database.strOrNone(row[4])
        entity.genus = paftol.database.strOrNone(row[5])
        entity.speciesHybrid = paftol.database.strOrNone(row[6])
        entity.species = paftol.database.strOrNone(row[7])
        entity.infraspecificRank = paftol.database.strOrNone(row[8])
        entity.infraspecies = paftol.database.strOrNone(row[9])
        entity.parentheticalAuthor = paftol.database.strOrNone(row[10])
        entity.primaryAuthor = paftol.database.strOrNone(row[11])
        entity.publicationAuthor = paftol.database.strOrNone(row[12])
        entity.placeOfPublication = paftol.database.strOrNone(row[13])
        entity.volumeAndPage = paftol.database.strOrNone(row[14])
        entity.firstPublished = paftol.database.strOrNone(row[15])
        entity.nomenclaturalRemarks = paftol.database.strOrNone(row[16])
        entity.geographicArea = paftol.database.strOrNone(row[17])
        entity.lifeformDescription = paftol.database.strOrNone(row[18])
        entity.climateDescription = paftol.database.strOrNone(row[19])
        entity.taxonName = paftol.database.strOrNone(row[20])
        entity.taxonAuthors = paftol.database.strOrNone(row[21])
        entity.acceptedPlantName = paftol.database.strOrNone(row[22])
        entity.basionymPlantNameId = paftol.database.strOrNone(row[23])
        entity.replacedSynonymAuthor = paftol.database.strOrNone(row[24])
        entity.homotypicSynonym = paftol.database.strOrNone(row[25])
        entity.parentPlantName = paftol.database.strOrNone(row[26])
        entity.powoId = paftol.database.strOrNone(row[27])
        entity.hybridFormula = paftol.database.strOrNone(row[28])
        entity.reviewed = paftol.database.strOrNone(row[29])
        # many to one: family
        entityId = paftol.database.intOrNone(row[30])
        if entityId is None:
            entity.family = None
        elif entityId not in productionDatabase.wCVPFamilyDict:
            raise StandardError, 'no WCVPFamily entity with idFamily = %d' % entityId
        else:
            entity.family = productionDatabase.wCVPFamilyDict[entityId]
            # type: int, name: idFamily, foreignTable: WCVPFamily, foreignColumn: idFamily
            entity.family.wcvpNameFamilyList.append(entity)
        entityDict[entity.powoId] = entity
    cursor.close()
    return entityDict


def loadWCVPOrderDict(connection, productionDatabase):
    cursor = connection.cursor()
    entityDict = {}
    sqlStatement = 'SELECT `idOrder`, `Order`, `APGIVOrderId`, `isAngiosperm`, `Subclass`, `Class`, `Phylum`, `Kingdom` FROM `WCVPOrder`'
    cursor.execute(sqlStatement)
    for row in cursor:
        entity = WCVPOrder()
        entity.idOrder = paftol.database.intOrNone(row[0])
        entity.order = paftol.database.strOrNone(row[1])
        entity.apgivOrderId = paftol.database.intOrNone(row[2])
        entity.isAngiosperm = paftol.database.intOrNone(row[3])
        entity.subclass = paftol.database.strOrNone(row[4])
        entity.className = paftol.database.strOrNone(row[5])
        entity.phylum = paftol.database.strOrNone(row[6])
        entity.kingdom = paftol.database.strOrNone(row[7])
        entityDict[entity.idOrder] = entity
    cursor.close()
    return entityDict


class ProductionDatabase(object):

    def __init__(self, connection):
        self.actionDict = {}
        self.additionalBaitKitDict = {}
        self.alternativeNameDict = {}
        self.blacklistedReasonDict = {}
        self.coordinatesDict = {}
        self.dNAVolumeDict = {}
        self.dataReleaseDict = {}
        self.dataSourceDict = {}
        self.dataSourceAndMethodInReleaseDict = {}
        self.decisionReasonDict = {}
        self.exemplarGeneDict = {}
        self.extractionTypeDict = {}
        self.fastqStatsDict = {}
        self.geneDict = {}
        self.geneStatsDict = {}
        self.geneTreeDict = {}
        self.geneTypeDict = {}
        self.herbariumDict = {}
        self.iSOCountryDict = {}
        self.indexesDict = {}
        self.libraryDict = {}
        self.locationDict = {}
        self.logDict = {}
        self.materialSourceDict = {}
        self.mergedSequenceDict = {}
        self.migrationsLogDict = {}
        self.nameSourceDict = {}
        self.platformDict = {}
        self.projectDict = {}
        self.qualityDict = {}
        self.rawFastaFileDict = {}
        self.rawFastqFileDict = {}
        self.recoveredContigDict = {}
        self.recoveryMethodDict = {}
        self.referenceTargetDict = {}
        self.sampleDict = {}
        self.sequenceDict = {}
        self.sequenceDataReleaseDict = {}
        self.sequenceGeneStatsDict = {}
        self.sequenceRawReadsDict = {}
        self.sequenceRecoveryDict = {}
        self.sequencingStrategyDict = {}
        self.sourceDict = {}
        self.sourceSpecimenDict = {}
        self.speciesTreeDict = {}
        self.specimenDict = {}
        self.specimenGeneStatsDict = {}
        self.specimenRawReadsDict = {}
        self.statusDict = {}
        self.targetSetDict = {}
        self.testResultDict = {}
        self.trimmedRawFastqFileDict = {}
        self.validationResultDict = {}
        self.wCVPFamilyDict = {}
        self.wCVPHigherRankDict = {}
        self.wCVPLoadDict = {}
        self.wCVPNameDict = {}
        self.wCVPOrderDict = {}
        self.actionDict = loadActionDict(connection, self)
        self.additionalBaitKitDict = loadAdditionalBaitKitDict(connection, self)
        self.wCVPOrderDict = loadWCVPOrderDict(connection, self)
        self.wCVPFamilyDict = loadWCVPFamilyDict(connection, self)
        self.wCVPNameDict = loadWCVPNameDict(connection, self)
        self.sourceSpecimenDict = loadSourceSpecimenDict(connection, self)
        self.projectDict = loadProjectDict(connection, self)
        self.iSOCountryDict = loadISOCountryDict(connection, self)
        self.materialSourceDict = loadMaterialSourceDict(connection, self)
        self.herbariumDict = loadHerbariumDict(connection, self)
        self.dataSourceDict = loadDataSourceDict(connection, self)
        self.specimenDict = loadSpecimenDict(connection, self)
        self.nameSourceDict = loadNameSourceDict(connection, self)
        self.alternativeNameDict = loadAlternativeNameDict(connection, self)
        self.blacklistedReasonDict = loadBlacklistedReasonDict(connection, self)
        self.coordinatesDict = loadCoordinatesDict(connection, self)
        self.dNAVolumeDict = loadDNAVolumeDict(connection, self)
        self.dataReleaseDict = loadDataReleaseDict(connection, self)
        self.targetSetDict = loadTargetSetDict(connection, self)
        self.recoveryMethodDict = loadRecoveryMethodDict(connection, self)
        self.dataSourceAndMethodInReleaseDict = loadDataSourceAndMethodInReleaseDict(connection, self)
        self.decisionReasonDict = loadDecisionReasonDict(connection, self)
        self.exemplarGeneDict = loadExemplarGeneDict(connection, self)
        self.extractionTypeDict = loadExtractionTypeDict(connection, self)
        self.fastqStatsDict = loadFastqStatsDict(connection, self)
        self.geneTypeDict = loadGeneTypeDict(connection, self)
        self.geneDict = loadGeneDict(connection, self)
        self.geneStatsDict = loadGeneStatsDict(connection, self)
        self.speciesTreeDict = loadSpeciesTreeDict(connection, self)
        self.geneTreeDict = loadGeneTreeDict(connection, self)
        self.indexesDict = loadIndexesDict(connection, self)
        self.qualityDict = loadQualityDict(connection, self)
        self.sampleDict = loadSampleDict(connection, self)
        self.statusDict = loadStatusDict(connection, self)
        self.libraryDict = loadLibraryDict(connection, self)
        self.locationDict = loadLocationDict(connection, self)
        self.logDict = loadLogDict(connection, self)
        self.platformDict = loadPlatformDict(connection, self)
        self.sequencingStrategyDict = loadSequencingStrategyDict(connection, self)
        self.sequenceDict = loadSequenceDict(connection, self)
        self.mergedSequenceDict = loadMergedSequenceDict(connection, self)
        self.migrationsLogDict = loadMigrationsLogDict(connection, self)
        self.rawFastaFileDict = loadRawFastaFileDict(connection, self)
        self.rawFastqFileDict = loadRawFastqFileDict(connection, self)
        self.sequenceRecoveryDict = loadSequenceRecoveryDict(connection, self)
        self.referenceTargetDict = loadReferenceTargetDict(connection, self)
        self.recoveredContigDict = loadRecoveredContigDict(connection, self)
        self.testResultDict = loadTestResultDict(connection, self)
        self.validationResultDict = loadValidationResultDict(connection, self)
        self.sequenceDataReleaseDict = loadSequenceDataReleaseDict(connection, self)
        self.sequenceGeneStatsDict = loadSequenceGeneStatsDict(connection, self)
        self.sequenceRawReadsDict = loadSequenceRawReadsDict(connection, self)
        self.sourceDict = loadSourceDict(connection, self)
        self.specimenGeneStatsDict = loadSpecimenGeneStatsDict(connection, self)
        self.specimenRawReadsDict = loadSpecimenRawReadsDict(connection, self)
        self.trimmedRawFastqFileDict = loadTrimmedRawFastqFileDict(connection, self)
        self.wCVPHigherRankDict = loadWCVPHigherRankDict(connection, self)
        self.wCVPLoadDict = loadWCVPLoadDict(connection, self)

    def __str__(self):
        s = ''
        s = s + 'action: %d\n' % len(self.actionDict)
        s = s + 'additionalBaitKit: %d\n' % len(self.additionalBaitKitDict)
        s = s + 'wCVPOrder: %d\n' % len(self.wCVPOrderDict)
        s = s + 'wCVPFamily: %d\n' % len(self.wCVPFamilyDict)
        s = s + 'wCVPName: %d\n' % len(self.wCVPNameDict)
        s = s + 'sourceSpecimen: %d\n' % len(self.sourceSpecimenDict)
        s = s + 'project: %d\n' % len(self.projectDict)
        s = s + 'iSOCountry: %d\n' % len(self.iSOCountryDict)
        s = s + 'materialSource: %d\n' % len(self.materialSourceDict)
        s = s + 'herbarium: %d\n' % len(self.herbariumDict)
        s = s + 'dataSource: %d\n' % len(self.dataSourceDict)
        s = s + 'specimen: %d\n' % len(self.specimenDict)
        s = s + 'nameSource: %d\n' % len(self.nameSourceDict)
        s = s + 'alternativeName: %d\n' % len(self.alternativeNameDict)
        s = s + 'blacklistedReason: %d\n' % len(self.blacklistedReasonDict)
        s = s + 'coordinates: %d\n' % len(self.coordinatesDict)
        s = s + 'dNAVolume: %d\n' % len(self.dNAVolumeDict)
        s = s + 'dataRelease: %d\n' % len(self.dataReleaseDict)
        s = s + 'targetSet: %d\n' % len(self.targetSetDict)
        s = s + 'recoveryMethod: %d\n' % len(self.recoveryMethodDict)
        s = s + 'dataSourceAndMethodInRelease: %d\n' % len(self.dataSourceAndMethodInReleaseDict)
        s = s + 'decisionReason: %d\n' % len(self.decisionReasonDict)
        s = s + 'exemplarGene: %d\n' % len(self.exemplarGeneDict)
        s = s + 'extractionType: %d\n' % len(self.extractionTypeDict)
        s = s + 'fastqStats: %d\n' % len(self.fastqStatsDict)
        s = s + 'geneType: %d\n' % len(self.geneTypeDict)
        s = s + 'gene: %d\n' % len(self.geneDict)
        s = s + 'geneStats: %d\n' % len(self.geneStatsDict)
        s = s + 'speciesTree: %d\n' % len(self.speciesTreeDict)
        s = s + 'geneTree: %d\n' % len(self.geneTreeDict)
        s = s + 'indexes: %d\n' % len(self.indexesDict)
        s = s + 'quality: %d\n' % len(self.qualityDict)
        s = s + 'sample: %d\n' % len(self.sampleDict)
        s = s + 'status: %d\n' % len(self.statusDict)
        s = s + 'library: %d\n' % len(self.libraryDict)
        s = s + 'location: %d\n' % len(self.locationDict)
        s = s + 'log: %d\n' % len(self.logDict)
        s = s + 'platform: %d\n' % len(self.platformDict)
        s = s + 'sequencingStrategy: %d\n' % len(self.sequencingStrategyDict)
        s = s + 'sequence: %d\n' % len(self.sequenceDict)
        s = s + 'mergedSequence: %d\n' % len(self.mergedSequenceDict)
        s = s + 'migrationsLog: %d\n' % len(self.migrationsLogDict)
        s = s + 'rawFastaFile: %d\n' % len(self.rawFastaFileDict)
        s = s + 'rawFastqFile: %d\n' % len(self.rawFastqFileDict)
        s = s + 'sequenceRecovery: %d\n' % len(self.sequenceRecoveryDict)
        s = s + 'referenceTarget: %d\n' % len(self.referenceTargetDict)
        s = s + 'recoveredContig: %d\n' % len(self.recoveredContigDict)
        s = s + 'testResult: %d\n' % len(self.testResultDict)
        s = s + 'validationResult: %d\n' % len(self.validationResultDict)
        s = s + 'sequenceDataRelease: %d\n' % len(self.sequenceDataReleaseDict)
        s = s + 'sequenceGeneStats: %d\n' % len(self.sequenceGeneStatsDict)
        s = s + 'sequenceRawReads: %d\n' % len(self.sequenceRawReadsDict)
        s = s + 'source: %d\n' % len(self.sourceDict)
        s = s + 'specimenGeneStats: %d\n' % len(self.specimenGeneStatsDict)
        s = s + 'specimenRawReads: %d\n' % len(self.specimenRawReadsDict)
        s = s + 'trimmedRawFastqFile: %d\n' % len(self.trimmedRawFastqFileDict)
        s = s + 'wCVPHigherRank: %d\n' % len(self.wCVPHigherRankDict)
        s = s + 'wCVPLoad: %d\n' % len(self.wCVPLoadDict)
        return s

