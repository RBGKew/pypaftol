import sys
import re
import copy
import os
import os.path
import logging
import unicodedata

import mysql.connector

import paftol
import paftol.tools
import paftol.database.analysis
import paftol.database.production


logger = logging.getLogger(__name__)


def strOrNone(x):
    if x is None:
        return None
    else:
        return unicodedata.normalize('NFKD', x).encode('ascii', 'ignore').strip()

    
def intOrNone(x):
    if x is None:
        return None
    else:
        return int(x)
    

def floatOrNone(x):
    if x is None:
        return None
    else:
        return float(x)
    

def findFastaFile(analysisDatabase, fastaFname):
    for fastaFile in analysisDatabase.fastaFileDict.values():
        if fastaFile.filename == fastaFname:
            return fastaFile
    return None


def findFastqFile(analysisDatabase, fastqFname):
    for fastqFile in analysisDatabase.fastqFileDict.values():
        if fastqFile.filename == fastqFname:
            return fastqFile
    return None


class PaftolDatabaseDetails(object):

    reDbusername = re.compile('username: *([^ ]+)')
    reDbpassword = re.compile('password: *([^ ]+)')
    reDbhost = re.compile('host: *([^ ]+)')
    reDbname = re.compile('dbname: *([^ ]+)')

    def __init__(self, detailsFile=None):
        if detailsFile is None:
            self.dbusername = None
            self.dbpassword = None
            self.dbhost = None
            self.dbname = None
        else:
            self.readFile(detailsFile)

    def readDetailsLine(self, detailsFile, detailsRe, errorMsg):
        line = detailsFile.readline()
        # sys.stderr.write('got line: "%s"\n' % repr(line))
        m = detailsRe.match(line.strip())
        if m is None:
            raise StandardError, errorMsg
        return m.group(1)

    def readFile(self, detailsFile):
        self.dbusername = self.readDetailsLine(detailsFile, self.reDbusername, 'malformed dbusername line')
        self.dbpassword = self.readDetailsLine(detailsFile, self.reDbpassword, 'malformed dbpassword line')
        self.dbhost = self.readDetailsLine(detailsFile, self.reDbhost, 'malformed dbhost line')
        self.dbname = self.readDetailsLine(detailsFile, self.reDbname, 'malformed dbname line')

    def makeConnection(self):
        return mysql.connector.connection.MySQLConnection(user=self.dbusername, password=self.dbpassword, host=self.dbhost, database=self.dbname)


def getDatabaseDetails(detailsFname):
    databaseDetails = None
    with open(detailsFname, 'r') as f:
        databaseDetails = PaftolDatabaseDetails(f)
    return databaseDetails


def getProductionDatabaseDetails(detailsFname=None):
    if detailsFname is None:
        detailsFname = os.path.join(os.environ['HOME'], '.paftol', 'productiondb.cfg')
    return getDatabaseDetails(detailsFname)


def getAnalysisDatabaseDetails(detailsFname=None):
    if detailsFname is None:
        detailsFname = os.path.join(os.environ['HOME'], '.paftol', 'analysisdb.cfg')
    return getDatabaseDetails(detailsFname)


def getProductionDatabase(detailsFname=None):
    productionDatabaseDetails = getProductionDatabaseDetails(detailsFname)
    connection = productionDatabaseDetails.makeConnection()
    productionDatabase = paftol.database.production.ProductionDatabase(connection)
    return productionDatabase


def getAnalysisDatabase(detailsFname=None):
    analysisDatabaseDetails = getAnalysisDatabaseDetails(detailsFname)
    connection = analysisDatabaseDetails.makeConnection()
    analysisDatabase = paftol.database.analysis.AnalysisDatabase(connection)
    return analysisDatabase


def matchesExpectedFastqFname(fastqFname, sequence):
    if sequence.r1FastqFile is None or sequence.r2FastqFile is None:
        return False
    fastqBasename = os.path.basename(fastqFname)
    return fastqBasename == sequence.r1FastqFile or fastqBasename == sequence.r2FastqFile


def matchesExpectedSequencingRun(sequence, sequencingPoolNumber):
    return sequence.sequencingRun is not None and sequence.sequencingRun.upper() == 'SP%04d' % sequencingPoolNumber

    
def findMatchingSequenceList(productionDatabase, fastqFname, sequencingPoolNumber):
    sequenceList = []
    for sequence in productionDatabase.sequenceDict.values():
        if matchesExpectedFastqFname(fastqFname, sequence):
            if matchesExpectedSequencingRun(sequence, sequencingPoolNumber):
                sequenceList.append(sequence)
    return sequenceList


class ExistingFastqFile(object):

    # FIXME: regular expressions used to search along entire path, spurious matches not impossible
    spNumberRe = re.compile('SP([0-9][0-9][0-9][0-9])((-|_+)[A-Z][A-Za-z0-9_ ()-[\\]]+)?$')
    paftolPrefixedFastqFnameRe = re.compile('PAFTOL[-_]([0-9]+)_R[12]_[0-9]+\\.fastq')

    def __init__(self, rawFastqFname):
        self.rawFastqFname = rawFastqFname

    def findSequencingPoolNumber(self):
        dirname, basename = os.path.split(self.rawFastqFname)
        if dirname != '':
            d, spDirname = os.path.split(dirname)
            sys.stderr.write('spDirname: %s\n' % spDirname)
            m = self.spNumberRe.match(spDirname)
            if m is not None:
                return int(m.group(1))
        return None
    
    def findPaftolPrefixedNumber(self):
        m = self.paftolPrefixedFastqFnameRe.search(self.rawFastqFname)
        if m is not None:
            return int(m.group(1))
        return None


def findSequenceForFastqFname(productionDatabase, fastqFname):
    logList = []
    existingFastqFile = ExistingFastqFile(fastqFname)
    paftolPrefixedNumber = existingFastqFile.findPaftolPrefixedNumber()
    sequencingPoolNumber = existingFastqFile.findSequencingPoolNumber()
    logList.append('fastqFname: %s' % fastqFname)
    logList.append('paftolPrefixedNumber: %s' % paftolPrefixedNumber)
    logList.append('sequencingPoolNumber: %s' % sequencingPoolNumber)
    # logger.debug('raw: %s, paftolPrefixedNumber: %s, spNumber: %s', existingFastqFile.rawFastqFname, paftolPrefixedNumber, sequencingPoolNumber)
    paftolPrefixedNumber = existingFastqFile.findPaftolPrefixedNumber()
    if paftolPrefixedNumber is not None:
        idSequencing = paftolPrefixedNumber
        if idSequencing in productionDatabase.sequenceDict:
            sequence = productionDatabase.sequenceDict[idSequencing]
            if matchesExpectedFastqFname(fastqFname, sequence):
                if sequencingPoolNumber is not None:
                    if matchesExpectedSequencingRun(sequence, sequencingPoolNumber):
                        return sequence
                    else:
                        logList.append('found sequencingRun %s, not consistent with sequencingPoolNumber %d' % (sequence.sequencingRun, sequencingPoolNumber))
            else:
                logList.append('fastqFname %s does not match expected names %s or %s' % (fastqFname, sequence.r1FastqFile, sequence.r2FastqFile))
        else:
            logList.append('no sequence with idSequencing %d' % idSequencing)
    else:
        if sequencingPoolNumber is not None:
            sequenceList = findMatchingSequenceList(productionDatabase, fastqFname, sequencingPoolNumber)
            if len(sequenceList) == 0:
                logList.append('no match by fname found')
            elif len(sequenceList) == 1:
                sequence = sequenceList[0]
                return sequence
            else:
                logList.append('multiple matches: %s' % ', '.join(['%d' % s.idSequencing for s in sequenceList]))
        else:
            logList.append('no sequencingPoolNumber')
    logList.append('unresolved')
    logger.debug(', '.join(logList))
    return None


def findSequenceForFastqFnameOld(productionDatabase, fastqFname):
    paftolPrefixedFastqFnameRe = re.compile('PAFTOL-([0-9]+)_R[12]_[0-9]+\\.fastq')
    spNumberRe = re.compile('(SP[0-9][0-9][0-9][0-9])([^/]*)/([0-9]+).*\\.fastq')
    m = paftolPrefixedFastqFnameRe.match(fastqFname)
    if m is not None:
        idSequencing = int(m.group(1))
        if idSequencing in productionDatabase.sequenceDict:
            sequence = productionDatabase.sequenceDict[idSequencing]
            if matchesExpectedFastqFname(fastqFname, sequence):
                return sequence
        return None
    m = spNumberRe.search(fastqFname)
    if m is not None:
        sequencingRun = m.group(1)
        i = int(m.group(3))
        if i in productionDatabase.sequenceDict:
            sequence = productionDatabase.sequenceDict[i]
            if sequence.sequencingRun == sequencingRun and matchesExpectedFastqFname(fastqFname, sequence):
                return sequence
        for sequence in productionDatabase.sequenceDict.values():
            if sequence.library is not None and sequence.library.sample is not None and sequence.library.sample.specimen is not None and sequence.library.sample.specimen.idPaftol is not None and sequence.library.sample.specimen.idPaftol == i and matchesExpectedFastqFname(fastqFname, sequence):
                return sequence
        return None
    # raise StandardError, 'malformed fastqFname: %s' % fastqFname
    return None


def canonicalSymlinkName(sequence, orientation, gzipped):
    # sys.stderr.write('sequence id=%d\n' % sequence.idSequencing)
    gzExt = ''
    if gzipped:
        gzExt = '.gz'
    return 'PAFTOL_%06d_R%1d.fastq%s' % (sequence.idSequencing, orientation, gzExt)


def parseCanonicalSymlink(symlinkName):
    print 'symlinkName: ', symlinkName
    symlinkRe = re.compile('PAFTOL_([0-9]+)_R([12]).fastq')
    m = symlinkRe.match(symlinkName)
    if m is not None:
        return int(m.group(1)), int(m.group(2))
    return None, None


def makeSymlink(symlinkDirname, sequence, fastqFname):
    """Set up a canonical symlink for a sequence in the specified directory.
    
If a symlink or other file with the computed canonical name already exists, no
symlink is created and a warning is issued.

@param symlinkDirname: directory in which to create the symlink
@type symlinkDirname: C{str}
@param sequence: the sequence for which to generate the symlink
@type sequence: C{paftol.database.production.Sequence} instance
@param fastqFname: name of the fastq file
@type fastqFname: C{str}
    """
    orientation = paftol.tools.fastqOrientation(fastqFname)
    gzipped = paftol.tools.isGzipped(fastqFname)
    symlinkName = canonicalSymlinkName(sequence, orientation, gzipped)
    symlinkPath = os.path.join(symlinkDirname, symlinkName)
    if os.path.lexists(symlinkPath) or os.path.exists(symlinkPath):
        logger.warning('sequence %d: link %s already exists', sequence.idSequencing, symlinkPath)
    else:
        os.symlink(fastqFname, symlinkPath)


def generateUnusedPrimaryKey(cursor, tableName, primaryKeyColumnName='id'):
    sqlStatement = 'SELECT max(%s) FROM %s' % (primaryKeyColumnName, tableName)
    cursor.execute(sqlStatement)
    row = cursor.fetchone()
    maxPk = 0
    if row is not None and row[0] is not None:
        maxPk = int(row[0])
    return maxPk + 1


def insertGene(connection, geneName, geneTypeId):
    lockCursor = connection.cursor()
    lockCursor.execute('LOCK TABLE PaftolGene WRITE')
    try:
        cursor = connection.cursor(prepared=True)
        try:
            paftolGeneId = generateUnusedPrimaryKey(cursor, 'PaftolGene')
            cursor.execute('INSERT INTO PaftolGene (id, geneName, geneTypeId) VALUES (%s, %s, %s)', (paftolGeneId, geneName, geneTypeId, ))
        finally:
            cursor.close()
    finally:
        lockCursor.execute('UNLOCK TABLES')
        lockCursor.close()
    return paftolGeneId


def insertFastaFile(connection, fastaFname, dirname=None):
    fastaPath = fastaFname
    if dirname is not None:
        fastaPath = os.path.join(dirname, fastaFname)
    md5 = paftol.tools.md5HexdigestFromFile(fastaPath)
    numSequences = len(paftol.tools.fastaSeqRecordList(fastaPath))
    lockCursor = connection.cursor()
    lockCursor.execute('LOCK TABLE FastaFile WRITE')
    try:
        cursor = connection.cursor(prepared=True)
        try:
            fastaFileId = generateUnusedPrimaryKey(cursor, 'FastaFile')
            cursor.execute('INSERT INTO FastaFile (id, filename, md5sum, numSequences) VALUES (%s, %s, %s, %s)', (fastaFileId, fastaFname, md5, numSequences, ))
        finally:
            cursor.close()
    finally:
        lockCursor.execute('UNLOCK TABLES')
        lockCursor.close()
    return fastaFileId


def insertFastaExternalAccession(connection, fastaFilename, dataOriginAcronym, accession):
    lockCursor = connection.cursor()
    lockCursor.execute('LOCK TABLE ExternalAccesion WRITE')
    try:
        cursor = connection.cursor(prepared=True)
        try:
            fastaFileId = generateUnusedPrimaryKey(cursor, 'FastaFile')
            cursor.execute('INSERT INTO FastaFile (id, filename, md5sum, numSequences) VALUES (%s, %s, %s, %s)', (fastaFileId, fastaFname, md5, numSequences, ))
        finally:
            cursor.close()
    finally:
        lockCursor.execute('UNLOCK TABLES')
        lockCursor.close()
    return fastaFileId


def addFastqStats(connection, fastqcStats):
    logger.debug('starting')
    fastqcSummaryStats = paftol.tools.FastqcSummaryStats(fastqcStats)
    cursor = connection.cursor(prepared=True)
    try:
        fastqStatsId = generateUnusedPrimaryKey(cursor, 'FastqStats')
        logger.debug('meanAdapterContent: %s, maxAdapterContent: %s', fastqcSummaryStats.meanAdapterContent, fastqcSummaryStats.maxAdapterContent)        
        cursor.execute('INSERT INTO FastqStats (id, numReads, qual28, meanA, meanC, meanG, meanT, stddevA, stddevC, stddevG, stddevT, meanN, stddevN, meanAdapterContent, maxAdapterContent) VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s)', (fastqStatsId, fastqcSummaryStats.numReads, fastqcSummaryStats.qual28, fastqcSummaryStats.meanA, fastqcSummaryStats.meanC, fastqcSummaryStats.meanG, fastqcSummaryStats.meanT, fastqcSummaryStats.stddevA, fastqcSummaryStats.stddevC, fastqcSummaryStats.stddevG, fastqcSummaryStats.stddevT, fastqcSummaryStats.meanN, fastqcSummaryStats.stddevN, fastqcSummaryStats.meanAdapterContent, fastqcSummaryStats.maxAdapterContent))
    finally:
        cursor.close()
    return fastqStatsId


def preInsertCheckPaftolFastqFile(paftolFastqFile):
    if paftolFastqFile.id is not None:
        raise StandardError, 'illegal state: PaftolFastqFile instance has id %d, not None' % paftolFastqFile.id
    if paftolFastqFile.fastqFile is None:
        raise StandardError, 'illegal state: PaftolFastqFile has fastqFile attribute set to None'
    if paftolFastqFile.fastqFile.fastqStats is not None and  paftolFastqFile.fastqFile.fastqStats.id is not None:
        raise StandardError, 'illegal state: new FastqFile %s has existing FastqStats %d' % (paftolFastqFile.fastqFile.filename, paftolFastqFile.fastqFile.fastqStats.id)


def insertPaftolFastqFileList(connection, paftolFastqFileList):
    for paftolFastqFile in paftolFastqFileList:
        preInsertCheckPaftolFastqFile(paftolFastqFile)
    insertedPaftolFastqFileList = []
    transactionSuccessful = False
    lockCursor = connection.cursor()
    lockCursor.execute('LOCK TABLE PaftolFastqFile WRITE, FastqFile WRITE, FastqStats WRITE')
    try:
        cursor = connection.cursor(prepared=True)
        try:
            for paftolFastqFile in paftolFastqFileList:
                insertedPaftolFastqFile = copy.deepcopy(paftolFastqFile)
                insertedPaftolFastqFile.id = generateUnusedPrimaryKey(cursor, 'PaftolFastqFile')
                insertedPaftolFastqFile.fastqFile.id = generateUnusedPrimaryKey(cursor, 'FastqFile')
                if insertedPaftolFastqFile.fastqFile.fastqStats is not None:
                    insertedPaftolFastqFile.fastqFile.fastqStats.id = generateUnusedPrimaryKey(cursor, 'FastqStats')
                    insertedPaftolFastqFile.fastqFile.fastqStats.insertIntoDatabase(cursor)
                insertedPaftolFastqFile.fastqFile.insertIntoDatabase(cursor)
                insertedPaftolFastqFile.insertIntoDatabase(cursor)
            connection.commit()
            transactionSuccessful = True
        finally:
            if not transactionSuccessful:
                connection.rollback()
            cursor.close()
    finally:
        if not transactionSuccessful:
            connection.rollback()
        lockCursor.execute('UNLOCK TABLES')
        lockCursor.close()
    return transactionSuccessful


def fastqStatsFromFastqcStats(fastqcStats):
    fastqcSummaryStats = paftol.tools.FastqcSummaryStats(fastqcStats)
    return paftol.database.analysis.FastqStats(None, numReads=fastqcSummaryStats.numReads, qual28=fastqcSummaryStats.qual28, meanA=fastqcSummaryStats.meanA, meanC=fastqcSummaryStats.meanC, meanG=fastqcSummaryStats.meanG, meanT=fastqcSummaryStats.meanT, stddevA=fastqcSummaryStats.stddevA, stddevC=fastqcSummaryStats.stddevC, stddevG=fastqcSummaryStats.stddevG, stddevT=fastqcSummaryStats.stddevT, meanN=fastqcSummaryStats.meanN, stddevN=fastqcSummaryStats.stddevN, meanAdapterContent=fastqcSummaryStats.meanAdapterContent, maxAdapterContent=fastqcSummaryStats.maxAdapterContent)

def addPaftolFastqFiles(fastqFnameList):
    productionDatabaseDetails = getProductionDatabaseDetails()
    connection = productionDatabaseDetails.makeConnection()
    productionDatabase = paftol.database.production.ProductionDatabase(connection)
    connection.close()
    analysisDatabaseDetails = getAnalysisDatabaseDetails()
    connection = analysisDatabaseDetails.makeConnection()
    analysisDatabase = paftol.database.analysis.AnalysisDatabase(connection)
    newPaftolFastqFileList = []
    for fastqFname in fastqFnameList:
        print 'fastqFname: ', fastqFname
        idSequencing, orientation = parseCanonicalSymlink(fastqFname)
        print 'idSequencing: ', idSequencing , 'orientation: ', orientation
        if idSequencing is not None:
            md5sum = paftol.tools.md5HexdigestFromFile(fastqFname)
            fastqFile = findFastqFile(analysisDatabase, fastqFname)
            if fastqFile is None:
                fastqcStats = paftol.tools.generateFastqcStats(fastqFname)
                newFastqStats = fastqStatsFromFastqcStats(fastqcStats)
                ### 1.10.2019 - I think script got up to here
                newFastqFile = paftol.database.analysis.FastqFile(None, filename=fastqFname, md5sum=md5sum, fastqStats=newFastqStats)
                
                ### Paul B. corrected from: newPaftolFastqFile = paftol.database.analysis.PaftolFastqFile(None, idSequencing, newFastqFile)
                newPaftolFastqFile = paftol.database.analysis.PaftolFastqFile(None, idSequencing=idSequencing, fastqFile=newFastqFile)
                ###print 'Filename in newPaftolFastqFile: ', newPaftolFastqFile.newFastqFile.filename  - NB - why does this not work?
                newPaftolFastqFileList.append(newPaftolFastqFile)
            else:
                if fastqFile.md5sum == md5sum:
                    logger.info('fastq file %s already in database, verified md5sum', fastqFname)
                else:
                    raise StandardError, 'fastq file %s in database with md5sum = %s, but found md5sum = %s' % (fastqFname, fastqFile.md5sum, md5sum)
        else:
            logger.warning('not a canonical PAFTOL fastq name: %s', fastqFname)
    transactionSuccessful = insertPaftolFastqFileList(connection, newPaftolFastqFileList)
    connection.close()
    return transactionSuccessful

    
def findGeneType(analysisDatabase, geneTypeName):
    for geneType in analysisDatabase.geneTypeDict.values():
        if geneType.geneTypeName == geneTypeName:
            return geneType
    return None


def addTargetsFile(targetsFname, description=None, insertGenes=False, geneTypeName=None):
    if insertGenes and geneTypeName is None:
        raise StandardError, 'illegal state: insertion of new genes requested but no gene type name given'
    md5sum = paftol.tools.md5HexdigestFromFile(targetsFname)        
    paftolTargetSet = paftol.PaftolTargetSet()
    paftolTargetSet.readFasta(targetsFname)
    numSequences = len(paftolTargetSet.getSeqRecordList())
    analysisDatabaseDetails = getAnalysisDatabaseDetails()
    connection = analysisDatabaseDetails.makeConnection()
    analysisDatabase = paftol.database.analysis.AnalysisDatabase(connection)
    geneType = findGeneType(analysisDatabase, geneTypeName)
    if geneType is None:
        connection.close()
        raise StandardError, 'unknown gene type: %s' % geneTypeName
    targetsFile = findFastaFile(analysisDatabase, targetsFname)
    if targetsFile is not None:
        connection.close()
        if targetsFile.md5sum == md5sum:
            logger.info('targets file %s already in database, verified md5sum', targetsFname)
            return
        else:
            raise StandardError, 'targets file %s in database with md5sum = %s, but found md5sum = %s' % (targetsFname, targetsFile.md5sum, md5sum)
    # connection.start_transaction(isolation_level='REPEATABLE READ', readonly=False)
    targetsFastaFile = paftol.database.analysis.FastaFile(None, targetsFname, md5sum, None, numSequences)
    
    knownGeneNameList = [paftolGene.geneName for paftolGene in analysisDatabase.paftolGeneDict.values()]
    missingGeneNameList = []
    for geneName in paftolTargetSet.paftolGeneDict.keys():
        if geneName not in knownGeneNameList:
            missingGeneNameList.append(geneName)
    if not insertGenes and len(missingGeneNameList) > 0:
        connection.close()
        raise StandardError, 'missing genes: %s' % ', '.join(missingGeneNameList)
    newPaftolGeneList = []
    for missingGeneName in missingGeneNameList:
        newPaftolGene = paftol.database.analysis.PaftolGene(None, missingGeneName, geneType)
        newPaftolGeneList.append(newPaftolGene)
    paftolGeneDict = {}
    for newPaftolGene in newPaftolGeneList:
        paftolGeneDict[newPaftolGene.geneName] = newPaftolGene
    for paftolGene in analysisDatabase.paftolGeneDict.values():
        paftolGeneDict[paftolGene.geneName] = paftolGene
    referenceTargetList = []
    for paftolGene in paftolTargetSet.paftolGeneDict.values():
        for paftolTarget in paftolGene.paftolTargetDict.values():
            referenceTargetList.append(paftol.database.analysis.ReferenceTarget(None, paftolGeneDict[paftolGene.name], paftolTarget.organism.name, len(paftolTarget.seqRecord), targetsFastaFile))
    transactionSuccessful = False
    lockCursor = connection.cursor()
    lockCursor.execute('LOCK TABLE FastaFile WRITE, FastqFile WRITE, FastqStats WRITE, GeneType WRITE, PaftolGene WRITE, ReferenceTarget WRITE')
    try:
        logger.info('adding new targets file %s', targetsFname)
        cursor = connection.cursor(prepared=True)
        try:
            for newPaftolGene in newPaftolGeneList:
                newPaftolGene.id = generateUnusedPrimaryKey(cursor, 'PaftolGene')
                newPaftolGene.insertIntoDatabase(cursor)
            targetsFastaFile.id = generateUnusedPrimaryKey(cursor, 'FastaFile')
            targetsFastaFile.insertIntoDatabase(cursor)
            for referenceTarget in referenceTargetList:
                referenceTarget.id = generateUnusedPrimaryKey(cursor, 'ReferenceTarget')
                referenceTarget.insertIntoDatabase(cursor)
            connection.commit()
            transactionSuccessful = True
        finally:
            if not transactionSuccessful:
                connection.rollback()
            cursor.close()
    finally:
        if not transactionSuccessful:
            connection.rollback()
        lockCursor.execute('UNLOCK TABLES')
        lockCursor.close()
    connection.close()
    return transactionSuccessful


def findFastqFiles(analysisDatabase, result):
    fwdFastqFname = os.path.basename(result.forwardFastq)
    revFastqFname = os.path.basename(result.reverseFastq)
    fwdFastqFile = None
    revFastqFile = None
    for fastqFile in analysisDatabase.fastqFileDict.values():
        if fastqFile.filename == fwdFastqFname:
            fwdFastqFile = fastqFile
        if fastqFile.filename == revFastqFname:
            revFastqFile = fastqFile
    return fwdFastqFile, revFastqFile


def findContigRecoveryForFastqFname(analysisDatabase, fastqFname):
    fastqFile = findFastqFile(analysisDatabase, fastqFname)
    if fastqFile is None:
        return None
    if len(fastqFile.contigRecoveryFwdFastqList) + len(fastqFile.contigRecoveryRevFastqList) > 1:
        raise StandardError, 'multiple ContigRecovery instances for %s: %s' % (fastqFname, ', '.join(['%d' % cr.id for cr in fastqFile.contigRecoveryFwdFastqList +  fastqFile.contigRecoveryRevFastqList]))
    if len(fastqFile.contigRecoveryFwdFastqList) == 1:
        return fastqFile.contigRecoveryFwdFastqList[0]
    if len(fastqFile.contigRecoveryRevFastqList) == 1:
        return fastqFile.contigRecoveryRevFastqList[0]
    return None


def preRecoveryCheck(forwardFastqFname, reverseFastqFname):
    msgList = []
    analysisDatabase = getAnalysisDatabase()
    contigRecovery = findContigRecoveryForFastqFname(analysisDatabase, forwardFastqFname)
    if contigRecovery is not None:
        msgList.append('recovery already done for %s (contigRecovery.id = %d)' % (forwardFastqFname, contigRecovery.id))
    contigRecovery = findContigRecoveryForFastqFname(analysisDatabase, reverseFastqFname)
    if contigRecovery is not None:
        msgList.append('recovery already done for %s (contigRecovery.id = %d)' % (reverseFastqFname, contigRecovery.id))
    if len(msgList) > 0:
        raise StandardError, ', '.join(msgList)


def findContigRecoveryForSequencing(analysisDatabase, idSequencing):
    fastqFileList = []
    for paftolFastqFile in analysisDatabase.paftolFastqFileDict.values():
        if paftolFastqFile.idSequencing is not None and paftolFastqFile.idSequencing == idSequencing:
            if paftolFastqFile.fastqFile is None:
                raise StandardError, 'illegal state: PaftolFastqFile instance %d has no fastqFile' % paftolFastqFile.id
            fastqFile.List.append(paftolFastqFile.fastqFile)
    contigRecoveryList = []
    for fastqFile in fastqFileList:
        for contigRecovery in fastqFile.contigRecoveryFwdFastqList:
            if contigRecovery not in contigRecoveryList:
                contigRecoveryList.append(contigRecovery)
        for contigRecovery in fastqFile.contigRecoveryRevFastqList:
            if contigRecovery not in contigRecoveryList:
                contigRecoveryList.append(contigRecovery)
    if len(contigRecoveryList) == 0:
        return None
    elif len(contigRecoveryList) == 1:
        return contigRecoveryList[0]
    else:
        raise StandardError, 'idSequencing %d: found multiple ContigRecovery instances: %s' % (idSequencing, ', '.join(['%d' % cr.id for cr in contigRecoveryList]))

    
def findReferenceTarget(analysisDatabase, geneName, paftolOrganism):
    logger.debug('searching for %s-%s', paftolOrganism, geneName)
    for referenceTarget in analysisDatabase.referenceTargetDict.values():
        logger.debug('checking %s-%s', referenceTarget.paftolOrganism, referenceTarget.paftolGene.geneName)
        if referenceTarget.paftolOrganism == paftolOrganism and referenceTarget.paftolGene.geneName == geneName:
            return referenceTarget
    return None


def addRecoveryResult(result):
    analysisDatabaseDetails = getAnalysisDatabaseDetails()
    connection = analysisDatabaseDetails.makeConnection()
    analysisDatabase = paftol.database.analysis.AnalysisDatabase(connection)
    targetsFastaFile = findFastaFile(analysisDatabase, result.paftolTargetSet.fastaHandleStr)
    numMappedReads = len(result.paftolTargetSet.getMappedReadNameSet())
    numUnmappedReads = result.paftolTargetSet.numOfftargetReads
    if targetsFastaFile is None:
        # raise StandardError, 'targets file "%s" not in database' % result.paftolTargetSet.fastaHandleStr
        logger.info('unknown targets file "%s" -- continuing')
    fwdFastqFile, revFastqFile = findFastqFiles(analysisDatabase, result)
    if fwdFastqFile is None:
        raise StandardError, 'forward fastq file "%s" not in database' % result.forwardFastq
    if revFastqFile is None:
        raise StandardError, 'reverse fastq file "%s" not in database' % result.reverseFastq
    trimmedForwardFastqStats = None
    if result.forwardTrimmedPairedFastqcStats is not None:
        trimmedForwardFastqStats = fastqStatsFromFastqcStats(result.forwardTrimmedPairedFastqcStats)
    trimmedReverseFastqStats = None
    if result.reverseTrimmedPairedFastqcStats is not None:
        trimmedReverseFastqStats = fastqStatsFromFastqcStats(result.reverseTrimmedPairedFastqcStats)
    paftolGeneDict = {}
    for paftolGene in analysisDatabase.paftolGeneDict.values():
        paftolGeneDict[paftolGene.geneName] = paftolGene
    for geneName in result.contigDict:
        if geneName not in paftolGeneDict:
            raise StandardError, 'found gene %s in result but it is not in the analysis database' % geneName
    contigFastaFile = None
    if result.contigFastaFname is not None:
        contigFastaFile = paftol.database.analysis.FastaFile(None, result.contigFastaFname, paftol.tools.md5HexdigestFromFile(result.contigFastaFname), None, len(paftol.tools.fastaSeqRecordList(result.contigFastaFname)))
    contigRecovery = paftol.database.analysis.ContigRecovery(id=None, fwdFastq=fwdFastqFile, revFastq=revFastqFile, fwdTrimmedFastqStats=trimmedForwardFastqStats, revTrimmedFastqStats=trimmedReverseFastqStats, contigFastaFile=contigFastaFile, targetsFastaFile=targetsFastaFile, numMappedReads=numMappedReads, numUnmappedReads=numUnmappedReads, softwareVersion=paftol.__version__, cmdLine=result.cmdLine)
    recoveredContigList = []
    for geneName in result.contigDict:
        if result.contigDict[geneName] is not None and len(result.contigDict[geneName]) > 0:
            for contig in result.contigDict[geneName]:
                representativeReferenceTarget = findReferenceTarget(analysisDatabase, geneName, result.representativePaftolTargetDict[geneName].organism.name)
                if representativeReferenceTarget is None:
                    raise StandardError, 'unknown reference target for geneName = %s, organismName = %s' % (geneName, result.representativePaftolTargetDict[geneName].organism.name)
                recoveredContig = paftol.database.analysis.RecoveredContig(None, contigRecovery, paftolGeneDict[geneName], len(contig), None, None, representativeReferenceTarget)
                recoveredContigList.append(recoveredContig)
    transactionSuccessful = False
    lockCursor = connection.cursor(prepared=False)
    lockCursor.execute('LOCK TABLE FastaFile WRITE, FastqFile WRITE, FastqStats WRITE, ContigRecovery WRITE, RecoveredContig WRITE')
    try:
        cursor = connection.cursor(prepared=True)
        try:
            if trimmedForwardFastqStats is not None:
                trimmedForwardFastqStats.id = generateUnusedPrimaryKey(cursor, 'FastqStats')
                trimmedForwardFastqStats.insertIntoDatabase(cursor) 
            if trimmedReverseFastqStats is not None:
                trimmedReverseFastqStats.id = generateUnusedPrimaryKey(cursor, 'FastqStats')
                trimmedReverseFastqStats.insertIntoDatabase(cursor) 
            if contigFastaFile is not None:
                contigFastaFile.id = generateUnusedPrimaryKey(cursor, 'FastaFile')
                contigFastaFile.insertIntoDatabase(cursor)
            contigRecovery.id = generateUnusedPrimaryKey(cursor, 'ContigRecovery')
            contigRecovery.insertIntoDatabase(cursor)
            for recoveredContig in recoveredContigList:
                recoveredContig.id = generateUnusedPrimaryKey(cursor, 'RecoveredContig')
                recoveredContig.insertIntoDatabase(cursor)
            connection.commit()
            transactionSuccessful = True
        finally:
            if not transactionSuccessful:
                connection.rollback()
            cursor.close()
    finally:
        if not transactionSuccessful:
            connection.rollback()
        lockCursor.execute('UNLOCK TABLES')
        lockCursor.close()
    connection.close()
    return transactionSuccessful
