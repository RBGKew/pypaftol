import sys
import re
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
    productionDatabase = paftol.database.production.AnalysisDatabase(connection)
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


def findSequenceForFastqFname(productionDatabase, fastqFname):
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
    symlinkRe = re.compile('PAFTOL_([0-9]+)_R([12])\\.fastq')
    m = symlinkRe.match(symlinkName)
    if m is not None:
        return int(m.group(1)), int(m.group(2))
    return None, None


def makeSymlink(symlinkDirname, sequence, fastqFname):
    orientation = paftol.tools.fastqOrientation(fastqFname)
    gzipped = paftol.tools.isGzipped(fastqFname)
    symlinkName = canonicalSymlinkName(sequence, orientation, gzipped)
    symlinkPath = os.path.join(symlinkDirname, symlinkName)
    if os.path.lexists(symlinkPath) or os.path.exists(symlinkPath):
        logger.warning('sequence %d: link %s already exists', sequence.idSequencing, symlinkPath)
    else:
        os.symlink(fastqFname, symlinkPath)


def generateUnusedPrimaryKey(connection, tableName, primaryKeyColumnName='id'):
    sqlStatement = 'SELECT max(%s) FROM %s' % (primaryKeyColumnName, tableName)
    cursor = connection.cursor()
    cursor.execute(sqlStatement)
    row = cursor.fetchone()
    maxPk = 0
    if row is not None and row[0] is not None:
        maxPk = int(row[0])
    cursor.close()
    return maxPk + 1


def insertGene(connection, geneName, geneTypeId):
    cursor = connection.cursor(prepared=True)
    paftolGeneId = generateUnusedPrimaryKey(connection, 'PaftolGene')
    cursor.execute('INSERT INTO PaftolGene (id, geneName, geneTypeId) VALUES (%s, %s, %s)', (paftolGeneId, geneName, geneTypeId, ))
    cursor.close()
    return paftolGeneId

    
def insertFastaFile(connection, fastaFname, dirname=None):
    fastaPath = fastaFname
    if dirname is not None:
        fastaPath = os.path.join(dirname, fastaFname)
    md5 = paftol.tools.md5HexdigestFromFile(fastaPath)
    numSequences = len(paftol.tools.fastaSeqRecordList(fastaPath))
    cursor = connection.cursor(prepared=True)
    fastaFileId = generateUnusedPrimaryKey(connection, 'PaftolGene')
    cursor.execute('INSERT INTO FastaFile (id, filename, md5sum, numSequences) VALUES (%s, %s, %s, %s)', (fastaFileId, fastaFname, md5, numSequences, ))
    cursor.close()
    return fastaFileId


def addFastqStats(connection, fastqcStats):
    fastqcSummaryStats = paftol.tools.FastqcSummaryStats(fastqcStats)
    cursor = connection.cursor(prepared=True)
    fastqStatsId = generateUnusedPrimaryKey(connection, 'FastqStats')
    cursor.execute('INSERT INTO FastqStats (id, numReads, qual28, meanA, meanC, meanG, meanT, stddevA, stddevC, stddevG, stddevT, meanN, stddevN, meanAdapterContent, maxAdapterContent) VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s)', (fastqStatsId, fastqcSummaryStats.numReads, fastqcSummaryStats.qual28, fastqcSummaryStats.meanA, fastqcSummaryStats.meanC, fastqcSummaryStats.meanG, fastqcSummaryStats.meanT, fastqcSummaryStats.stddevA, fastqcSummaryStats.stddevC, fastqcSummaryStats.stddevG, fastqcSummaryStats.stddevT, fastqcSummaryStats.meanN, fastqcSummaryStats.stddevN, None, None))
    cursor.close()
    return fastqStatsId


def addPaftolFastqFiles(fastqFnameList):
    productionDatabaseDetails = getProductionDatabaseDetails()
    connection = productionDatabaseDetails.makeConnection()
    productionDatabase = paftol.database.production.ProductionDatabase(connection)
    connection.close()
    analysisDatabaseDetails = getAnalysisDatabaseDetails()
    connection = analysisDatabaseDetails.makeConnection()
    analysisDatabase = paftol.database.analysis.AnalysisDatabase(connection)
    cursor = connection.cursor(prepared=True)
    for fastqFname in fastqFnameList:
        idSequencing, orientation = parseCanonicalSymlink(fastqFname)
        if idSequencing is not None:
            md5sum = paftol.tools.md5HexdigestFromFile(fastqFname)
            fastqFile = findFastqFile(analysisDatabase, fastqFname)
            if fastqFile is None:
                fastqcStats = paftol.tools.generateFastqcStats(fastqFname)
                fastqStatsId = addFastqStats(connection, fastqcStats)
                fastqFileId = generateUnusedPrimaryKey(connection, 'FastqFile')
                paftolFastqFileId = generateUnusedPrimaryKey(connection, 'PaftolFastqFile')
                cursor.execute('INSERT INTO FastqFile (id, filename, md5sum, enaAccession, description, fastqStatsId) VALUES (%s, %s, %s, %s, %s, %s)', (fastqFileId, fastqFname, md5sum, None, None, fastqStatsId, ))
                cursor.execute('INSERT INTO PaftolFastqFile (id, idSequencing, fastqFileId) VALUES (%s, %s, %s)', (paftolFastqFileId, idSequencing, fastqFileId, ))
                logger.info('added new fastq file %s', fastqFname)
            else:
                if fastqFile.md5sum == md5sum:
                    logger.info('fastq file %s already in database, verified md5sum', fastqFname)
                else:
                    raise StandardError, 'fastq file %s in database with md5sum = %s, but found md5sum = %s' % (fastqFname, fastqFile.md5sum, md5sum)
        else:
            logger.warning('not a canonical PAFTOL fastq name: %s', fastqFname)
    connection.commit()
    connection.close()


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
    # connection.start_transaction(isolation_level='REPEATABLE READ', readonly=False)
    targetsFile = findFastaFile(analysisDatabase, targetsFname)
    if targetsFile is not None:
        if targetsFile.md5sum == md5sum:
            logger.info('targets file %s already in database, verified md5sum', targetsFname)
            connection.close()
            return
        else:
            raise StandardError, 'targets file %s in database with md5sum = %s, but found md5sum = %s' % (targetsFname, targetsFile.md5sum, md5sum)
    logger.info('adding new targets file %s', targetsFname)
    newFastaFileId = generateUnusedPrimaryKey(connection, 'FastaFile')
    cursor = connection.cursor(prepared=True)
    geneTypeId = None
    if geneTypeName is not None:
        cursor.execute('SELECT id FROM GeneType WHERE geneTypeName = %s', (geneTypeName, ))
        row = cursor.fetchone()
        if row is None:
            cursor.close()
            connection.close()
            raise StandardError, 'unknown gene type: %s' % geneTypeName
        if row[0] is not None:
            geneTypeId = int(row[0])
    geneNameList = paftolTargetSet.paftolGeneDict.keys()
    missingGeneNameList = []
    for geneName in geneNameList:
        cursor.execute('SELECT id, geneName FROM `PaftolGene` WHERE geneName = %s', (geneName, ))
        row = cursor.fetchone()
        if row is None:
            missingGeneNameList.append(geneName)
    if not insertGenes and len(missingGeneNameList) > 0:
        cursor.close()
        connection.close()
        raise StandardError, 'missing genes: %s' % ', '.join(missingGeneNameList)
    # sys.stderr.write('%s\n' % str(missingGeneNameList))
    for geneName in missingGeneNameList:
        insertGene(connection, geneName, geneTypeId)
    geneIdDict = {}
    cursor.execute('SELECT id, geneName FROM PaftolGene')
    for row in cursor:
        geneId = int(row[0])
        geneName = str(row[1])
        if geneName in geneNameList:
            geneIdDict[geneName] = geneId
    valueTuple = (newFastaFileId, targetsFname, md5sum, description, numSequences, )
    sqlStatement = 'INSERT INTO FastaFile (id, filename, md5sum, description, numSequences) VALUES (%s, %s, %s, %s, %s)'
    # sys.stderr.write('%s\n' % sqlStatement)
    # sys.stderr.write('%s\n' % str(valueTuple))
    cursor.execute(sqlStatement, valueTuple)
    for paftolGene in paftolTargetSet.paftolGeneDict.values():
        for paftolTarget in paftolGene.paftolTargetDict.values():
            newReferenceTargetId = generateUnusedPrimaryKey(connection, 'ReferenceTarget')
            cursor.execute('INSERT INTO ReferenceTarget (id, paftolGeneId, paftolOrganism, paftolTargetLength, targetsFastaFileId) VALUES (%s, %s, %s, %s, %s)', (newReferenceTargetId, geneIdDict[paftolTarget.paftolGene.name], paftolTarget.organism.name, len(paftolTarget.seqRecord), newFastaFileId, ))
    cursor.close()
    connection.commit()
    connection.close()


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


def addRecoveryResult(result):
    analysisDatabaseDetails = getAnalysisDatabaseDetails()
    connection = analysisDatabaseDetails.makeConnection()
    analysisDatabase = paftol.database.analysis.AnalysisDatabase(connection)
    targetsFastaFile = findFastaFile(analysisDatabase, result.paftolTargetSet.fastaHandleStr)
    numMappedReads = len(result.paftolTargetSet.getMappedReadNameSet())
    numUnmappedReads = result.paftolTargetSet.numOfftargetReads
    targetsFastaFileId = None
    if targetsFastaFile is None:
        # raise StandardError, 'targets file "%s" not in database' % result.paftolTargetSet.fastaHandleStr
        pass
    else:
        targetsFastaFileId = targetsFastaFile.id
    fwdFastqFile, revFastqFile = findFastqFiles(analysisDatabase, result)
    if fwdFastqFile is None:
        raise StandardError, 'forward fastq file "%s" not in database' % result.forwardFastq
    if revFastqFile is None:
        raise StandardError, 'reverse fastq file "%s" not in database' % result.reverseFastq
    trimmedForwardFastqStatsId = None
    if result.forwardTrimmedPairedFastqcStats is not None:
        trimmedForwardFastqStatsId = addFastqStats(connection, result.forwardTrimmedPairedFastqcStats)
    trimmedReverseFastqStatsId = None
    if result.reverseTrimmedPairedFastqcStats is not None:
        trimmedReverseFastqStatsId = addFastqStats(connection, result.reverseTrimmedPairedFastqcStats)
    paftolGeneEntityDict = {}
    for paftolGeneEntity in analysisDatabase.paftolGeneDict.values():
        paftolGeneEntityDict[paftolGeneEntity.geneName] = paftolGeneEntity
    for geneName in result.contigDict:
        if geneName not in paftolGeneEntityDict:
            raise StandardError, 'found gene %s in result but it is not in the analysis database' % geneName
    contigFastaFileId = None
    if result.contigFastaFname is not None:
        contigFastaFileId = insertFastaFile(connection, result.contigFastaFname)
    contigRecoveryId = generateUnusedPrimaryKey(connection, 'ContigRecovery')
    cursor = connection.cursor(prepared=True)
    cursor.execute('INSERT INTO ContigRecovery (id, fwdFastqId, revFastqId, fwdTrimmedFastqStatsId, revTrimmedFastqStatsId, contigFastaFileId, targetsFastaFileId, numMappedReads, numUnmappedReads, softwareVersion, cmdLine) VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s)', (contigRecoveryId, fwdFastqFile.id, revFastqFile.id, trimmedForwardFastqStatsId, trimmedReverseFastqStatsId, contigFastaFileId, targetsFastaFileId, numMappedReads, numUnmappedReads, paftol.__version__, result.cmdLine))
    # FIXME: should check result
    for geneName in result.contigDict:
        if result.contigDict is not None and len(result.contigDict[geneName]) > 0:
            for contig in result.contigDict[geneName]:
                recoveredContigId = generateUnusedPrimaryKey(connection, 'RecoveredContig')
                cursor.execute('INSERT INTO RecoveredContig (id, contigRecoveryId, paftolGeneId, seqLength, fwdPaftolFastqId, revPaftolFastqId, representativeReferenceTargetId) VALUES (%s, %s, %s, %s, %s, %s, %s)', (recoveredContigId, None, paftolGeneEntityDict[geneName].id, len(contig), None, None, None))
    connection.commit()
    connection.close()
