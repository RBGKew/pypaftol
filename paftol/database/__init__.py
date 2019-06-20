import sys
import re
import os
import os.path

import mysql.connector

import paftol
import paftol.tools


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
    productionDatabaseDetails = None
    with open(detailsFname, 'r') as f:
        productionDatabaseDetails = PaftolDatabaseDetails(f)
    return productionDatabaseDetails


def getProductionDatabaseDetails(detailsFname=None):
    if detailsFname is None:
        detailsFname = os.path.join(os.environ['HOME'], '.paftol', 'productiondb.cfg')
    return getDatabaseDetails(detailsFname)


def getAnalysisDatabaseDetails(detailsFname=None):
    if detailsFname is None:
        detailsFname = os.path.join(os.environ['HOME'], '.paftol', 'analysisdb.cfg')
    return getDatabaseDetails(detailsFname)


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
    return 'XPAFTOL_%06d_R%1d.fastq%s' % (sequence.idSequencing, orientation, gzExt)


def parseCanonicalSymlink(symlinkName):
    symlinkRe = re.compile('XPAFTOL_([0-9]+)_R([12])\\.fastq')
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


def addTargetsFile(targetsFname, description=None, insertGenes=False, geneTypeName=None):
    if insertGenes and geneTypeName is None:
        raise StandardError, 'illegal state: insertion of new genes requested but no gene type name given'
    md5sum = paftol.tools.md5HexdigestFromFile(targetsFname)
    paftolTargetSet = paftol.PaftolTargetSet()
    paftolTargetSet.readFasta(targetsFname)
    numSequences = len(paftolTargetSet.getSeqRecordList())
    analysisDatabaseDetails = getAnalysisDatabaseDetails()
    connection = analysisDatabaseDetails.makeConnection()
    # connection.start_transaction(isolation_level='REPEATABLE READ', readonly=False)
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
    sys.stderr.write('%s\n' % sqlStatement)
    sys.stderr.write('%s\n' % str(valueTuple))
    cursor.execute(sqlStatement, valueTuple)
    for paftolGene in paftolTargetSet.paftolGeneDict.values():
        for paftolTarget in paftolGene.paftolTargetDict.values():
            newReferenceTargetId = generateUnusedPrimaryKey(connection, 'ReferenceTarget')
            cursor.execute('INSERT INTO ReferenceTarget (id, paftolGeneId, paftolOrganism, paftolTargetLength, targetsFastaFileId) VALUES (%s, %s, %s, %s, %s)', (newReferenceTargetId, geneIdDict[paftolTarget.paftolGene.name], paftolTarget.organism.name, len(paftolTarget.seqRecord), newFastaFileId, ))
    cursor.close()
    connection.commit()
    connection.close()
