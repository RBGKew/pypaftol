# Copyright (c) 2020 The Board of Trustees of the Royal Botanic Gardens, Kew

import sys
import re
import copy
import os
import os.path
import logging
import unicodedata
import datetime
import time             ### Paul B. added to sleep after commiting
import random           # Paul B. added
import hashlib          ### Paul B. added - md5 module is deprecated       

import mysql.connector
import Bio.SeqIO

import paftol
import paftol.tools
import paftol.database.analysis
import paftol.database.production


logger = logging.getLogger(__name__)


def strOrNone(x):
    if x is None:
        return None
    # FIXME: coercing datetime to strings, ORM generator should really handle these separately
    elif isinstance(x, datetime.datetime):
        return str(x)
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
    ''' Finds the targets file from the ReferenceTarget table.

    
    Returns a ReferenceTarget row object. 
    '''
    # Paul B. - now getting file from the ReferenceTarget table
    #for fastaFile in analysisDatabase.fastaFileDict.values():   # fastaFile is a row object containing column headers as values
    for fastaFile in analysisDatabase.referenceTargetDict.values(): 
        # Paul B. filename variable has changed:
        #if fastaFile.filename == fastaFname:
        if fastaFile.targetsFastaFile == fastaFname:
            return fastaFile
    return None


def findFastqFile(analysisDatabase, fastqFname):
    # Paul B. - now using inputSequenceDict:
    # NB - each value of the dict is a row object containing column headers as values
    # which gets returned, either from FastqFile (old db) or InputSequence (new db).
    #for fastqFile in analysisDatabase.fastqFileDict.values():
    for inputSequence in analysisDatabase.inputSequenceDict.values():
        if inputSequence.filename == fastqFname:
            return inputSequence
    return None


class PaftolDatabaseDetails(object):

    reDbusername = re.compile('user: *([^ ]+)')     # Paul B - Renamed to 'user' to fit with the actual keyword required
    reDbpassword = re.compile('password: *([^ ]+)') 
    reDbhost = re.compile('host: *([^ ]+)')
    reDbname = re.compile('database: *([^ ]+)')     # Paul B - Renamed to 'database' to fit with the actual keyword required
    reDbport = re.compile('port: *([^ ]+)')         # Paul B - added provision to name a mysql port to use

    def __init__(self, detailsFile=None):
        if detailsFile is None:
            self.dbusername = None
            self.dbpassword = None
            self.dbhost = None
            self.dbname = None
            self.dbport = None  # Paul B. added
        else:
            self.readFile(detailsFile)

    def readDetailsLine(self, detailsFile, detailsRe, errorMsg):
        line = detailsFile.readline()
        # sys.stderr.write('got line: "%s"\n' % repr(line))
        m = detailsRe.match(line.strip())
        if m is None:
            ### Need to compare line with ^port:''
            rePortFlag = re.compile('port:')        
            if rePortFlag.match(line.strip()) is None:  # Paul B. added - optional parameter - but what happens if this parameter is empty but shouldn't be?! 
                return None                             # Paul B. added
            else:                                       # Paul B. added
                raise StandardError, errorMsg
        return m.group(1)

    def readFile(self, detailsFile):
        self.dbusername = self.readDetailsLine(detailsFile, self.reDbusername, 'malformed dbusername line')
        self.dbpassword = self.readDetailsLine(detailsFile, self.reDbpassword, 'malformed dbpassword line')
        self.dbhost = self.readDetailsLine(detailsFile, self.reDbhost, 'malformed dbhost line')
        self.dbname = self.readDetailsLine(detailsFile, self.reDbname, 'malformed dbname line')
        self.dbport = self.readDetailsLine(detailsFile, self.reDbport, 'malformed dbport line') # Paul B. added - NB - the order these lines needs to match that of the db config file
                                                                                                # NB - the db test code works but lacks the test for the port    

    def makeConnection(self):
        if self.dbport is None:    # Paul B. added
            return mysql.connector.connection.MySQLConnection(user=self.dbusername, password=self.dbpassword, host=self.dbhost, database=self.dbname)
        else: 
            return mysql.connector.connection.MySQLConnection(user=self.dbusername, password=self.dbpassword, host=self.dbhost, database=self.dbname, port=self.dbport) # Paul B. added


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
    #paftolPrefixedFastqFnameRe = re.compile('PAFTOL[-_]([0-9]+)_R[12]_[0-9]+\\.fastq')    # Paul B. changed to also match e.g. PAFTOL_007767_1.fastq.gz
    paftolPrefixedFastqFnameRe = re.compile('PAFTOL[-_]([0-9]+)_R?[12](_[0-9]+)?\\.fastq')  

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
    ''' Paul B. added:
        Finds sequence for fastq filenme using two strategies:
        1. Tries to extract what should be the idSequencing id just after the PAFTOL prefix
           and use that id to look up the correct fastq file.
           NB - some of these ids are idPaftol ids - these will fail
           Many post pilot fastq files have this format, mostly later sequencing runs. 
        2. If paftolPrefixedNumber is None, the real fastq filename and its path will be
           used to obtain the sequencing pool. This pool id will be used to create a list of
           all fastq filenames from that pool. If one is found, that's OK and a symlink will be made.'
    '''

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
    ''' Paul B. - 25.5.2020
        Doesn't look like this method is used any more - insert to db done via the analysis.py API
    '''
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
    ''' Paul B. added:
        The FastaFile table doesn't exist anymore.
    '''

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


# Paul B. - made a copy of this function just below to use for the new table structure.
def preInsertCheckPaftolFastqFile(paftolFastqFile):
    if paftolFastqFile.id is not None:
        raise StandardError, 'illegal state: PaftolFastqFile instance has id %d, not None' % paftolFastqFile.id
    if paftolFastqFile.fastqFile is None:
        raise StandardError, 'illegal state: PaftolFastqFile has fastqFile attribute set to None'
    if paftolFastqFile.fastqFile.fastqStats is not None and  paftolFastqFile.fastqFile.fastqStats.id is not None:
        raise StandardError, 'illegal state: new FastqFile %s has existing FastqStats %d' % (paftolFastqFile.fastqFile.filename, paftolFastqFile.fastqFile.fastqStats.id)


def preInsertCheckInputSequence(inputSequence):
    ''' Paul B - new method for new db schema - not actually sure why we need to do these checks.
        Why no other checks on other variables and objects?
        Actually does catch this situation: e.g. 5853_R1.fastq if --dataOrigin == PAFTOL

        Note that the names e.g. inputSequence.sraRunSequence don't always exactly match the table names. 
        It's because these are the names of the objects defined in the InputSequence class object. Yes, but
        these names can only come from the db then written out by the API!!
    '''

    if inputSequence.id is not None:
        raise StandardError, 'illegal state: InputSequence instance has id %d, not None' % inputSequence.id
    # Paul B. - added and conditionals - one of these foreign keys needs to be defined: 
    #if inputSequence.paftolSequence is None
    if inputSequence.paftolSequence is None and inputSequence.OneKP_Sequence is None and inputSequence.sraRunSequence is None and inputSequence.annotatedGenome is None and inputSequence.GAP_Sequence is None and inputSequence.UnannotatedGenome is None:
        #raise StandardError, 'illegal state: inputSequence has paftolSequence attribute set to None'
        raise StandardError, 'illegal state: inputSequence has to have a value for one of these attributes: paftolSequence, OneKP_Sequence, sraRunSequence, annotatedGenome, GAP_Sequence, UnannotatedGenome'
    if inputSequence.fastqStats is not None and  inputSequence.fastqStats.id is not None:
        raise StandardError, 'illegal state: new inputSequence %s has existing FastqStats %d' % (inputSequence.filename, inputSequence.fastqStats.id)


### Paul B. added this method - redundant with addPaftolFastqFiles() above - can delete
def insertPaftolFastqFileList(connection, paftolFastqFileList):
    for paftolFastqFile in paftolFastqFileList:
        preInsertCheckPaftolFastqFile(paftolFastqFile)
    insertedPaftolFastqFileList = []        ### Paul B. - this never seems to be used
    transactionSuccessful = False
    # Paul B. - removed table locking and introduced auto_increment for each primary key:
    #lockCursor = connection.cursor()
    #lockCursor.execute('LOCK TABLE PaftolFastqFile WRITE, FastqFile WRITE, FastqStats WRITE')
    try:
        cursor = connection.cursor(prepared=True)
        try:
            for paftolFastqFile in paftolFastqFileList:
                insertedPaftolFastqFile = copy.deepcopy(paftolFastqFile)
                # Paul B - altered for auto_increment:
                #insertedPaftolFastqFile.id = generateUnusedPrimaryKey(cursor, 'PaftolFastqFile')
                #insertedPaftolFastqFile.fastqFile.id = generateUnusedPrimaryKey(cursor, 'FastqFile')
                if insertedPaftolFastqFile.fastqFile.fastqStats is not None:
                    # Paul B. removed - now using auto_increment
                    #insertedPaftolFastqFile.fastqFile.fastqStats.id = generateUnusedPrimaryKey(cursor, 'FastqStats')
                    insertedPaftolFastqFile.fastqFile.fastqStats.insertIntoDatabase(cursor)
                    insertedPaftolFastqFile.fastqFile.fastqStats.id = cursor.lastrowid
                insertedPaftolFastqFile.fastqFile.insertIntoDatabase(cursor)
                insertedPaftolFastqFile.fastqFile.id = cursor.lastrowid
                insertedPaftolFastqFile.insertIntoDatabase(cursor)
                insertedPaftolFastqFile.id = cursor.lastrowid
            connection.commit()
            transactionSuccessful = True
        finally:
            if not transactionSuccessful:
                connection.rollback()
                print "ERROR: commit unsucessful for insertedPaftolFastqFile.fastqFile.fastqStats.id: ", insertedPaftolFastqFile.fastqFile.fastqStats.id
                print "ERROR: commit unsucessful for insertedPaftolFastqFile.fastqFile.id: ", insertedPaftolFastqFile.fastqFile.id
                print "ERROR: commit unsucessful for insertedPaftolFastqFile.id: ", insertedPaftolFastqFile.id
            cursor.close()
    finally:
        if not transactionSuccessful:
            connection.rollback()
        #lockCursor.execute('UNLOCK TABLES')
        #lockCursor.close()
    return transactionSuccessful


# Paul B. - added a new method equivalent to and to replace insertPaftolFastqFileList():
def insertInputSequenceList(connection, analysisDatabase, inputSequenceList, newPaftolSequence, externalGenesFile):
    for inputSequence in inputSequenceList:
        # Paul B.: preInsertCheckPaftolFastqFile(paftolFastqFile)
        # Now checking entries from the InputSequence table:
        preInsertCheckInputSequence(inputSequence)
    transactionSuccessful = False
    # Paul B. - removed table locking and introduced auto_increment for each primary key:
    #lockCursor = connection.cursor()
    #lockCursor.execute('LOCK TABLE PaftolFastqFile WRITE, FastqFile WRITE, FastqStats WRITE')
    try:
        cursor = connection.cursor(prepared=True)
        # Inserting into paftolSequence once, then id is available for both fastq files.
        # An alternative would be to bring in the newPaftolSequence object to avoid 
        # accessing the PaftolSequence table via an inputSequence object - now doing this.
        #inputSequenceList[0].paftolSequence.insertIntoDatabase(cursor)
        #inputSequenceList[0].paftolSequence.id = cursor.lastrowid
        newPaftolSequence.insertIntoDatabase(cursor)
        newPaftolSequence.id = cursor.lastrowid
        if newPaftolSequence.id is not None:
            print "newPaftolSequence.id:", newPaftolSequence.id
        #print "inputSequenceList[0].paftolSequence.id: ", inputSequenceList[0].paftolSequence.id
        #print "inputSequenceList[1].paftolSequence.id: ", inputSequenceList[1].paftolSequence.id  # proves paftolSequence object is the same instance.
        try:
            for inputSequence in inputSequenceList:
                # Paul B. - don't understand the purpose of the deepcopy - deepcopy is a copy whose items are different instances I think
                #           It now causes an error so will removed it: RuntimeError: maximum recursion depth exceeded while calling a Python object
                #insertedInputSequence = copy.deepcopy(inputSequence)
                #print "insertedInputSequence.fastqStats.numReads: ", insertedInputSequence.fastqStats.numReads
                #print "insertedInputSequence.paftolSequence.idSequencing: ", insertedInputSequence.paftolSequence.idSequencing
                # Paul B - altered for auto_increment:
                #insertedPaftolFastqFile.id = generateUnusedPrimaryKey(cursor, 'PaftolFastqFile')
                #insertedPaftolFastqFile.fastqFile.id = generateUnusedPrimaryKey(cursor, 'FastqFile')
                ###if insertedInputSequence.fastqStats is not None:
                if inputSequence.fastqStats is not None:
                    # Paul B. removed - now using auto_increment
                    #insertedPaftolFastqFile.fastqFile.fastqStats.id = generateUnusedPrimaryKey(cursor, 'FastqStats')
                    ###insertedInputSequence.fastqStats.insertIntoDatabase(cursor)
                    ###insertedInputSequence.fastqStats.id = cursor.lastrowid
                    ###print 'insertedInputSequence.fastqStats.id: '           #, insertedInputSequence.fastqStats.id
                    inputSequence.fastqStats.insertIntoDatabase(cursor)
                    inputSequence.fastqStats.id = cursor.lastrowid
                    if inputSequence.fastqStats.id is not None:
                        print 'inputSequence.fastqStats.id: ', inputSequence.fastqStats.id      # inputSequence.fastqStats.id
                # InputSequence requires the primary keys from FastqStats and PaftolSequence tables (retrieved above) 
                ###insertedInputSequence.insertIntoDatabase(cursor)
                ###insertedInputSequence.id = cursor.lastrowid
                inputSequence.insertIntoDatabase(cursor)
                inputSequence.id = cursor.lastrowid
                if inputSequence.id is not None:
                    print 'insertedInputSequence.id: ', inputSequence.id                    # insertedInputSequence.id
                    print 'insertedInputSequence.filename: ', inputSequence.filename        # insertedInputSequence.filename
                    if inputSequence.dataOrigin.acronym == 'OneKP_Transcripts' or inputSequence.dataOrigin.acronym == 'AG' or inputSequence.dataOrigin.acronym == 'UG':
                        if externalGenesFile is not None:
                            ### NB - 8.9.2020 - should only be one file, so only tolerating a single file here - would be good to break out of loop after
                            addExternalGenes(cursor=cursor, analysisDatabase=analysisDatabase, inputSequence=inputSequence, newDatasetSequence=newPaftolSequence, externalGenesFile=externalGenesFile)
                            ### Consider to break out of loop after proceesing [0] element
            connection.commit()
            transactionSuccessful = True
        finally:
            if not transactionSuccessful:
                connection.rollback()
                print "ERROR: commit unsucessful for insertedInputSequence.fastqStats.id: "             #, insertedInputSequence.fastqStats.id
                print "ERROR: commit unsucessful for insertedInputSequence.id: "                        #, insertedInputSequence.id
                                                                                                        # NB - these variables may not exist if commit fails
                print "ERROR: commit unsucessful for contigRecovery.id, if --addExternalGenes Flag in use"                          #, contigRecovery.id
                print "ERROR: commit unsucessful for recoveredContig.id (for last row created), if --addExternalGenes Flag in use"  #, recoveredContig.id
            cursor.close()
    finally:
        if not transactionSuccessful:
            connection.rollback()
            print "ERROR: commit unsucessful for newPaftolSequence"      #, inputSequenceList[0].paftolSequence.id
        #lockCursor.execute('UNLOCK TABLES')
        #lockCursor.close()
    return transactionSuccessful


def fastqStatsFromFastqcStats(fastqcStats):
    fastqcSummaryStats = paftol.tools.FastqcSummaryStats(fastqcStats)
    # Paul B. - removed None in 1st position for auto_increment:
    return paftol.database.analysis.FastqStats(numReads=fastqcSummaryStats.numReads, qual28=fastqcSummaryStats.qual28, meanA=fastqcSummaryStats.meanA, meanC=fastqcSummaryStats.meanC, meanG=fastqcSummaryStats.meanG, meanT=fastqcSummaryStats.meanT, stddevA=fastqcSummaryStats.stddevA, stddevC=fastqcSummaryStats.stddevC, stddevG=fastqcSummaryStats.stddevG, stddevT=fastqcSummaryStats.stddevT, meanN=fastqcSummaryStats.meanN, stddevN=fastqcSummaryStats.stddevN, meanAdapterContent=fastqcSummaryStats.meanAdapterContent, maxAdapterContent=fastqcSummaryStats.maxAdapterContent)


def rawFilenameStats(filename=None):
    ''' Paul B. added method:
    Gets stats for a raw fasta file - specific to OneKP_Transcripts and AG/UG data only
    
    @param filename: fasta file name
    @type filename: C{str}
    @return: C{numbrSequences}, C{sumLengthOfContigs}

    '''
    numbrSequences = 0
    sumLengthOfContigs = 0
    seqRecords = {}     #

    # The fasta file has to be unzipped for the Bio.SeqIO to work - it does not throw an error
    # if file is zipped and produces gobbledeegook values for numbrSequences, sumLengthOfContigs !
    # Attempting to check if file is zipped here:
    if re.search('.gz$|.bz2$', filename) is not None:
        raise StandardError, 'Raw data fasta file is zipped, needs to be unzipped: %s' % filename

    for record in Bio.SeqIO.parse(filename, "fasta"):     # returns a SeqRecord object, includes a Seq object called seq
        #print record.id, "\n", record.seq, len(record)
        sumLengthOfContigs = sumLengthOfContigs + len(record)
        numbrSequences += 1
        seqRecords[record.id] = record
    logger.info('Raw data file stats  numbrSequences: %s; sumLengthOfContigs: %s ', numbrSequences, sumLengthOfContigs)
    return  numbrSequences, sumLengthOfContigs 


def findSequence(productionDatabase, sampleId):
    ''' Paul B added - Finds a row in the PAFTOL db Sequence table corresponding to the sample Id in the ExternalSequenceID column

    NB - 29.10.2020 - what happens if there are > 1 Sequence entries for the same ExternalSequenceID?
    Then this loop woud have to return a dict of them all then exit if keys > 1 OR pick the one with the largest idSequencing 
    which would be the most recent sequencing. OK for now though.

    NB - 'ExternalSequenceID' db table field name == 'externalSequenceId' in Python Sequence object!
    NB - If the production database is updated with extra rows, unless these rows are required, it seems that the db API doesn't need to be remade.

    Returns a matching sequence table row object if one exists or None
    '''
    for sequence in productionDatabase.sequenceDict.values():   # returns a copy of all dict VALUES i.e. a Sequence row object
        if sequence.externalSequenceId == sampleId:
            return sequence
    return None


def addPaftolFastqFiles(fastqFnameList=None, dataOriginAcronym=None, fastqPath=None, sampleId=None, externalGenesFile=None):    # Paul B. changed to include path to fastq file(s) and sampleId (for use with non-paftol data)
    ''' Adds input sequence file(s) info into the paftol_da database e.g. fastq, fasta files

    Was specific to PAFTOL only data, now can handle more data set types as defined in the DataOrigin table.
    '''
   
    # Code to connect to the production db already existed (21.7.2020) but then never did anything.
    # Now using this connection to look up the idSequencing with the sampleId for other non-paftol data sets, 
    # now recorded in paftol.sequence.ExternalSequenceID
    productionDatabaseDetails = getProductionDatabaseDetails()
    connection = productionDatabaseDetails.makeConnection()
    productionDatabase = paftol.database.production.ProductionDatabase(connection)
    connection.close()
    

    analysisDatabaseDetails = getAnalysisDatabaseDetails()
    connection = analysisDatabaseDetails.makeConnection()
    analysisDatabase = paftol.database.analysis.AnalysisDatabase(connection)
    # Paul B. added to check dataOriginName with name in DataOrigin table:
    # NB - this check means that the --dataOrigin flag is no longer strictly an 'option'
    # but I followed the logic for --geneType flag in addTargetsFile 
    dataOrigin = findDataOrigin(analysisDatabase, dataOriginAcronym)     # Paul B. - returns a DataOrigin object
    if dataOrigin is None:
        connection.close()
        raise StandardError, 'Data origin entry is incorrect. Allowed values are: \'PAFTOL\', \'OneKP_Transcripts\', \'OneKP_Reads\'  \'SRA\', \'GAP\', \'AG\', \'UG\'  '
        ### NB - I think argparse can handle this if required= is set - if so can delete above code
    # Paul B. - need to create a new PaftolSequence object with info from read1 file only.
    # First, get the idSequencing from the R1 file ('paftol' data set only)
    newDatasetSequence = None           # Paul B. added - new table object for any data set type
    idSequencing = None                 # Paul B. added - used later in method so needs to be global (?)
    if dataOriginAcronym == 'PAFTOL':   # Paul B. added
        print 'idSequencing: ', idSequencing , 'orientation: ', orientation
        #newPaftolSequence = paftol.database.analysis.PaftolSequence(idSequencing=idSequencing, replicate=None)
        if idSequencing is not None:    # idSequencing may not be found
            newDatasetSequence = paftol.database.analysis.PaftolSequence(idSequencing=idSequencing, replicate=None)
        else:
            logger.warning('not a canonical PAFTOL fastq name, can\'t obtain idSequencing identifier: %s', fastqFnameList[0])
    else:
        # For non-paftol data set types, obtaining idSequencing from the paftol.Sequence table instead using sampleId from the sampleId input flag:
        ### NB - might be able to handle data set options better by using word matches rather than conditonals
        if sampleId is None:
            raise StandardError, 'sampleId from the --sampleId option is empty.'
        #print 'sampleId before findSequence: ', sampleId
        sequenceRow = findSequence(productionDatabase, sampleId)
        if sequenceRow is None:
            raise StandardError, 'Sequence row object is empty - externalSequenceID not found for %s' % sampleId
            
        idSequencing = sequenceRow.idSequencing # idSequencing is the primary key for table Sequence so it will always exist 
        logger.info('Found idSequencing for sampleId %s: %s', sampleId, idSequencing)   # Visible with this paftools flag:  --loglevel INFO

        if dataOriginAcronym == 'OneKP_Transcripts' or dataOriginAcronym == 'OneKP_Reads':  # sampleId should exist here
            newDatasetSequence = paftol.database.analysis.OneKP_Sequence(sampleId=sampleId, idSequencing=idSequencing)
        elif dataOriginAcronym == 'SRA': 
            newDatasetSequence = paftol.database.analysis.SRA_RunSequence(accessionId=sampleId, idSequencing=idSequencing)
        elif dataOriginAcronym == 'AG':
            newDatasetSequence = paftol.database.analysis.AnnotatedGenome(accessionId=sampleId, idSequencing=idSequencing)
        elif dataOriginAcronym == 'GAP':
            newDatasetSequence = paftol.database.analysis.GAP_Sequence(sampleId=sampleId, idSequencing=idSequencing)
        elif dataOriginAcronym == 'UG':
            newDatasetSequence = paftol.database.analysis.UnannotatedGenome(accessionId=sampleId, idSequencing=idSequencing)    
        else:
            raise StandardError, 'unknown data origin: %s' % dataOriginAcronym

    # Paul B. - removed: newPaftolFastqFileList = []
    newInputSequenceList = []
    for fastqFname in fastqFnameList:
        print 'fastqFname: ', fastqFname
        # Paul B. - moved above outside of loop:
        #idSequencing, orientation = parseCanonicalSymlink(fastqFname)
        #print 'idSequencing: ', idSequencing , 'orientation: ', orientation
        if idSequencing is not None or sampleId is not None:    # Paul B added sampleId
            md5sum = paftol.tools.md5HexdigestFromFile(fastqFname)
            #fastqFile = findFastqFile(analysisDatabase, fastqFname)
            # Paul - now returns an inputSequence object
            inputSequence = findFastqFile(analysisDatabase, fastqFname)
            # Paul B. - changed: if fastqFile is None:
            # NB - if you insert the same fastq file twice i.e. two R1 files by mistake, this will not
            #      be detected here because the first fastq file is not yet uploaded until the insert
            #      method is called.
            if inputSequence is None:

                # Paul B. recreated the full path to the raw fastq files (have to be unzipped for the above commands!):
                if fastqPath is not None:
                    # Only add .gz ending for fastq files (assumed to be zipped for the raw files):
                    if re.search('.fastq$|.fq$', fastqFname) is not None: 
                        fastqPathName = fastqPath + '/' + fastqFname + '.gz'
                    # Now need to extend the match to fasta files - NB - this time endings are specific to OneKP (.bz2) and AG (.gz) data sets:
                    elif re.search('.fasta$|.fa$', fastqFname) is not None and dataOriginAcronym == 'OneKP_Transcripts':
                        fastqPathName = fastqPath + '/' + fastqFname + '.bz2'
                    elif re.search('.fasta$|.fa$', fastqFname) is not None and dataOriginAcronym == 'AG' or dataOriginAcronym == 'UG':
                        fastqPathName = fastqPath + '/' + fastqFname + '.gz'
                    else:
                        fastqPathName = fastqPath + '/' + fastqFname
                    logger.info('fastqPathName: %s', fastqPathName)
                else:
                    fastqPathName = None
                    logger.info('fastqPathName: %s', fastqPathName)

                # Paul B added:
                if dataOriginAcronym == 'PAFTOL' or dataOriginAcronym == 'OneKP_Reads' or dataOriginAcronym == 'SRA' or dataOriginAcronym == 'GAP':  
                    fastqcStats = paftol.tools.generateFastqcStats(fastqFname)
                    newSeqStats = fastqStatsFromFastqcStats(fastqcStats)      # Paul B. NB - this is a database table object
                elif dataOriginAcronym == 'OneKP_Transcripts' or dataOriginAcronym == 'AG' or dataOriginAcronym == 'UG':  
                    newSeqStats = paftol.database.analysis.FastqStats()


                # Paul B - altered to fit with auto_increment + to add the full path to the fastq file:
                #newFastqFile = paftol.database.analysis.FastqFile(filename=fastqFname, pathName=fastqPathName, md5sum=md5sum, fastqStats=newFastqStats)            
                # Paul B. corrected from: newPaftolFastqFile = paftol.database.analysis.PaftolFastqFile(None, idSequencing, newFastqFile)
                # Paul B - also altered for auto_increment:
                #newPaftolFastqFile = paftol.database.analysis.PaftolFastqFile(idSequencing=idSequencing, fastqFile=newFastqFile)
                ###print 'Filename in newPaftolFastqFile: ', newPaftolFastqFile.newFastqFile.filename  - NB - why does this not work?
                #newPaftolFastqFileList.append(newPaftolFastqFile)
                ### NBNB - it may be that  newPaftolFastqFile goes inside newFastqFile rather than the other way round due to which table
                ### has the foreign key - it might be clear from the DDL.
                # Paul B - 18.2.2020
                # Updated to use the new tables of the new db structure:
                #newInputSequence = paftol.database.analysis.InputSequence(dataOrigin=dataOrigin, sequenceType=None, filename=fastqFname, pathName=fastqPathName, md5sum=md5sum, fastqStats=newFastqStats, paftolSequence=newPaftolSequence, sraRunSequence=None, OneKP_Sequence=None, annotatedGenome=None)
                #newInputSequence = paftol.database.analysis.InputSequence(dataOrigin=dataOrigin, filename=fastqFname, pathName=fastqPathName, md5sum=md5sum, fastqStats=newFastqStats, paftolSequence=newPaftolSequence)
                # Paul B - modified above line to be able to input a generic data set table object from above and add newFastqStats but only for data sets with fastq files.
                if dataOriginAcronym == 'PAFTOL':
                    newInputSequence = paftol.database.analysis.InputSequence(dataOrigin=dataOrigin, filename=fastqFname, pathName=fastqPathName, md5sum=md5sum, fastqStats=newSeqStats, paftolSequence=newDatasetSequence)
                elif dataOriginAcronym == 'OneKP_Reads':
                    newInputSequence = paftol.database.analysis.InputSequence(dataOrigin=dataOrigin, filename=fastqFname, pathName=fastqPathName, md5sum=md5sum, fastqStats=newSeqStats, OneKP_Sequence=newDatasetSequence)
                elif dataOriginAcronym == 'SRA':
                    newInputSequence = paftol.database.analysis.InputSequence(dataOrigin=dataOrigin, filename=fastqFname, pathName=fastqPathName, md5sum=md5sum, fastqStats=newSeqStats, sraRunSequence=newDatasetSequence)
                elif dataOriginAcronym == 'GAP':
                    newInputSequence = paftol.database.analysis.InputSequence(dataOrigin=dataOrigin, filename=fastqFname, pathName=fastqPathName, md5sum=md5sum, fastqStats=newSeqStats, GAP_Sequence=newDatasetSequence)
                elif dataOriginAcronym == 'OneKP_Transcripts':
                    # Also calculate number of genes and sum length of contigs in the raw fasta file and add to the newSeqsStats object:
                    numbrSequences, sumLengthOfContigs = rawFilenameStats(filename=fastqFname)
                    newSeqStats.numbrRecords = numbrSequences
                    newSeqStats.sumLengthOfSeqs = sumLengthOfContigs
                    newInputSequence = paftol.database.analysis.InputSequence(dataOrigin=dataOrigin, filename=fastqFname, pathName=fastqPathName, md5sum=md5sum, fastqStats=newSeqStats, OneKP_Sequence=newDatasetSequence)
                elif dataOriginAcronym == 'AG':
                    # Also calculate number of genes and sum length of contigs in the raw fasta file and add to the newSeqsStats object:
                    numbrSequences, sumLengthOfContigs = rawFilenameStats(filename=fastqFname)
                    newSeqStats.numbrRecords = numbrSequences
                    newSeqStats.sumLengthOfSeqs = sumLengthOfContigs
                    newInputSequence = paftol.database.analysis.InputSequence(dataOrigin=dataOrigin, filename=fastqFname, pathName=fastqPathName, md5sum=md5sum, fastqStats=newSeqStats, annotatedGenome=newDatasetSequence)
                elif dataOriginAcronym == 'UG':
                    # Also calculate number of genes and sum length of contigs in the raw fasta file and add to the newSeqsStats object:
                    numbrSequences, sumLengthOfContigs = rawFilenameStats(filename=fastqFname)
                    newSeqStats.numbrRecords = numbrSequences
                    newSeqStats.sumLengthOfSeqs = sumLengthOfContigs
                    newInputSequence = paftol.database.analysis.InputSequence(dataOrigin=dataOrigin, filename=fastqFname, pathName=fastqPathName, md5sum=md5sum, fastqStats=newSeqStats, UnannotatedGenome=newDatasetSequence)
                #print dir(newInputSequence)
                #print "1.Looking at newInputSequence contents: ", newInputSequence.filename
                newInputSequenceList.append(newInputSequence)
            else:
                # Paul B. - if fastqFile.md5sum == md5sum:
                if inputSequence.md5sum == md5sum:
                    logger.info('fastq file %s already in database, verified md5sum', fastqFname)
                else:
                    # Paul B. - raise StandardError, 'fastq file %s in database with md5sum = %s, but found md5sum = %s' % (fastqFname, fastqFile.md5sum, md5sum)
                    raise StandardError, 'fastq file %s in database with md5sum = %s, but found md5sum = %s' % (fastqFname, inputSequence.md5sum, md5sum)
        else:
### Need to change message - for other datatypes
            logger.warning('No sample identifier obtainable for %s', fastqFname)
    # Paul B. - transactionSuccessful = insertPaftolFastqFileList(connection, newPaftolFastqFileList)
    # NB - as the code was, it seems that the insertInputSequenceList() method was entered even though the files were
    #      already in the db. If so, newInputSequenceList would be empty so I thought a conditional was required.
    #      However the method will just fall silent if list is empty.
    #      Can't really test the array in the method though, will get e.g. index out of range error.
    # However, newPaftolSequence (now generically newDatasetSequence) is now included in the method so it will not fall silent if list is empty, so need a conditional now so as not to enter insertInputSequenceList method if list is empty.
    # i.e. if there are no filenames to include, the insertInputSequenceList method will just add the sampleIds and nothign else - OK?
    if newInputSequenceList:
        ###print "2.Looking at array contents: ", newInputSequenceList[0].filename    # gives error unless occupied
        ###dir(newInputSequenceList[0])     # gives error unless occupied
        transactionSuccessful = insertInputSequenceList(connection, analysisDatabase, newInputSequenceList, newDatasetSequence, externalGenesFile)
                                                                                                                                ### check whether addExternalGenes can be empty
        return transactionSuccessful
    else:
        logger.warning('Not attempting to insert filename info into database')


def findGeneType(analysisDatabase, geneTypeName):
    for geneType in analysisDatabase.geneTypeDict.values():
        if geneType.geneTypeName == geneTypeName:
            return geneType
    return None


def findDataOrigin(analysisDatabase, dataOriginAcronym):
    for dataOrigin in analysisDatabase.dataOriginDict.values():
        if dataOrigin.acronym == dataOriginAcronym:
            return dataOrigin
    return None


def addTargetsFile(targetsFname, fastaPath=None, description=None, insertGenes=False, geneTypeName=None):
    if insertGenes and geneTypeName is None:
        raise StandardError, 'illegal state: insertion of new genes requested but no gene type name given'
     # Paul B. - recreated the full path to the fasta target file (assumed to be within the paftol dir):
    if fastaPath is not None:       
        fastaPath = fastaPath + '/' + targetsFname
        print 'fastqPath: ', fastaPath
    else:
        fastaPath = None
        print 'fastaPath: ', fastaPath
    md5sum = paftol.tools.md5HexdigestFromFile(targetsFname)        
    paftolTargetSet = paftol.PaftolTargetSet()         # A PaftolTargetSet object
    paftolTargetSet.readFasta(targetsFname)
    numSequences = len(paftolTargetSet.getSeqRecordList())
    analysisDatabaseDetails = getAnalysisDatabaseDetails()
    connection = analysisDatabaseDetails.makeConnection()
    analysisDatabase = paftol.database.analysis.AnalysisDatabase(connection)
    geneType = findGeneType(analysisDatabase, geneTypeName)     # Paul B. - a geneType object
    if geneType is None:
        connection.close()
        raise StandardError, 'unknown gene type: %s' % geneTypeName
    targetsFile = findFastaFile(analysisDatabase, targetsFname)
    if targetsFile is not None:
        connection.close()
        if targetsFile.md5sum == md5sum:
            logger.info('targets file %s already in database, verified md5sum', targetsFname)
            #print 'ok'
            return
        else:
            raise StandardError, 'targets file %s in database with md5sum = %s, but found md5sum = %s' % (targetsFname, targetsFile.md5sum, md5sum)
    # connection.start_transaction(isolation_level='REPEATABLE READ', readonly=False)
    # Paul B. - changed to fit with the auto_increment
    #targetsFastaFile = paftol.database.analysis.FastaFile(None, targetsFname, md5sum, None, numSequences)
    # Paul B. - now sending to the ReferenceTarget table:
    #targetsFastaFile = paftol.database.analysis.FastaFile(targetsFname, None, md5sum, None, numSequences)
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
        # Paul B. - changed to fit with the auto_increment - 25.5.2020 - also added the exemplarGeneId:
        #newPaftolGene = paftol.database.analysis.PaftolGene(None, missingGeneName, geneType)
        newPaftolGene = paftol.database.analysis.PaftolGene(missingGeneName, geneType, None)
        newPaftolGeneList.append(newPaftolGene)
    paftolGeneDict = {}
    for newPaftolGene in newPaftolGeneList:
        paftolGeneDict[newPaftolGene.geneName] = newPaftolGene
    for paftolGene in analysisDatabase.paftolGeneDict.values():
        paftolGeneDict[paftolGene.geneName] = paftolGene
    referenceTargetList = []
    for paftolGene in paftolTargetSet.paftolGeneDict.values():
        for paftolTarget in paftolGene.paftolTargetDict.values():
            # Paul B. - changed to fit with the auto_increment:
            #referenceTargetList.append(paftol.database.analysis.ReferenceTarget(None, paftolGeneDict[paftolGene.name], paftolTarget.organism.name, len(paftolTarget.seqRecord), targetsFastaFile))
            # Paul B. - now changing to add the targets file info:
            #referenceTargetList.append(paftol.database.analysis.ReferenceTarget(paftolGeneDict[paftolGene.name], paftolTarget.organism.name, len(paftolTarget.seqRecord), targetsFastaFile))
#### 27.2.2020 - ACTUALLY NOW planning TO CHANGE AGAIN AND CREATE A SEPARATE TABLE TO HOUSE THE FASTA FILE
            referenceTargetList.append(paftol.database.analysis.ReferenceTarget(paftolGene=paftolGeneDict[paftolGene.name], paftolOrganism=paftolTarget.organism.name, paftolTargetLength=len(paftolTarget.seqRecord), \
            targetsFastaFile=targetsFname, targetsFastaFilePathName=fastaPath, md5sum=md5sum, numTargetSequences=numSequences))
    transactionSuccessful = False
    # Paul B. - removed table locking and introduced auto_increment for each primary key:
    #lockCursor = connection.cursor()
    #lockCursor.execute('LOCK TABLE FastaFile WRITE, FastqFile WRITE, FastqStats WRITE, GeneType WRITE, PaftolGene WRITE, ReferenceTarget WRITE')
    try:
        logger.info('adding new targets file %s', targetsFname)
        cursor = connection.cursor(prepared=True)
        try:
            for newPaftolGene in newPaftolGeneList:
                #newPaftolGene.id = generateUnusedPrimaryKey(cursor, 'PaftolGene')
                newPaftolGene.insertIntoDatabase(cursor)
                newPaftolGene.id = cursor.lastrowid
                #print "newPaftolGene.id: ", newPaftolGene.id
            #targetsFastaFile.id = generateUnusedPrimaryKey(cursor, 'FastaFile')
            # Paul B. - changed to send to ReferenceTarget table instead:
            #targetsFastaFile.insertIntoDatabase(cursor)
            #targetsFastaFile.id = cursor.lastrowid
            #print "targetsFastaFile.id: ", targetsFastaFile.id
            for referenceTarget in referenceTargetList:
                #referenceTarget.id = generateUnusedPrimaryKey(cursor, 'ReferenceTarget')
                referenceTarget.insertIntoDatabase(cursor)
                referenceTarget.id = cursor.lastrowid
                #print "referenceTarget.id: ", referenceTarget.id
            connection.commit()
            transactionSuccessful = True
        finally:
            if not transactionSuccessful:
                connection.rollback()
                 # Paul B added:
                print "ERROR: commit unsucessful for newPaftolGene.id: ", newPaftolGene.id
                #print "ERROR: commit unsucessful for targetsFastaFile.id: ", targetsFastaFile.id
                print "ERROR: commit unsucessful for referenceTarget.id: ", referenceTarget.id
            cursor.close()
    finally:
        if not transactionSuccessful:
            connection.rollback()
        #lockCursor.execute('UNLOCK TABLES')
        #lockCursor.close()
    connection.close()
    return transactionSuccessful


def findFastqFiles(analysisDatabase, result):
    fwdFastqFname = os.path.basename(result.forwardFastq)
    if result.reverseFastq is not None:
        revFastqFname = os.path.basename(result.reverseFastq)
    fwdFastqFile = None
    revFastqFile = None
    # Paul B. - changed to use inputSequenceDict - each element contains a row object:
    #for fastqFile in analysisDatabase.fastqFileDict.values():
    for fastqFile in analysisDatabase.inputSequenceDict.values():
        if fastqFile.filename == fwdFastqFname:
            fwdFastqFile = fastqFile
        if result.reverseFastq is not None:
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
    analysisDatabaseDetails = getAnalysisDatabaseDetails()      ### PaulB - returns a mysql.connector connection object
    
    connectionSuccessful = False                              # Paul B added
    connectionCountr = 1                                      # Paul B added
    while connectionSuccessful == False:                      # Paul B added
        try:                                                  # Paul B added
            logger.warning('Connection attempt %s ...', connectionCountr) # Paul B added
            connection = analysisDatabaseDetails.makeConnection()
            #connection.autocommit = True                               ### Paul B. - tried autocommit
            analysisDatabase = paftol.database.analysis.AnalysisDatabase(connection)
            # Paul B. - 25.2.2020 - now acesses the ReferenceTarget table instead:
            ### NBNB - this is not good I think but it works - targetsFastaFile needs to have its own table
            targetsFastaFile = findFastaFile(analysisDatabase, result.paftolTargetSet.fastaHandleStr)
            numMappedReads = len(result.paftolTargetSet.getMappedReadNameSet())
            numUnmappedReads = result.paftolTargetSet.numOfftargetReads
            if targetsFastaFile is None:
                # raise StandardError, 'targets file "%s" not in database' % result.paftolTargetSet.fastaHandleStr
                logger.info('unknown targets file "%s" -- continuing')
            # Paul B. - alter method find fastq files in the InputSequence table:
            fwdFastqFile, revFastqFile = findFastqFiles(analysisDatabase, result)
            if fwdFastqFile is None:
                raise StandardError, 'forward fastq file "%s" not in database' % result.forwardFastq
            if revFastqFile is None:
                logger.info('reverse fastq file not in database')
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
            ### Paul B - changed use auto_increment and to add in the CDS fasta filename rather than all the contigs (required a new Paftol.HybpiperResult object variable)
            #contigFastaFile = None
            reconstructedCdsFastaFname = None
            #if result.contigFastaFname is not None:
            #if result.reconstructedCdsFastaFname is not None:
                #contigFastaFile = paftol.database.analysis.FastaFile(None, result.contigFastaFname, paftol.tools.md5HexdigestFromFile(result.contigFastaFname), None, len(paftol.tools.fastaSeqRecordList(result.contigFastaFname)))
                # Paul B. - removed FastaFile, fasta contig file now goes into ContigRecovery table:
                #contigFastaFile = paftol.database.analysis.FastaFile(result.reconstructedCdsFastaFname, result.reconstructedCdsFastaFnamePath, paftol.tools.md5HexdigestFromFile(result.reconstructedCdsFastaFname), None, len(paftol.tools.fastaSeqRecordList(result.reconstructedCdsFastaFname)))
            ### Paul B - removed id=None first parameter to fit with the auto_increment change:
            #print "testing targetsFastaFile: ", result.paftolTargetSet.fastaHandleStr
            #contigRecovery = paftol.database.analysis.ContigRecovery(fwdFastq=fwdFastqFile, revFastq=revFastqFile, fwdTrimmedFastqStats=trimmedForwardFastqStats, revTrimmedFastqStats=trimmedReverseFastqStats, contigFastaFile=contigFastaFile, targetsFastaFile=targetsFastaFile, numMappedReads=numMappedReads, numUnmappedReads=numUnmappedReads, softwareVersion=paftol.__version__, cmdLine=result.cmdLine)
            # Paul B. - contig file info now goes into ContigRecovery table:
            contigRecovery = paftol.database.analysis.ContigRecovery(fwdFastq=fwdFastqFile, revFastq=revFastqFile, \
            fwdTrimmedFastqStats=trimmedForwardFastqStats, revTrimmedFastqStats=trimmedReverseFastqStats, \
            contigFastaFileName=result.reconstructedCdsFastaFname, contigFastaFilePathName=result.reconstructedCdsFastaFnamePath, contigFastaFileMd5sum=paftol.tools.md5HexdigestFromFile(result.reconstructedCdsFastaFname), \
            referenceTarget=targetsFastaFile, \
            numMappedReads=numMappedReads, numUnmappedReads=numUnmappedReads, softwareVersion=paftol.__version__, cmdLine=result.cmdLine)
            recoveredContigList = []
            ### Paul B. - Added info from result.reconstructedCdsDict to RecoveredContig table instead.
            ### NB - result.contigDict[geneName] is a list of contig BioSeqRecord objects but
            ### but result.reconstructedCdsDict[geneName] just contains a single supercontig BioSeqRecord object NOT in a list
            ### so don't need the contig for loop.
            #for geneName in result.contigDict:
            #   if result.contigDict[geneName] is not None and len(result.contigDict[geneName]) > 0:  # I think this is the numbr  of Seq records for this gene!!!
            #        for contig in result.contigDict[geneName]:
                        # I think this accesses the id and seq values (not tested)
                        #print "Contig.id ", contig.id   
                        #print "contig.seq", contig.seq

            RC_Countr = 0   # Counting the number of recovered contigs so I can compare and check it with the value calculated by the db (to check that auto_increment is working)
            for geneName in result.reconstructedCdsDict:
                if result.reconstructedCdsDict[geneName] is not None and len(result.reconstructedCdsDict[geneName]) > 0:    # 8.9.2020 - This is the length of the seq, not the eqivalent of len(result.contigDict[geneName]) > 0!!!! Still OK I think?
                    #print "result.reconstructedCdsDict[geneName].id: ", result.reconstructedCdsDict[geneName].id
                    #print "Length of seq:", len(result.reconstructedCdsDict[geneName])
                    #print "Seq:", result.reconstructedCdsDict[geneName].seq
                    representativeReferenceTarget = findReferenceTarget(analysisDatabase, geneName, result.representativePaftolTargetDict[geneName].organism.name)
                    if representativeReferenceTarget is None:
                        raise StandardError, 'unknown reference target for geneName = %s, organismName = %s' % (geneName, result.representativePaftolTargetDict[geneName].organism.name)

                    # Store the md5sum of each sequence for sample-gene version control:
                    md5Seq = hashlib.md5()
                    md5Seq.update(str(result.reconstructedCdsDict[geneName].seq))
                    md5sum = md5Seq.hexdigest()
                    #print 'seq record md5sum:', md5sum

                    ### Paul B - removed 'None' first parameter to fit with the auto_increment change; changed from len(contig) to len(result.reconstructedCdsDict[geneName])
                    ###          Also started to use argument=value format
                    #recoveredContig = paftol.database.analysis.RecoveredContig(contigRecovery, paftolGeneDict[geneName], len(result.reconstructedCdsDict[geneName]), representativeReferenceTarget)
                    recoveredContig = paftol.database.analysis.RecoveredContig(contigRecovery=contigRecovery, paftolGene=paftolGeneDict[geneName], seqLength=len(result.reconstructedCdsDict[geneName]), md5sum=md5sum, representativeReferenceTarget=representativeReferenceTarget)
                    recoveredContigList.append(recoveredContig)
                    RC_Countr += 1
            contigRecovery.numRecoveredContigsCheck = RC_Countr
            transactionSuccessful = False
            ### Paul B. - can now remove table locking because now using auto_increment
            #lockCursor = connection.cursor(prepared=False)
            #lockCursor.execute('LOCK TABLE FastaFile WRITE, FastqFile WRITE, FastqStats WRITE, ContigRecovery WRITE, RecoveredContig WRITE')
            try:
                cursor = connection.cursor(prepared=True)
                try:
                    ###  Paul B. - making changes to use auto_increment:
                    if trimmedForwardFastqStats is not None:
                        #trimmedForwardFastqStats.id = generateUnusedPrimaryKey(cursor, 'FastqStats')
                        trimmedForwardFastqStats.insertIntoDatabase(cursor)
                        trimmedForwardFastqStats.id = cursor.lastrowid
                        if trimmedForwardFastqStats.id is not None:
                            print "trimmedForwardFastqStats.id: ", trimmedForwardFastqStats.id
                    if trimmedReverseFastqStats is not None:
                        #trimmedReverseFastqStats.id = generateUnusedPrimaryKey(cursor, 'FastqStats')
                        trimmedReverseFastqStats.insertIntoDatabase(cursor)
                        trimmedReverseFastqStats.id = cursor.lastrowid
                        if trimmedReverseFastqStats.id is not None:
                            print "trimmedReverseFastqStats.id: ", trimmedReverseFastqStats.id
                    # Paul B. - contigFastaFile now goes into ContigRecovery table (no need for conditional either? contigFastaFile value should just remain NULL
                    #if contigFastaFile is not None:
                        #contigFastaFile.id = generateUnusedPrimaryKey(cursor, 'FastaFile')
                        #contigFastaFile.insertIntoDatabase(cursor)
                        #contigFastaFile.id = cursor.lastrowid
                        #print "contigFastaFile.id: ", contigFastaFile.id
                    #contigRecovery.id = generateUnusedPrimaryKey(cursor, 'ContigRecovery')
                    contigRecovery.insertIntoDatabase(cursor)
                    contigRecovery.id = cursor.lastrowid
                    if contigRecovery.id is not None:
                        print "contigRecovery.id: ", contigRecovery.id
                    for recoveredContig in recoveredContigList:
                        #recoveredContig.id = generateUnusedPrimaryKey(cursor, 'RecoveredContig')
                        recoveredContig.insertIntoDatabase(cursor)
                        recoveredContig.id = cursor.lastrowid
                        if recoveredContig.id is not None:
                            print "recoveredContig.id: ", recoveredContig.id
                        #time.sleep(0.06)
                    #### time delay - 1 second doen to 60milsec
                    connection.commit()
                    transactionSuccessful = True
                finally:
                    if not transactionSuccessful:
                        connection.rollback()
                        # Paul B added:                                                                         # NB - these variables may not exist if commit fails
                        print "ERROR: commit unsucessful for contigRecovery.id: "                               #, contigRecovery.id
                        print "ERROR: commit unsucessful for trimmedForwardFastqStats.id: "                     #, trimmedForwardFastqStats.id
                        print "ERROR: commit unsucessful for trimmedReverseFastqStats.id: "                     #, trimmedReverseFastqStats.id
                        #print "ERROR: commit unsucessful for contigFastaFile.id: ", contigFastaFile.id
                        print "ERROR: commit unsucessful for recoveredContig.id (for last row created): "       #, recoveredContig.id
                    cursor.close()
            finally:
                if not transactionSuccessful:
                    connection.rollback()
                    # Paul B. removed this again- it should appear above:
                    #print "ERROR: commit unsucessful for contigRecovery.id: ", contigRecovery.id
                #lockCursor.execute('UNLOCK TABLES')
                #lockCursor.close()
            connection.close()
            connectionSuccessful = True    # Paul B added
            logger.warning('connectionSuccessful == %s', connectionSuccessful) # Paul B added
        except (mysql.connector.errors.InterfaceError, mysql.connector.errors.InternalError, mysql.connector.errors.DatabaseError): # Paul B added
            timeToSleep = random.randint(1,25)                                                                                      # Paul B added
            print "Sample unable to connect to the database - retrying in ", timeToSleep, "seconds..."                              # Paul B added
            time.sleep(timeToSleep)                                                 # Paul B added
            connectionCountr += 1                                                   # Paul B added
            if connectionCountr > 10:                                               # Paul B added 
                logger.warning('Tried to connect to database 10 times, giving up!') # Paul B added 
                return  transactionSuccessful                                       # Paul B added
    return transactionSuccessful


def addExternalGenes(cursor=None, analysisDatabase=None, inputSequence=None, newDatasetSequence=None, externalGenesFile=None):

    ''' Adds required info for a recovered genes fasta file into the paftol_da.ContigRecovery and paftol_da.RecoveredContig tables.
        Relevant to OneKP_Transcript and AnnotatedGenome samples ONLY

        Input parameters: cursor from the database connection object, analysisDatabase object, input sequence row object and the name of an external recovered genes file

        Note: There is no check to see whether a recovered genes fasta file present (c.f. gene recovery from fastq files - preRecoveryCheck method).
        However in this case, the recovered genes fasta file is being loaded at the same time as the sample is being added for the first time (which is checked),
        so the recovered genes file can't already be there. 
    '''
    #print paftol.tools.md5HexdigestFromFile(externalGenesFile)


    # Get seq records from externalGenesFile into a dict.
    ### NB - numbrSequences and sumLengthOfContigs not required here - can remove
    numbrSequences = 0
    sumLengthOfContigs = 0
    seqRecords = {}     #

    # The fasta file has to be unzipped for the Bio.SeqIO to work - it does not throw an error
    # if file is zipped and produces gobbledeegook values for numbrSequences, sumLengthOfContigs !
    # Attempting to check if file is zipped here:
    if re.search('.gz$|.bz2$', externalGenesFile) is not None:
        raise StandardError, 'Raw data fasta file is zipped, needs to be unzipped: %s' % externalGenesFile
    for record in Bio.SeqIO.parse(externalGenesFile, "fasta"):     # returns a SeqRecord object, includes a Seq object called seq
        #print record.id, "\n", record.seq, len(record)
        sumLengthOfContigs = sumLengthOfContigs + len(record)
        numbrSequences += 1
        seqRecords[record.id] = record
    logger.info('Recovered genes fasta file: numbrSequences=%s; sumLengthOfContigs=%s', numbrSequences, sumLengthOfContigs)


    # Get all gene names present in paftol_da database:
    paftolGeneDict = {}
    for paftolGene in analysisDatabase.paftolGeneDict.values():     # returns a PaftolGene row object
        paftolGeneDict[paftolGene.geneName] = paftolGene
    # Cross-check gene names input gene contigs fasta file:
    for geneName in seqRecords:
        if geneName not in paftolGeneDict:
            raise StandardError, 'found gene %s in result but it is not in the analysis database' % geneName


    # Populate the ContigRecovery object.
    contigFastaFileName = os.path.basename(externalGenesFile)
    #contigFastaFilePathName = os.path.dirname(externalGenesFile)   # Using the full path for the moment
    contigRecovery = paftol.database.analysis.ContigRecovery( \
    fwdFastq=inputSequence, \
    revFastq=None, \
    fwdTrimmedFastqStats=None, revTrimmedFastqStats=None, \
    contigFastaFileName=contigFastaFileName, \
    contigFastaFilePathName=externalGenesFile, \
    contigFastaFileMd5sum=paftol.tools.md5HexdigestFromFile(externalGenesFile), \
    referenceTarget=None, \
    numMappedReads=None, numUnmappedReads=None, softwareVersion=None, cmdLine=None)


    # Populate the RecoveredContig object
    recoveredContigList = []
    RC_Countr = 0   # Counting the number of recovered contigs so I can compare and check it with the value calculated by the db (to check that auto_increment is working)
    for geneName in seqRecords:
        if seqRecords[geneName] is not None and len(seqRecords[geneName]) > 0:
            #print "geneName: ", geneName
            #print 'Sequence: ', seqRecords[geneName].seq
            #print 'Sequence: ', seqRecords[geneName].description

            # Parse the organism name from the fasta header (2nd field in record.description in the foramat returned by retrieveTargets):
            headrFields = seqRecords[geneName].description.split()
            #print 'headrFields:', headrFields[1]
            # Get the representative reference target:
            representativeReferenceTarget = findReferenceTarget(analysisDatabase, geneName, headrFields[1])
            if representativeReferenceTarget is None:
                raise StandardError, 'unknown reference target for geneName = %s, organismName = %s' % (geneName, headrFields[1])

            # Store the md5sum of each sequence for sample-gene version control:
            md5Seq = hashlib.md5()
            md5Seq.update(str(seqRecords[geneName].seq))
            md5sum = md5Seq.hexdigest()
            #print 'seq record md5sum:', md5sum

            recoveredContig = paftol.database.analysis.RecoveredContig(contigRecovery=contigRecovery, paftolGene=paftolGeneDict[geneName], seqLength=len(seqRecords[geneName]), md5sum=md5sum, representativeReferenceTarget=representativeReferenceTarget)
            recoveredContigList.append(recoveredContig)
            RC_Countr += 1
    contigRecovery.numRecoveredContigsCheck = RC_Countr
    #print 'contigRecovery.numRecoveredContigsCheck: ', contigRecovery.numRecoveredContigsCheck


    # Start upload to paftol_da db:
    contigRecovery.insertIntoDatabase(cursor)
    contigRecovery.id = cursor.lastrowid
    if contigRecovery.id is not None:
        print "contigRecovery.id: ", contigRecovery.id
        for recoveredContig in recoveredContigList:
            recoveredContig.insertIntoDatabase(cursor)
            recoveredContig.id = cursor.lastrowid
            if recoveredContig.id is not None:
                print "recoveredContig.id: ", recoveredContig.id

