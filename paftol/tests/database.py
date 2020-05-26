import sys
import unittest
import copy
import cStringIO

import Bio
import Bio.Seq
import Bio.SeqRecord

import paftol
import paftol.tools
import paftol.database


class DatabaseTestCase(unittest.TestCase):
    
        
    def test_database(self):
        paftolDatabaseDetails = paftol.database.PaftolDatabaseDetails()
        self.assertIsNone(None, paftolDatabaseDetails.dbusername)
        self.assertIsNone(None, paftolDatabaseDetails.dbpassword)
        self.assertIsNone(None, paftolDatabaseDetails.dbhost)
        self.assertIsNone(None, paftolDatabaseDetails.dbname)
        f = cStringIO.StringIO("""user: paftol
password: topsecret
host: localhost
database: paftol
""")
        paftolDatabaseDetailsFromFile = paftol.database.PaftolDatabaseDetails(detailsFile=f)
        self.assertEqual('paftol', paftolDatabaseDetailsFromFile.dbusername)
        self.assertEqual('topsecret', paftolDatabaseDetailsFromFile.dbpassword)
        self.assertEqual('localhost', paftolDatabaseDetailsFromFile.dbhost)
        self.assertEqual('paftol', paftolDatabaseDetailsFromFile.dbname)
        # f = open(...)

    def test_ExistingFastqFile(self):
        paftolPrefixedFname = './miseq/post_pilot_run/SP0052/PAFTOL-006707_R1_001.fastq.gz'
        invalidPaftolPrefixedFname = './miseq/post_pilot_run/SP00521/PAFTOL-006707_R1_001.fastq.gz'
        orchidaceaeFname = '../miseq/post_pilot_run/SP0018-Orchidaceae/Myoxanthus-ceratothallis_S59_L001_R1_001.fastq.gz'
        e = paftol.database.ExistingFastqFile(paftolPrefixedFname)
        self.assertEqual(52, e.findSequencingPoolNumber())
        e = paftol.database.ExistingFastqFile(invalidPaftolPrefixedFname)
        self.assertIsNone(e.findSequencingPoolNumber())
        e = paftol.database.ExistingFastqFile(orchidaceaeFname)
        self.assertEqual(18, e.findSequencingPoolNumber())
        

# paftolDatabaseDetails = paftol.database.PaftolDatabaseDetails(detailsFile=open('~/.paftol/dpdatabase.cfg'))
# connection = paftolDatabaseDetails.makeConnection()
