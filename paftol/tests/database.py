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
        f = cStringIO.StringIO("""username: paftol
password: topsecret
host: localhost
dbname: paftol
""")
        paftolDatabaseDetailsFromFile = paftol.database.PaftolDatabaseDetails(detailsFile=f)
        self.assertEqual('paftol', paftolDatabaseDetailsFromFile.dbusername)
        self.assertEqual('topsecret', paftolDatabaseDetailsFromFile.dbpassword)
        self.assertEqual('localhost', paftolDatabaseDetailsFromFile.dbhost)
        self.assertEqual('paftol', paftolDatabaseDetailsFromFile.dbname)
        # f = open(...)

# paftolDatabaseDetails = paftol.database.PaftolDatabaseDetails(detailsFile=open('~/.paftol/dpdatabase.cfg'))
# connection = paftolDatabaseDetails.makeConnection()
