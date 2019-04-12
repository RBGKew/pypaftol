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
        paftolDpDetails = paftol.database.PaftolDpDetails()
        self.assertIsNone(None, paftolDpDetails.dbusername)
        self.assertIsNone(None, paftolDpDetails.dbpassword)
        self.assertIsNone(None, paftolDpDetails.dbhost)
        self.assertIsNone(None, paftolDpDetails.dbname)
        f = cStringIO.StringIO("""username: paftol
password: topsecret
host: localhost
dbname: paftol
""")
        paftolDpDetailsFromFile = paftol.database.PaftolDpDetails(detailsFile=f)
        self.assertEqual('paftol', paftolDpDetailsFromFile.dbusername)
        self.assertEqual('topsecret', paftolDpDetailsFromFile.dbpassword)
        self.assertEqual('localhost', paftolDpDetailsFromFile.dbhost)
        self.assertEqual('paftol', paftolDpDetailsFromFile.dbname)
        # f = open(...)
        # paftolDbDetails.readFile(f)
