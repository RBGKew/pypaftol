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
        self.assertEqual(None, paftolDpDetails.dbusername)
        self.assertEqual(None, paftolDpDetails.dbpassword)
        self.assertEqual(None, paftolDpDetails.dbhost)
        self.assertEqual(None, paftolDpDetails.dbname)
        # f = open(...)
        # paftolDbDetails.readFile(f)
