import unittest

import paftol
import paftol.msarunner


class MultipleSequenceAlignmentRunnerTestCase(unittest.TestCase):
    
    def test_instantiation(self):
        multipleSequenceAlignmentRunner = paftol.msarunner.MultipleSequenceAlignmentRunner()
        self.assertIsNotNone(multipleSequenceAlignmentRunner)

        
class MafftRunnerTestCase(unittest.TestCase):
    
    def test_subprocess(self):
        mafftRunner = paftol.msarunner.MafftRunner()
        subprocess = mafftRunner.makeSubprocess()


