import sys
import unittest

import paftol
import paftol.tests


testLoader = unittest.TestLoader()
suite = testLoader.loadTestsFromModule(paftol.tests)
testResult = unittest.TestResult()
suite.run(testResult)
sys.stdout.write('%s\n' % str(testResult))
