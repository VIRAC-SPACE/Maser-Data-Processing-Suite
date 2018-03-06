import unittest

from experimentsLogReader import ExperimentLogReader


class logreaderTest(unittest.TestCase):
    
    def test_getParametrs(self):
        testLogReader = ExperimentLogReader("logs/" + "test.log", "prettyLogs/")
        testLogsDict = testLogReader.getLogs()
        print testLogsDict
        self.assertEqual(testLogsDict[1]["source"], 'abc')

if __name__ == '__main__':
    unittest.main()