import unittest

from experimentsLogReader import ExperimentLogReader


class logreaderTest(unittest.TestCase):
    
    def test_getParametrs(self):
        testLogReader = ExperimentLogReader("logs/" + "test.log", "prettyLogs/")
        testLogsDict = testLogReader.getLogs()
        self.assertEqual(testLogsDict["location"], 'IRBENE')
        self.assertEqual(testLogReader.dates, '05 Nov 2017')
        self.assertEqual(testLogReader.scan_names, ['1', '2'])
        self.assertEqual(testLogReader.scanLines[1], ['2017.309.08:17:01.04:scan_name=no0001,m86,ir,420,420\n', '2017.309.08:17:01.04:source=abc,061437.05,133936.2,2000.0,\n'])
        self.assertEqual(testLogReader.scanLines[2], ['2017.309.08:17:01.04:scan_name=no0002,m86,ir,420,420\n', '2017.309.08:17:01.04:source=qwe,061437.05,133936.2,2000.0,\n'])

if __name__ == '__main__':
    unittest.main()