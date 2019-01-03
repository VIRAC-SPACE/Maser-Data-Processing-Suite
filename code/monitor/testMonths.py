import unittest

from months import Months

class TestConfigParser(unittest.TestCase):
    
    def test_getMonthsNumber(self):
        month = Months()
        self.assertEqual(month.getMonthNumber("Jan"), '1')
        self.assertEqual(month.getMonthNumber("Jun"), '6')
        self.assertEqual(month.getMonthNumber("Dec"), '12')
        
if __name__ == '__main__':
    unittest.main()