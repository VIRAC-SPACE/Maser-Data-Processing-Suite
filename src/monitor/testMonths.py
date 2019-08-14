import unittest

from months import Months

class TestConfigParser(unittest.TestCase):
    
    def test_getMonthsNumber(self):
        month = Months()
        self.assertEqual(month.getMonthNumber("Jan"), '1')
        self.assertEqual(month.getMonthNumber("Feb"), '2')
        self.assertEqual(month.getMonthNumber("Mar"), '3')
        self.assertEqual(month.getMonthNumber("Apr"), '4')
        self.assertEqual(month.getMonthNumber("May"), '5')
        self.assertEqual(month.getMonthNumber("Jun"), '6')
        self.assertEqual(month.getMonthNumber("Jul"), '7')
        self.assertEqual(month.getMonthNumber("Aug"), '8')
        self.assertEqual(month.getMonthNumber("Sep"), '9')
        self.assertEqual(month.getMonthNumber("Oct"), '10')
        self.assertEqual(month.getMonthNumber("Nov"), '11')
        self.assertEqual(month.getMonthNumber("Dec"), '12')
       
    def test_MonoState(self):
        m1 = Months()
        m2 = Months()
        
        self.assertEqual(m1.getMonthNumber("Dec"), m2.getMonthNumber("Dec"))
                
if __name__ == '__main__':
    unittest.main()