class Months():
    __months = dict()
    __initialized = False
    
    def __init__(self):
        if Months.__initialized == False:
            Months.__initialized = True
            Months.__months = {"Jan":"1", "Feb":"2", "Mar":"3", "Apr":"4", "May":"5", "Jun":"6", "Jul":"7", "Aug":"8", "Sep":"9", "Oct":"10", "Nov":"11", "Dec":"12"}
            
    def getMonthNumber(self, month):
        return  Months.__months[month]
