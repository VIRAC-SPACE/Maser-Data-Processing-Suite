import datetime

class MonitoringViewHelper():
        
        @staticmethod 
        def formatDate(xdata, index):
            date = xdata[index][0].strftime("%H %M %S %d %m %Y").split()
            month = datetime.date(1900, int(date[-2]), 1).strftime('%B')[0:3].title().replace("ū", "u").replace("i", "y").replace("k", "c")
            date[-2] = month
            date = "_".join(date)
            return date
        
        @staticmethod 
        def getIteration(iteration_list, index):
            return str(iteration_list[index])
        
        @staticmethod 
        def getLocation(location_list, index):
            return location_list[index]