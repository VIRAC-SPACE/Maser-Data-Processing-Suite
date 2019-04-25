class Result(object):
    __slots__ = ('location', 'date', 'polarizationU1', 'polarizationU9', 'polarizationAVG', 'iteration_number', 'specie', 'type', 'modifiedJulianDays', 'areas')
    def __init__(self, location, date, polarizationU1, polarizationU9, polarizationAVG, iteration_number, specie, type, modifiedJulianDays, areas):
        self.location = location
        self.date = date
        self.polarizationU1 = polarizationU1
        self.polarizationU9 = polarizationU9
        self.polarizationAVG = polarizationAVG
        self.iteration_number = iteration_number
        self.specie = specie
        self.type = type
        self.modifiedJulianDays = modifiedJulianDays
        self.areas = areas
        
    def __iter__(self):
        yield 'location', self.location
        yield 'date', self.date
        yield 'polarizationU1', self.polarizationU1
        yield 'polarizationU9', self.polarizationU9
        yield 'polarizationUAVG', self.polarizationAVG
        yield 'iteration_number', self.iteration_number
        yield 'specie', self.specie
        yield 'type', self.type
        yield 'modifiedJulianDays', self.modifiedJulianDays
        yield 'areas', self.areas
