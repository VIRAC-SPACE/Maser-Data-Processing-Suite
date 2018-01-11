To RUN SPELI

Observation is prepared using sched.
For methanol observation the following code name is used 
m<NE> - where 'm' stands for methanol, and <NE> is the experiment number, e.g. m74 
The scan number is designated by
m<NS> - where 'n' stands for word number, and <NS> is the corresponding scan number, e.g. n37

The data file name is required to be in form (e.g. m74_n37.dat)
m<NE>_n<NS>.dat

1) Before running the speli collect information from the corresponding experiment station (Field System) log file *ir.log file  (e.g. m74ir.log)
by using experimentsLogReader.py, run: python experimentsLogReader.py m<NS>ir.log
Output will be m<NS>irsch.dat that will be read by speli automatically.
Note: This requires dopsetpy_v1.5 (compiled from f77 dopsetpy_v1.5.f)

EXAMPLE
$ python experimentsLogReader.py m74ir.log
OUTPUT m74irsch.dat 

2) Run python speli.py m74_n37.dat and do analyses using GUI.
Note: Using GUI with mouse and keybord - select bad points manually, cut out sygnal when modeling background etc

EXAMPLE
$ python speli.py m74_n37.dat
Output: m74_n37.dat.out, graphs, etc
