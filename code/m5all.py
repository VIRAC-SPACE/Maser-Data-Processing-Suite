#!/usr/bin/python

import os
import subprocess
import os.path


M5A_PATH = "/datadisk/new_procid/14mar"
PATH = "/datadisk/new_procid/14mar"
SCH_PATH = "/datadisk/new_procid/14mar"

for file in os.listdir(M5A_PATH):
	if file.endswith(".m5a"):
		m5aFile = os.path.join(M5A_PATH, file)
		fileSplit = file.split("_")
		expName = fileSplit[0]
		scanNo = fileSplit[2].split(".")[0][-2:]
		datFile = SCH_PATH+"%s_n%s.dat" % (expName, scanNo)
		if os.path.isfile(datFile):
			print "File %s exists, m5spec will not be called!" % (datFile)
			continue
		else:
			print "m5spec will be called for file %s!" % (datFile)
		p = subprocess.Popen(["m5spec", m5aFile, "VDIF_8000-32-4-2", "4096", "100000000000000", datFile, "0", "-nopol"])
		p.wait()
		
