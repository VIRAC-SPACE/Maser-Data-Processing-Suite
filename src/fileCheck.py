from  inotify import adapters
import os

notifier = adapters.Inotify()
notifier.add_watch("/home/janis/Documents/workspace-sts/DataProcessingForMaserObservation/dataFiles")

for event in notifier.event_gen():
    if event is not None:
        # print event      # uncomment to see all events generated
        if 'IN_CREATE' in event[1]:
            dataFile = event[3]
            logFile = dataFile.split("_")[0] + ".log"
            print "file '{0}' created in '{1}'".format(event[3], event[2])
            print  dataFile, logFile
            os.system("touch " + " /home/janis/Documents/workspace-sts/DataProcessingForMaserObservation/" + dataFile + logFile)