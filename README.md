# MDPS -  Maser Data Processing Suite

# Table of Contents

## Dependencies
# Table of Contents
* [Capabilities of MDPS](#capabilities-of-MDPS)
* [Dependencies](#dependencies)
* [Processing SDR output](#Processing-SDR-output)
* [Changelog](#changelog)
* [Getting Help](#getting-help)
* [Acknowledgements](#acknowledgements)

## Capabilities of MDPS
MDPS is maser data processing suite, with user friendly GUI. Currently it allows users to perform:

- Process SDR output
- Display monitoring  
- Visualize data for publications and presentation 

## Dependencies
- python 3
  - PyQt5 5.14.0 
  - Astropy 4.0.1.
  - coloredlogs 14.0
  - Matplotlib 3.2.1 
  - Numpy 1.18.2
  - H5py 2.10.0 
  - ExperimentsLogReader 2.2.2
  - SciPy 1.4.1
  - PeakUtils 1.3.3
  - pandas 1.0.3 
  - mplcursors 0.3
  - jplephem 2.14
  - configparser 5.0.0
  
## Processing SDR output
SDR for each scan creates four files **r0** **r1** **s0** **s1**. File name is calibrator_<SASID>_SURIs.txt

| **Scripts** | **Description** |
| --- | --- |
| main.py | Automatically call sdr_fs.py and total_spectrum_analyzer_qt5.py |
| sdr_fs.py | Process four output files from SDR |
| total_spectrum_analyzer_qt5.py | Process sdr_fs.py output|

## Changelog

Version 1.0

First release. Changes since version pre1.0 are:

## Getting Help

Bug reports, feature requests and make contributions (e.g. code patches) can be reported by opening a &quot;new issue&quot; ticket on GitHub. Please give as much information (e.g. the software component, version) as you can in the ticket. For bugs, it is extremely useful if a small self-contained code snippet that reproduces the problem is provided.

## Acknowledgements
This software was written by Janis Steinbergs. If you make use of this software to get results that appear in a publication or presentation please include this acknowledgement: &quot;We have made use of MDPS, a tool developed by Janis Steinbergs.&quot;
