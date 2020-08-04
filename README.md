# MDPS -  Maser Data Processing Suite

# Table of Contents

- [Capabilities of MDPS](#capabilities-of-mdps)
- [Dependencies](#dependencies)
- [Configuration of MDPS](#configuration-of-mdps)
- [Directory structure](#directory-structure)
- [Processing SDR output](#processing-sdr-output)
- [Monitoring](#monitoring)
- [Publication and presentation of observations](#publication-and-presentation-of-observations)
- [Changelog](#changelog)
- [Getting Help](#getting-help)
- [Acknowledgements](#acknowledgements)

## Capabilities of MDPS
MDPS is maser data processing suite, with user-friendly GUI. Currently, it allows users to perform:

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
  - tabulate 0.8.7
- de435.bsp(https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/de435.bsp)

## Configuration of MDPS
MDPS consist of two configuration files: 1) config.cfg, 2) plot.cfg, both are located in the directory config. Configuration file plot.cfg have only one section main that contain matplotlib configuration see more in https://matplotlib.org/3.2.1/tutorials/introductory/customizing.html. Configuration file config.cfg have these sections paths, parameters, velocities, sources, cuts, base_frequencies_SDR, base_frequencies_DBBC, stations, gauss_lines, Full_source_name. The paths sections contain all of the data input and output paths. The parameters section are a collection of hardcoded parameters used in the data processing. The velocities section is to use to find the local maximum that is monitored. The section sources contain RA, DEC and epoch for observed source. The cuts section is signal regions that are used to compute signal to noise. The sections base_frequencies_SDR and base_frequencies_DBBC is used to compute the Dopler effect. The section stations contain stations coordinates. The section gauss_lines is used to compute Gauss approximation of spectre. Section Full_source_name is used for visualizing data.

## Directory structure
MDPS use 5 (**dataFilePath**, **logPath**, **outputFilePath**, **resultFilePath**, **prettyLogsPath**) different directories. The directory dataFilePath contains SDR output, directory logPath contain SDR logs, directory outputFilePath contain script _sdr_fs.py_ and script _total_spectrum_analyzer_qt5.py_ outputs, directory resultFilePath contains monitoring files, directory prettyLogsPath contains ExperimentsLogReader output products.

## Processing SDR output
SDR for each scan creates four files **r0** **r1** **s0** **s1**. File name is &lt;source&gt; __f&lt;frequency&gt; _&lt;station label&gt; _&lt;iteration&gt; _no&lt;scan number&gt;&lt;r0, r1, s0, s1&gt;.dat file type is ASCII. 
Script _sdr_fs.py_ us frequency shifting algorithm described by publication: 
* Winkel, B., Kraus, A. and Bach, U., 2012. Unbiased flux calibration methods for spectral-line radio observations. Astronomy & Astrophysics, 540, p.A140. (https://www.aanda.org/articles/aa/pdf/2012/04/aa18092-11.pdf)

VIRAC SDR back-end is described in this publication:
* M. Bleiders et al., Spectral Line Registration Back-end based on USRP X300 Software Defined Radio, Journal of Astronomical Instrumentation, doi: 10.1142/S22511717205000099, 2020.
(https://www.worldscientific.com/doi/abs/10.1142/S2251171720500099)

It creates output file with name <source> _<MJD> _<station name> _<iteration>.h5 that file, has hdf5 format. With table amplitude, that have colons velocity and amplitude for left and right polarization.

Script total_spectrum_analyzer_qt5.py smooth data and create output for monitoring. It appends sdr_fs.py output with these tables amplitude_corrected, amplitude_corrected_not_smooht.

The main.py script has source and line mandatory parameters. The main.py script can be run with additional options:

-v or --version to display the current version

-h or --help for help

-c or --config to point to configuration file. Default path is: config/config.cfg

| **Scripts** | **Description** |
| --- | --- |
| main.py | Automatically call sdr_fs.py and total_spectrum_analyzer_qt5.py |
| sdr_fs.py | Process four output files from SDR |
| total_spectrum_analyzer_qt5.py | Process sdr_fs.py output|

## Monitoring

To view monitoring script _monitoring.py_ must be executed, this script has no parameters. After executing this script users will see GUI form that will ask to set observed source and frequency, after user will see flux density (Jy) of average polarization as function of time. By clicking on each point user will see spectre of that observation. User have options to change polarizations from average to left or right. User also two other options _plot periods_ and _plot maps_. The option _plot periods_ show _Lomb-Scargle Periodogram_ using https://docs.astropy.org/en/stable/api/astropy.timeseries.LombScargle.html. The option _plot maps_ create map of monitoring x axis is velocity, y axis is time color is flux density (Jy).   
An additional option is to use script _multiple_line_monitoring.py_, this script show monitoring of multiple stable sources, currently these source and thay compunents are selected - g32p745 (30.49, 39.18), w51 (59.29), g59p783 (19.2), on1 (14.64), s252 (10.84), ngc7538 (-58.04), w3oh (-44.6).

| **Scripts** | **Description** |
| --- | --- |
| monitoring.py | Show how amplitude change in time for given source|
| multiple_line_monitoring.py | Show how amplitude change in time for calibrator source|

## Publication and presentation of observations

This package contains various scripts that allow creating visualisations for presentations and publications, these scripts should be executed from publications_and_presentations folder.


| **Scripts** | **Description** |
| --- | --- |
| maser_hist.py | Has parameter _spectr_files_for_all_sources.py_, that contain list for each source, in that list are three output file names(for observation with lowest, middle and high amplitude) for all the source spectral low middle and high density is computed and histogram computed.|
| monitoring_for_publication.py | Create LATEX table for source with variability statistics|
| spectr_monitoring.py | Create a plot for publications where spectra and monitoring plot are viewed side by side|
| spectr_movie.py | Show how spectr changes over time and allow to create movie|

## Changelog

Version 1.0

First release. Changes since version pre1.0 are:
- Enable SDR data processing
- Enable to view monitoring
- Enable to create visualisations for presentations and publications

## Getting Help

Bug reports, feature requests and make contributions (e.g. code patches) can be reported by opening a &quot;new issue&quot; ticket on GitHub. Please give as much information (e.g. the software component, version) as you can in the ticket. For bugs, it is extremely useful if a small self-contained code snippet that reproduces the problem is provided.

## Acknowledgements
This software was written by Jānis Šteinbergs under the supervision of Artis Aberfelds. If you make use of this software to get results that appear in a publication or presentation please include this acknowledgement: &quot;We have made use of MDPS, a tool developed by Janis Steinbergs.&quot;

This work is the result of project implementation: «Physical and chemical processes in the interstellar medium», No 1.1.1.1/16/A/213 supported by ERDF​.

This work was supported by Latvian Council of Science Project “Research of Galactic Masers” Nr.: lzp-2018/1-0291.
