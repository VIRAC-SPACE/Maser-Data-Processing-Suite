# Utils

| **Scripts** | **Description** |
| --- | --- |
| convert_hdf5_to_ascii.py | Convert output file (format hdf5) to text file (ASCII format) script has three parameters input file (hdf5), output file (ASCII) and table what table to store.. |
| delete_flag_observation.py | Delete flagged observations. Search for flagged observation in result file for source, has two parameters source and frequency. |
| fix_result_amplitudes.py | If the output file is changed manually, run this script to change result file, has three parameters source, frequency and back end type. |
| observation_correction.py | Correct observation by a given factor, has six parameters source, frequency, factor, station, back end type and iteration list. This script will multiply observations from iteration list with factor. |
| help.py | Common used functions. |
| ploting_qt5.py | Plotting class to embed matplotlib to pyqt5. |
| vlsr.py | Compute local standard of rest. |