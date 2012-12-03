fits-process
============

Process a set of astronomical FITS multidimensional data cubes into separate 
slices, and record the header data into a spreadsheet.

Usage
-----

    usage: fits-process.py [-h] -f FILES -l LOGFILE -o OUT [-s] [-D]

    Slice fits cubes and record header info.

    optional arguments:
      -h, --help            show this help message and exit
      -f FILES, --files FILES
                            Search string for the fits files you want to ingest.
                            e.g. "data_dir/*.fits"
     -l LOGFILE, --logfile LOGFILE
                            Where to write the log file (CSV format)
      -o OUT, --out OUT     Directory for output
      -s, --slice           Slice fits into separate files
      -D, --debug           Debugging

Example:
`python fits-process.py -s -f '/home/joe/data/*.fits' -l /home/joe/test.csv -o /home/joe/output/`

Basic Functionality
-------------------

* Grab filenames of specified FITS files
* Parse header information
* Use filename to create directory in target data output directory
* Slice up the data into separate images
* Save header informations to CVC log (add more output types)
* Add data analysis

Dependencies
------------

* PyFITS http://www.stsci.edu/institute/software_hardware/pyfits]

References
----------

joe-at-lowmagnitude.com

* Homepage [https://github.com/viyh/fits-process]
* NumPy docs [http://docs.scipy.org/doc/numpy/reference/arrays.ndarray.html]
* AstroPython [http://www.astropython.org/]
* PyFITS docs [http://packages.python.org/pyfits/]
