fits-process
============

Process a set of astronomical FITS multidimensional data cubes into separate 
slices, and record the header data into a log file. Aligns double star images
and attempts to guess the position angle and separation.

An example of what this script can do:

http://www.lowmagnitude.com/data/reduced_example.jpg

The first two plots are a slice of the original data. The last four are different
types of plots of the data after it has been reduced by this script.

A sample FITS data cube can be downloaded here:

http://www.lowmagnitude.com/data/001_BU533-100.fits.zip

joe-at-lowmagnitude.com

Usage
-----

usage: fits-process.py [-h] -f FILES [-o OUT] [-l LOGFILE] [-O] [-n NORMALIZE]
                       [-S SIGMA] [-t THRESHOLD] [-i] [-P PLATESCALE]
                       [-T THETA] [-s] [-c CHOP] [-b BSIZE] [-D]

    Slice fits cubes, align slices, and record header info.

    optional arguments:
      -h, --help            show this help message and exit
      -f FILES, --files FILES
                            Search string for the fits files you want to ingest.
                            e.g. "data_dir/*.fits"
      -o OUT, --out OUT     Directory for output files
      -l LOGFILE, --logfile LOGFILE
                            Where to write the log file (CSV format)
      -O, --overwrite       Overwrite existing output files
      -n NORMALIZE, --normalize NORMALIZE
                            Normalization range (2^n, where n is this value,
                            default = 8)
      -S SIGMA, --sigma SIGMA
                            Sigma value for Laplace filter (default = 2.0)
      -t THRESHOLD, --threshold THRESHOLD
                            Threshold of usable data values, only data points
                            above this will be considered (Multiplied by the
                            standard deviation of the data per slice, default =
                            4.0)
      -i, --images          Display images
      -P PLATESCALE, --platescale PLATESCALE
                            Calibration plate scale value (arc-seconds per pixel)
      -C CAMANGLE, --camangle CAMANGLE
                            Calibration camera angle value (degrees)
      -s, --slice           Slice fits file into separate files
      -c CHOP, --chop CHOP  Chop fits file into separate groups of this many
                            slices each
      -b BSIZE, --bsize BSIZE
                            Box size
      -D, --debug           Debugging
    
Example:

To align images and calculate a relative theta and pixel distance, and then show the image plots:

`python fits-process.py --logfile output.csv --file '001_BU533-100.fits' -t 4 -S 2 -C 0.0538 -T -22.46 -i`

A sample FITS cube to play with is linked above. The camera angle is -22.46 degrees and the
platescale was 0.0538 arcseconds per pixel.

It is usually good to start with the threshold set fairly low (around 2) and the sigma set
low as well (around 1). Turn these up as needed. If the image is not getting centered properly,
turning up the threshold will help find the center of the reference image better. If the stars
are not showing up well, then turning up the sigma value will help, but decrease the precision
of the results. These two numbers are floating point, so they can be turned up slowly. I usually
increment by one until I have a desired result, then fine tune as needed.

If files are named in the format "ID_Target.fits" where "ID" is some arbitrary ID number and 
"Target" is the name of the target star, the code will automatically add this to the CVS log 
file and output of the graphs.

Basic Functionality
-------------------

* Grab filenames of specified FITS files
* Parse header information
* Use filename to create directory in target data output directory
* Slice up the data into separate images
* Save header informations to CVS log (todo: add more output types)
* Alignment of slices and analysis

To Do
-----

* Fix the rho/theta calibration stuff
* Add a "calibration" mode that will calculate the plate scale and camera angle
* Write out the calculated rho/theta value to the log
* Make bsize a CLI option
* Clean up function names

Dependencies
------------

* PyFITS [http://www.stsci.edu/institute/software_hardware/pyfits]
* NumPy [http://numpy.scipy.org/]
* SciPy [http://www.scipy.org/]
* Matplotlib [http://matplotlib.org/]

References
----------

* Homepage [https://github.com/viyh/fits-process]
* NumPy docs [http://docs.scipy.org/doc/numpy/reference/arrays.ndarray.html]
* AstroPython [http://www.astropython.org/]
* PyFITS docs [http://packages.python.org/pyfits/]
* Matplotlib docs [http://matplotlib.org/contents.html]
* SciPy docs [http://docs.scipy.org/doc/scipy/reference/ndimage.html]
