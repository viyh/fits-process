# fits-process
# https://github.com/viyh/fits-process

import pkg_resources
import pyfits
import sys
import csv
import os
import errno
import glob
import argparse

print "fits-process v0.0.1\n"
print "Using pyfits version [%s]\n" % pkg_resources.get_distribution('pyfits').version

log = {}

# -----------------------------------------------------------------------------

def parse_args():
    '''
    Parse the command line arguemnts.
    '''
    parser = argparse.ArgumentParser(
        description = 'Slice fits cubes and record header info.' )
    parser.add_argument('-f', '--files', required = True, help = 'Search string for the fits files you want to ingest.  e.g. "data_dir/*.fits"')
    parser.add_argument('-l', '--logfile', required = True, help = 'Where to write the log file (CSV format)')
    parser.add_argument('-o', '--out', required = True, help = 'Directory for output')
    parser.add_argument('-s', '--slice', action='store_true', help = 'Slice fits into separate files')
    parser.add_argument('-D', '--debug', action='store_true', help = 'Debugging')
    args = parser.parse_args()
    return args

# -----------------------------------------------------------------------------

def get_header(filename):
    '''
    Get the header information.
    '''
    f = pyfits.open(filename)
    header = f[0].header
    f.close()
    header.update('filename', os.path.basename(filename))

    return header

# -----------------------------------------------------------------------------

def ingest_header(args,filename):
    '''
    Ingest the header information and add it to the log.
    '''
    if args.debug == True:
        print filename
    header = get_header(filename)
    log[filename] = header
    return True

# -----------------------------------------------------------------------------

def ingest_data(args,filename):
    '''
    Ingest the data from the FITS file
    '''
    data_cube, header_data_cube = pyfits.getdata(filename, 0, header=True)
    numslices = data_cube.shape[0]

    imfiledir = os.path.join( args.out, os.path.splitext(os.path.split(filename)[1])[0] )

    mkdir_p( imfiledir )

    #import pylab
    #import mpl_toolkits.mplot3d.axes3d as p3

    #for i in range(numslices):
    for i in range(10):
        imfileout = os.path.join( imfiledir, os.path.splitext(os.path.split(filename)[1])[0] + '_' + ( "%05d" % i ) + '.fits' )
        hdu = pyfits.PrimaryHDU( data_cube[i] )
        hdu.writeto(imfileout)
        #pylab.plot(data_cube[i], height=32767)
        #pylab.imshow(data_cube[i], interpolation="nearest")
        #pylab.contour(data_cube[i])
        #pylab.show()

    return True

# -----------------------------------------------------------------------------

def mkdir_p(path):
    try:
        os.makedirs(path)
    except OSError as exc:
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else: raise

# -----------------------------------------------------------------------------

def get_file_list(search_string):
    '''
    Get the file list based on the seach string.
    '''
    search_string = os.path.abspath(search_string)
    print search_string
    file_list = glob.glob(search_string)

    assert file_list != [], 'No files found.'
    return file_list

# -----------------------------------------------------------------------------

def write_log(args,logfile):
    '''
    Write the header info for each file to a CSV output log
    '''
    with open(logfile, 'wb') as csvfile:
        logwriter = csv.writer(csvfile, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
        logwriter.writerow( log[log.keys()[0]].keys() )
        for logrow in log.keys():
            logwriter.writerow(log[logrow])
    return True

# -----------------------------------------------------------------------------

def main(args):
    '''
    The main controller.
    '''
    file_list = get_file_list(args.files)

    for filename in file_list:
        ingest_header(args,filename)
        ingest_data(args,filename)

    write_log(args,args.logfile)

# -----------------------------------------------------------------------------

if __name__ == '__main__':
    args = parse_args()
    main(args)

