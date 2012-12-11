# fits-process
# https://github.com/viyh/fits-process

import pkg_resources, sys, os, errno, glob, argparse, csv
import pyfits, numpy, scipy, scipy.ndimage
import matplotlib.pyplot as plt
from math import *

print "fits-process v0.0.2\n"
print "Using pyfits version [%s]\n" % pkg_resources.get_distribution('pyfits').version

# -----------------------------------------------------------------------------

def parse_args():
    '''
    Parse the command line arguments.
    '''
    parser = argparse.ArgumentParser(description = 'Slice fits cubes, align slices, and record header info.' )
    parser.add_argument('-f', '--files', required = True, help = 'Search string for the fits files you want to \
                         ingest.  e.g. "data_dir/*.fits"')
    parser.add_argument('-o', '--out', required = False, help = 'Directory for output files')
    parser.add_argument('-l', '--logfile', required = True, help = 'Where to write the log file (CSV format)')
    parser.add_argument('-O', '--overwrite', action='store_true', help = 'Overwrite existing output files')
    parser.add_argument('-n', '--normalize', required = False, type=int, default = 8, help = 'Normalization \
                         range (2^n, where n is this value, default = 8)')
    parser.add_argument('-S', '--sigma', required = False, type=float, default = 2.0, help = 'Sigma value for \
                         Laplace filter (default = 2.0)')
    parser.add_argument('-t', '--threshold', required = False, type=float, default = 4.0, help = 'Threshold of \
                         usable data values, only data points above this will be considered (Multiplied by the \
                         standard deviation of the data per slice, default = 4.0)')
    parser.add_argument('-i', '--images', action='store_true', help = 'Display images')
    parser.add_argument('-P', '--platescale', required = False, type=float, help = 'Calibration plate scale \
                         value (arc-seconds per pixel)')
    parser.add_argument('-T', '--theta', required = False, type=float, help = 'Calibration theta value (degrees)')
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

def ingest_header(args, log, filename):
    '''
    Ingest the header information and add it to the log.
    '''
    if args.debug == True:
        print "Found file [%f]" % filename
    header = get_header(filename)

    print filename
    log[filename] = header

    # For filenames that are in the format "IDNumber_TargetName", parse and add these to the log
    import re
    m = re.match(r'.+?(\d+)_(.+)\.fits', filename)
    if m and m.group():
        log[filename]['ID'] = m.group(1)
        log[filename]['TARGET'] = m.group(2)

    return True

# -----------------------------------------------------------------------------

def ingest_data(args, filename):
    '''
    Ingest the data from the FITS file
    '''
    data_cube, header_data_cube = pyfits.getdata(filename, 0, header=True)
    numslices = data_cube.shape[0]

    if args.out:
        split_cube(args,data_cube,filename)

    align(args,data_cube,os.path.splitext(os.path.split(filename)[1])[0])

    return True

# -----------------------------------------------------------------------------

def split_cube(args, data_cube, filename):
    '''
    Split up the slices of a FITS cube into individual files
    '''
    imfiledir = os.path.join( args.out, os.path.splitext(os.path.split(filename)[1])[0] )
    mkdir_p( imfiledir )

    for i in range(numslices):
        imfileout = os.path.join( imfiledir, os.path.splitext(os.path.split(filename)[1])[0] + \
            '_' + ( "%05d" % i ) + '.fits' )
        if args.overwrite and (os.path.exists(imfileout) or os.path.islink(imfileout)):
            os.remove(imfileout)
        hdu = pyfits.PrimaryHDU( data_cube[i] )
        if args.out:
            hdu.writeto(imfileout)
    return True

# -----------------------------------------------------------------------------

def align(args, data_cube, title):

    # Set the normalization range
    norm = 2**int(args.normalize) - 1

    # Grab first slice from the cube
    slice_ref = data_cube[0]
    slice_ref *= norm / numpy.amax(slice_ref, axis=None, out=None)

    # Set the threshold to the standard deviation of the data in the
    # first slice multiplied by the supplied argument
    threshold_ref = args.threshold * slice_ref.std( dtype=numpy.float32 )
    if threshold_ref > norm: threshold_ref = norm
    print "threshold: ", threshold_ref, "\nstddev of first slice: ", slice_ref.std( dtype=numpy.float32 )

    # Find useful data that exceeds the threshold value to determine the COM of the image
    lbl_ref, num_ref = scipy.ndimage.measurements.label( slice_ref >= threshold_ref, numpy.ones((3,3)) )
    com_ref = scipy.ndimage.measurements.center_of_mass( slice_ref, lbl_ref, range(1, num_ref + 1) )
    print "Centers_ref:\n", com_ref
    print "Num_ref: ", num_ref

    # This is the center of the reference image
    x_ref = numpy.array(com_ref)[:,0]
    y_ref = numpy.array(com_ref)[:,1]
    print "x =\n",x_ref,"\ny =\n",y_ref

    # copy data cube to throw processed slices into
    stacked_data_cube = numpy.zeros_like(data_cube, dtype=None, order='K', subok=True)

    for i in range( data_cube.shape[0] ):
        # Grab the next slice to process
        slice = data_cube[i]

        # Run a multidimensional Laplacian filter
        slice = scipy.ndimage.filters.gaussian_laplace( \
                    data_cube[i], args.sigma, output=None, mode='reflect', cval=0.0 )
        threshold = args.threshold * slice.std( dtype=numpy.float32 )
        if threshold > norm: threshold = norm

        # Find useful data that exceeds the threshold value to determine the COM of the image
        lbl, num = scipy.ndimage.measurements.label( slice >= threshold, numpy.ones((3,3)) )
        com = scipy.ndimage.measurements.center_of_mass( slice, lbl, range(1,num+1) )

        # If a center of mass is found, set x and y values to it.
        # If not, skip this slice (try setting the threshold lower if this happens)
        if com:
            x = numpy.array(com)[:,0]
            y = numpy.array(com)[:,1]
        else: continue
        if args.debug: print "Image offset: ", i, "\t", round( x_ref[0]-x[0],2), round(y_ref[0]-y[0], 2 )
        sys.stdout.write(".")
        sys.stdout.flush()

        # Normalize slice data
        slice = slice * ( norm / numpy.amax(slice, axis=None, out=None) )

        # align center of mass for slice to reference image and add to stacked data cube
        stacked_data_cube[i] = scipy.ndimage.interpolation.shift(slice,[x_ref[0]-x[0],y_ref[0]-y[0]])

    # This is the fully reduced image.
    stacked_image = numpy.mean(stacked_data_cube, axis=0)

    # An attempt to determine primary and secondary stars. Kind of a hack, should be worked on (k-means?) 
    idx1 = numpy.unravel_index( stacked_image.argmax(), stacked_image.shape )
    print "\nIndex of primary star:\t\t", idx1
    stacked_image_tmp = stacked_image.copy()
    bsize = 10
    for i in range(bsize):
        for j in range(bsize):
            stacked_image_tmp[ (idx1[0]-(bsize/2))+i, (idx1[1]-(bsize/2))+j ] = 0
    idx2 = numpy.unravel_index( stacked_image_tmp.argmax(), stacked_image_tmp.shape )
    print "Index of secondary star:\t", idx2

    # Figure out the relative pixel distance and theta
    x = idx2[1] - idx1[1]
    y = idx2[0] - idx1[0]
    theta = degrees( atan2( y, x ) ) % 360
    r = sqrt( x**2 + y**2 )
    print "\ntheta:\t%f degrees" % theta
    print "r:\t%f pixels" % r

    # If calibration offsets are given, output an adjusted rho/theta value
    if args.platescale:
        sep = args.platescale * r
        print "\nAdjusted separation:\t\t%f" % sep
    if args.theta:
        # Need to fix this!
        pa = (args.theta + 90 + theta)
        if pa < 0:
            pa += 360.0
        print "Adjusted position angle:\t%f" % pa

    if args.images:
        plot_data(slice_ref, stacked_image, title, norm, x_ref, y_ref, idx1, idx2)

    return True

# -----------------------------------------------------------------------------

def plot_data(slice_ref, stacked_image, title, norm, x_ref, y_ref, idx1, idx2):
    '''
    Graph the data
    '''
    plt.rcParams.update({'font.size': 9})

    plt.subplot(2,3,1), plt.imshow(slice_ref, origin='lower')
    plt.title(title + " first slice")

    plt.subplot(2,3,2), plt.imshow(slice_ref, cmap="gray", origin='lower')
    plt.title(title + " first slice grayscale")

    # Plots a green square marker over guess for the primary star, 
    # and a red square marker over guess for the secondary star
    plt.subplot(2,3,3), plt.imshow(stacked_image, norm=plt.Normalize(0,norm), origin='lower', vmin=0)
    plt.plot(idx1[1], idx1[0],'gs',markersize=7)
    plt.plot(idx2[1], idx2[0],'rs',markersize=7)
    plt.plot(x_ref,y_ref,'ro',alpha=0.2, markersize=3)
    plt.title(title + " normalized over 0-" + str( norm ) )
    plt.axis([0.0, 128.0, 0.0, 128.0])

    plt.subplot(2,3,4), plt.imshow(stacked_image, cmap="gray", norm=plt.Normalize(0,norm), origin='lower')
    plt.set_cmap('spectral')
    plt.title(title + " spectral grayscale normalized over 0-" + str( norm ) )

    plt.subplot(2,3,5), plt.imshow(stacked_image, interpolation="nearest", cmap="spectral", origin='lower')
    plt.title(title + " interpolated spectral")

    plt.subplot(2,3,6), plt.imshow(stacked_image, cmap="gray", norm=plt.Normalize(0,norm), origin='lower')
    plt.title(title + " grayscale normalized over 0-" + str( norm ) )

    plt.show()

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
    print "Finding files [%s]...." % search_string
    file_list = glob.glob(search_string)

    assert file_list != [], 'No files found.'
    return file_list

# -----------------------------------------------------------------------------

def write_log(args, log):
    '''
    Write the header info for each file to a CSV output log
    '''
    with open(args.logfile, 'wb') as csvfile:
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
    log = {}
    file_list = get_file_list(args.files)

    if not args.out:
        print "No output directory specified. Running in read-only mode."

    for filename in file_list:
        ingest_header(args,log,filename)
        ingest_data(args,filename)

    write_log(args, log)

# -----------------------------------------------------------------------------

if __name__ == '__main__':
    args = parse_args()
    main(args)

