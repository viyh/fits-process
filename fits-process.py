# fits-process
# https://github.com/viyh/fits-process

import pkg_resources, sys, os, errno, glob, argparse, csv
import pyfits, numpy, scipy, scipy.ndimage
import matplotlib.pyplot as plt
from math import *

print "fits-process v0.0.3\n"
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
    parser.add_argument('-l', '--logfile', required = False, help = 'Where to write the log file (CSV format)')
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
    parser.add_argument('-C', '--camangle', required = False, type=float, help = 'Calibration camera angle value (degrees)')
    parser.add_argument('-s', '--slice', action='store_true', help = 'Slice fits file into separate files')
    parser.add_argument('-c', '--chop', required=False, type=int, help = 'Chop fits file into separate groups of this many slices each')
    # dark subtraction not fully implemented yet
    parser.add_argument('-d', '--dark', required=False, help = 'Dark calibration file to subtract from data')
    parser.add_argument('-b', '--bsize', required=False, default = 10, type=int, help = 'Box size')
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
        print "Found file [%s]" % filename
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

def ingest_data(args, log, filename):
    '''
    Ingest the data from the FITS file
    '''
    data_cube, header = pyfits.getdata(filename, 0, header=True)
    numslices = data_cube.shape[0]

    dark = numpy.zeros_like(data_cube, dtype=None, order='K', subok=True)
    dark = dark[0]

    if args.out and not args.chop:
        split_cube(args,data_cube,filename)

    if args.dark:
        dark_cube, dark_exp = ingest_dark(args)
        dark_cube *= (header['EXPOSURE'] / dark_exp)
        dark = numpy.mean(dark_cube, axis=0)

    if args.chop:
        # figure out the number of pieces to chop the fits file into
        avg_r = 0.0
        avg_t = 0.0
        pieces = numslices / args.chop
        if pieces < 1: pieces = 1
        for p in range(pieces):
            new_data_cube = numpy.split(data_cube, pieces)[p]
            start = args.chop * p
            end = args.chop * ( p + 1 ) - 1
            if args.out:
                chop_cube(args,new_data_cube,filename,start,end)
            (r, t, a, b) = align(args,new_data_cube,dark,os.path.splitext(os.path.split(filename)[1])[0] + ' (slice ' + str(start) + ' to ' + str(end) + ')' )
            if p == 0:
                avg_r = r
                avg_t = t
            else:
                avg_r = (avg_r + r) / 2
                avg_t = (avg_t + t) / 2

            log[filename]['raw_rho'] = str(avg_r)
            log[filename]['raw_theta'] = str(avg_t)
            log[filename]['adj_rho'] = str(avg_r * args.platescale)
            log[filename]['adj_theta'] = str( 270 - (args.camangle + avg_t) )

        print "\nAverage rho:\t%.4f" % avg_r
        print "Average theta:\t%.4f\n" % avg_t
        if args.platescale:
            print "\nAdjusted average rho:\t%.4f" % (avg_r * args.platescale)
        if args.camangle:
            print "Adjusted average theta:\t%.4f\n" % ( 270 - (args.camangle + avg_t) )

    else:
        (raw_rho, raw_theta, adj_rho, adj_theta) = align(args,data_cube,dark,os.path.splitext(os.path.split(filename)[1])[0])
        log[filename]['rawrho'] = "%.4f" % raw_rho
        log[filename]['rawtheta'] = "%.4f" % raw_theta
        log[filename]['adjrho'] = "%.4f" % adj_rho
        log[filename]['adjtheta'] = "%.4f" % adj_theta


    return True

# -----------------------------------------------------------------------------

def chop_cube(args,data_cube,filename,start,end):
    '''
    Split up the slices of a FITS cube into sets
    '''
    imfiledir = os.path.join( args.out, os.path.splitext(os.path.split(filename)[1])[0] )
    mkdir_p( imfiledir )

    imfileout = os.path.join( imfiledir, os.path.splitext(os.path.split(filename)[1])[0] + '_' + str(start) + '-' + str(end) + '.fits' )
    if args.overwrite and (os.path.exists(imfileout) or os.path.islink(imfileout)):
        os.remove(imfileout)
    hdu = pyfits.PrimaryHDU( data_cube )
    hdu.writeto(imfileout)
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

def ingest_dark(args):
    dark_cube, dark_header = pyfits.getdata(args.dark, 0, header=True)
    dark_exp = dark_header['EXPOSURE']
    return dark_cube, dark_exp

# -----------------------------------------------------------------------------

def align(args, data_cube, dark, title):
    # Set the normalization range
    norm = 2**int(args.normalize) - 1

    # Grab first slice from the cube
    slice_ref = data_cube[0] - dark

    mult = 0

    # attempt to automatically figure out the proper threshold
    for i in range(0, 50):
        #slice_ref *= norm / numpy.amax(slice_ref, axis=None, out=None)
        mult = .2 * i

        # Set the threshold to the standard deviation of the data in the
        # first slice multiplied by the supplied argument
        threshold_ref = slice_ref.std( dtype=numpy.float32 ) * (args.threshold + mult)
        if threshold_ref > norm: threshold_ref = norm
        if args.debug:
            print "threshold: ", threshold_ref, "\nstddev of first slice: ", slice_ref.std( dtype=numpy.float32 )

        # Find useful data that exceeds the threshold value to determine the COM of the image
        lbl_ref, num_ref = scipy.ndimage.measurements.label( slice_ref >= threshold_ref, numpy.ones((3,3)) )
        com_ref = scipy.ndimage.measurements.center_of_mass( slice_ref, lbl_ref, range(1, num_ref + 1) )
        if args.debug:
            #print "Centers_ref:\n", com_ref
            print "Num_ref: ", num_ref

        # This is the center of the reference image
        x_ref_tmp = numpy.array(com_ref)[:,0]
        y_ref_tmp = numpy.array(com_ref)[:,1]
        if args.debug:
            print "x =\n",x_ref_tmp[0],"\ny =\n",y_ref_tmp[0]

        # test to see if the centers of mass are anywhere near the center of the image
        if x_ref_tmp[0] < (slice_ref.shape[1] / 8) or x_ref_tmp[0] > slice_ref.shape[1] - (slice_ref.shape[1] / 8):
            continue
        elif y_ref_tmp[0] < (slice_ref.shape[0] / 8) or y_ref_tmp[0] > slice_ref.shape[0] - (slice_ref.shape[0] / 8):
            continue
        else:
            # we've found a reasonable threshold, so stop
            break

    x_ref = slice_ref.shape[1] / 2
    y_ref = slice_ref.shape[0] / 2

    # copy data cube to throw processed slices into
    stacked_data_cube = numpy.zeros_like(data_cube, dtype=None, order='K', subok=True)

    for i in range( data_cube.shape[0] ):
        # Grab the next slice to process
        slice = data_cube[i] - dark

        # Run a multidimensional Laplacian filter
        slice = scipy.ndimage.filters.gaussian_laplace( \
                    data_cube[i], args.sigma, output=None, mode='reflect', cval=0.0 )

        # Set the threshold to the standard deviation of the data in the
        # first slice multiplied by the supplied threshold argument
        threshold = slice.std( dtype=numpy.float32 ) * (args.threshold + mult)
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
        if args.debug: print "Image offset: ", i, "\t", round( x_ref-x[0],2), round(y_ref-y[0], 2 )
        sys.stdout.write(".")
        sys.stdout.flush()

        # Normalize slice data
        #slice = slice * ( norm / numpy.amax(slice, axis=None, out=None) )

        # align center of mass for slice to reference image and add to stacked data cube
        stacked_data_cube[i] = scipy.ndimage.interpolation.shift(slice,[x_ref-x[0],y_ref-y[0]])

    # This is the fully reduced image.
    stacked_image = numpy.mean(stacked_data_cube, axis=0)

    # An attempt to determine primary and secondary stars. Kind of a hack, should be worked on (k-means?) 
    idx1 = numpy.unravel_index( stacked_image.argmax(), stacked_image.shape )
    print "\nIndex of primary star:\t\t", idx1
    stacked_image_tmp = stacked_image.copy()
    bsize = args.bsize
    for i in range(bsize):
        for j in range(bsize):
            stacked_image_tmp[ (idx1[0]-(bsize/2))+i, (idx1[1]-(bsize/2))+j ] = 0
    idx2 = numpy.unravel_index( stacked_image_tmp.argmax(), stacked_image_tmp.shape )
    print "Index of secondary star:\t", idx2

    print "Delta magnitude:\t\t%.5f" % (2.5 * log10( numpy.amax(stacked_image) / numpy.amax(stacked_image_tmp) ) )

    # Figure out the relative pixel distance and theta
    x = idx2[1] - idx1[1]
    y = idx2[0] - idx1[0]
    raw_theta = degrees( atan2( y, x ) ) % 360
    raw_rho = sqrt( x**2 + y**2 )
    print "\ntheta:\t%f degrees" % raw_theta
    print "r:\t%f pixels" % raw_rho

    # If calibration offsets are given, output an adjusted rho/theta value
    adj_rho = 0.0
    adj_theta = 0.0
    if args.platescale:
        adj_rho = args.platescale * raw_rho
        print "\nAdjusted separation:\t\t%f" % adj_rho
    if args.camangle:
        # Need to fix this!
        adj_theta = 270 - (args.camangle + raw_theta)
        if adj_theta < 0:
            adj_theta += 360.0
        print "Adjusted position angle:\t%f" % adj_theta

    if args.images:
        plot_data(slice_ref, stacked_image, title, norm, x_ref, y_ref, idx1, idx2)

    return raw_rho, raw_theta, adj_rho, adj_theta

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
    plt.axis([0.0, stacked_image.shape[1], 0.0, stacked_image.shape[0]])

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
    file_list = get_file_list( args.files )
    file_list.sort()

    if not args.out:
        print "No output directory specified. Running in read-only mode."

    for filename in file_list:
        ingest_header(args,log,filename)
        ingest_data(args,log,filename)

    if args.logfile:
        write_log(args, log)

# -----------------------------------------------------------------------------

if __name__ == '__main__':
    args = parse_args()
    main(args)

