#!/usr/bin/env python

# Library for treating FITS files as Python Imaging Library objects
# Copyright (c) 2005, 2006, 2007, Jeremy Brewer
# 
# All rights reserved.
# 
# Redistribution and use in source and binary forms, with or without 
# modification, are permitted provided that the following conditions are 
# met:
# 
#     * Redistributions of source code must retain the above copyright 
#       notice, this list of conditions and the following disclaimer.
#     * Redistributions in binary form must reproduce the above copyright 
#       notice, this list of conditions and the following disclaimer in 
#       the documentation and/or other materials provided with the
#       distribution.
#     * The names of the contributors may not be used to endorse or 
#       promote products derived from this software without specific prior 
#       written permission.
# 
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
# "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
# LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
# A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT 
# OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
# EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
# PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
# PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
# LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
# NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#
# Changelog:
#
# 9/25/07   Added call to fits_simple_verify() to verify input file is FITS.
#           Removed kwargs from FitsImage() because pyfits doesn't use them.
#
# 9/14/07   Change array usage from Numeric to numpy.  Changed underlying
#           FITS I/O library from fitslib to pyfits.  Modifications made
#           by Christopher Hanley.
#
# 8/20/07   Write arcsinh scaling algorithm and adding scaling options.
#           Updated documentation.  Dropped number of channels check on
#           color -- PIL should handle this instead.
#
# 8/17/07   Wrote new scaling algorithm, percentile_range(), that determines
#           the range to use from a configurable percentile cut.  Now
#           FitsImage() takes optional named arguments to configure which
#           contrast algorithm to use.  In addition, keyword arguments are
#           passed on to Fits() to configure how minor errors are handled.
#
# 7/4/07    Updated to use Numeric.  Improved speed of zscale_range().
#
# 10/10/06  Increased accuracy of draw_circle().
#
# 2/7/06    Updated documentation.
#
# 1/4/06    Fixed bug in zscale_range() where num_points and num_pixels
#           sometimes differed, resulting in the sigma iteration failing because
#           the arrays would differ in length.  Now the arrays are both of
#           size num_pixels.  Some additional checks were also added.
#
# 12/10/05  Updated documentation.
#
# 12/8/05   Now draw_circle will not draw points that lie outside of the image.
#
# 12/7/05   Wrote zscale_range() function which implements the ds9 zscale
#           autocontrast algorithm for FITs images.  Wrote a new version of
#           asImage(), now called FitsImage(), that returns a PIL Image object
#           without use of the convert commandline utility.  Rewrote convert()
#           and resize() methods so that they do not have to use the convert
#           command externally.  Removed all of the other asImage() methods
#           that weren't working.

"""
Module for treating a FITS image as a Python Imaging Library (PIL) object.

This is extremely useful if you want to convert a FITS image to jpeg and/or
perform various operations on it (such as drawing).

The contrast for FITS images is determined using the zscale algorithm by
default, but this can be configured with various options.  See the
documentation for zscale_range() and percentile_range() for more information.

Example Usage:

# image is a full PIL object
image = fitsimage.FitsImage("foo.fits")
image.save("foo.jpg")
"""

__author__ = "Jeremy Brewer (jeremy.d.brewer@gmail.com)"
__copyright__ = "Copyright 2005, 2006, 2007 Jeremy Brewer"
__license__ = "BSD"
__version__ = "1.1"

import os
import sys
import cmath
import fitslib
import pyfits
import pointarray
import Image
import ImageDraw
import numpy

def zscale_range(image_data, contrast=0.25, num_points=600, num_per_row=120):

    """
    Computes the range of pixel values to use when adjusting the contrast
    of FITs images using the zscale algorithm.  The zscale algorithm
    originates in Iraf.  More information about it can be found in the help
    section for DISPLAY in Iraf.

    Briefly, the zscale algorithm uses an evenly distributed subsample of the
    input image instead of a full histogram.  The subsample is sorted by
    intensity and then fitted with an iterative least squares fit algorithm.
    The endpoints of this fit give the range of pixel values to use when
    adjusting the contrast.

    Input:  image_data  -- the array of data contained in the FITs image
                           (must have 2 dimensions)
            contrast    -- the contrast parameter for the zscale algorithm
            num_points  -- the number of points to use when sampling the
                           image data
            num_per_row -- number of points per row when sampling
    
    Return: 1.) The minimum pixel value to use when adjusting contrast
            2.) The maximum pixel value to use when adjusting contrast
    """

    # check input shape
    if len(image_data.shape) != 2:
        raise ValueError("input data is not an image")

    # check contrast
    if contrast <= 0.0:
        contrast = 1.0

    # check number of points to use is sane
    if num_points > numpy.size(image_data) or num_points < 0:
        num_points = 0.5 * numpy.size(image_data)

    # determine the number of points in each column
    num_per_col = int(float(num_points) / float(num_per_row) + 0.5)

    # integers that determine how to sample the control points
    xsize, ysize = image_data.shape
    row_skip = float(xsize - 1) / float(num_per_row - 1)
    col_skip = float(ysize - 1) / float(num_per_col - 1)

    # create a regular subsampled grid which includes the corners and edges,
    # indexing from 0 to xsize - 1, ysize - 1
    data = []
   
    for i in xrange(num_per_row):
        x = int(i * row_skip + 0.5)
        for j in xrange(num_per_col):
            y = int(j * col_skip + 0.5)
            data.append(image_data[x, y])

    # actual number of points selected
    num_pixels = len(data)

    # sort the data by intensity
    data.sort()

    # check for a flat distribution of pixels
    data_min = min(data)
    data_max = max(data)
    center_pixel = (num_pixels + 1) / 2
    
    if data_min == data_max:
        return data_min, data_max

    # compute the median
    if num_pixels % 2 == 0:
        median = data[center_pixel - 1]
    else:
        median = 0.5 * (data[center_pixel - 1] + data[center_pixel])

    # compute an iterative fit to intensity
    pixel_indeces = map(float, xrange(num_pixels))
    points = pointarray.PointArray(pixel_indeces, data, min_err=1.0e-4)
    fit = points.sigmaIterate()

    num_allowed = 0
    for pt in points.allowedPoints():
        num_allowed += 1

    if num_allowed < int(num_pixels / 2.0):
        return data_min, data_max

    # compute the limits
    z1 = median - (center_pixel - 1) * (fit.slope / contrast)
    z2 = median + (num_pixels - center_pixel) * (fit.slope / contrast)

    if z1 > data_min:
        zmin = z1
    else:
        zmin = data_min

    if z2 < data_max:
        zmax = z2
    else:
        zmax = data_max

    # last ditch sanity check
    if zmin >= zmax:
        zmin = data_min
        zmax = data_max

    return zmin, zmax

def percentile_range(image_data, min_percent=3.0, max_percent=99.0,
                     num_points=5000, num_per_row=250):
    
    """
    Computes the range of pixel values to use when adjusting the contrast
    of FITs images using a simple percentile cut.  For efficiency reasons,
    only a subsample of the input image data is used.

    Input:  image_data  -- the array of data contained in the FITs image
                           (must have 2 dimensions)
            min_percent -- min percent value between (0, 100)
            max_percent -- max percent value between (0, 100)
            num_points  -- the number of points to use when sampling the
                           image data
            num_per_row -- number of points per row when sampling

    Return: 1.) The minimum pixel value to use when adjusting contrast
            2.) The maximum pixel value to use when adjusting contrast
    """

    if not 0 <= min_percent <= 100:
        raise ValueError("invalid value for min percent '%s'" % min_percent)
    elif not 0 <= max_percent <= 100:
        raise ValueError("invalid value for max percent '%s'" % max_percent)

    min_percent = float(min_percent) / 100.0
    max_percent = float(max_percent) / 100.0

    # check input shape
    if len(image_data.shape) != 2:
        raise ValueError("input data is not an image")

    # check number of points to use is sane
    if num_points > numpy.size(image_data) or num_points < 0:
        num_points = 0.5 * numpy.size(image_data)

    # determine the number of points in each column
    num_per_col = int(float(num_points) / float(num_per_row) + 0.5)

    # integers that determine how to sample the control points
    xsize, ysize = image_data.shape
    row_skip = float(xsize - 1) / float(num_per_row - 1)
    col_skip = float(ysize - 1) / float(num_per_col - 1)

    # create a regular subsampled grid which includes the corners and edges,
    # indexing from 0 to xsize - 1, ysize - 1
    data = []
    
    for i in xrange(num_per_row):
        x = int(i * row_skip + 0.5)
        for j in xrange(num_per_col):
            y = int(j * col_skip + 0.5)
            data.append(image_data[x, y])

    # perform a simple percentile cut
    data.sort()
    zmin = data[int(min_percent * len(data))]
    zmax = data[int(max_percent * len(data))]
    
    return zmin, zmax

def FitsImage(fitsfile, contrast="zscale", contrast_opts={}, scale="linear",
              scale_opts={}):

    """
    Constructor-like function that returns a Python Imaging Library (PIL)
    Image object.  This allows extremely easy and powerful manipulation of
    FITs files as images.  The contrast is automatically adjusted using the
    zscale algorithm (see zscale_range() above).

    Input:  fitsfile      -- a FITS image filename
            contrast      -- the algorithm for determining the min/max
                             values in the FITS pixel data to use when
                             compressing the dynamic range of the FITS
                             data to something visible by the eye, either
                             "zscale" or "percentile"
            contrast_opts -- options for the contrast algorithm, see
                             the optional args of [contrast]_range()
                             for what to name the keys
            scale         -- how to scale the pixel values between the
                             min/max values from the contrast
                             algorithm when converting to a a raster
                             format, either "linear" or "arcsinh"
            scale_opts    -- options for the scaling algorithm, currently
                             only "nonlinearity" is supported for arcsinh,
                             which has a default value of 3
    """

    if contrast not in ("zscale", "percentile"):
        raise ValueError("invalid contrast algorithm '%s'" % contrast)
    if scale not in ("linear", "arcsinh"):
        raise ValueError("invalid scale value '%s'" % scale)

    # open the fits file and read the image data and size
    fitslib.fits_simple_verify(fitsfile)
    fits = pyfits.open(fitsfile)

    try:
        hdr = fits[0].header
        xsize = hdr["NAXIS1"]
        ysize = hdr["NAXIS2"]
        fits_data = fits[0].data
    finally:
        fits.close()
    
    # compute the proper scaling for the image
    if contrast == "zscale":
        contrast_value = contrast_opts.get("contrast", 0.25)
        num_points = contrast_opts.get("num_points", 600)
        num_per_row = contrast_opts.get("num_per_row", 120)
        zmin, zmax = zscale_range(fits_data, contrast=contrast_value,
                                  num_points=num_points,
                                  num_per_row=num_per_row)
    elif contrast == "percentile":
        min_percent = contrast_opts.get("min_percent", 3.0)
        max_percent = contrast_opts.get("max_percent", 99.0)
        num_points = contrast_opts.get("num_points", 5000)
        num_per_row = contrast_opts.get("num_per_row", 250)
        zmin, zmax = percentile_range(fits_data, min_percent=min_percent,
                                      max_percent=max_percent,
                                      num_points=num_points,
                                      num_per_row=num_per_row)

    # set all points less than zmin to zmin and points greater than
    # zmax to zmax
    fits_data = numpy.where(fits_data > zmin, fits_data, zmin)
    fits_data = numpy.where(fits_data < zmax, fits_data, zmax)

    if scale == "linear":
        scaled_data = (fits_data - zmin) * (255.0 / (zmax - zmin)) + 0.5
    elif scale == "arcsinh":
        # nonlinearity sets the range over which we sample values of the
        # asinh function; values near 0 are linear and values near infinity
        # are logarithmic
        nonlinearity = scale_opts.get("nonlinearity", 3.0)
        nonlinearity = max(nonlinearity, 0.001)
        max_asinh = cmath.asinh(nonlinearity).real
        scaled_data = (255.0 / max_asinh) * \
                      (numpy.arcsinh((fits_data - zmin) * \
                                     (nonlinearity / (zmax - zmin))))

    # convert to 8 bit unsigned int ("b" in numpy)
    #scaled_data = scaled_data.astype("b")
    
    # create the image
    print 'hi'
    image = Image.frombuffer("L", (xsize, ysize), scaled_data, "raw", "L", 0, 0)
    return image

def draw_circle(image, x, y, radius, color):

    """
    Draws a circle on image at position x, y with the given radius and
    color.

    Input:  image  -- the image object to draw the circle on
            x      -- the x position of the center of the circle
            y      -- the y position of the center of the circle
            radius -- the radius of the circle in pixels
            color  -- a tuple containing the color of the border of the
                      circle, ranging from 0 to 255 for each channel
    """

    # arc takes the upper left and lower right corners of a box bounding the
    # circle as arguments. Here (x1, y1) gives the coordinates of the upper left
    # corner and (x2, y2) gives the lower right corner of the bounding box.
    x1 = int(x - radius + 0.5)
    y1 = int(y - radius + 0.5)
    x2 = int(x + radius + 0.5)
    y2 = int(y + radius + 0.5)

    xsize, ysize = image.size

    # draw the circle
    draw = ImageDraw.Draw(image)
    draw.arc((x1, y1, x2, y2), 0, 360, fill=color)

def main(argv):
    import time

    if len(argv) != 2:
        print "Usage: %s <fits-file>" % os.path.basename(argv[0])
        print "Input file will be converted to JPEG"
        sys.exit(2)
   
    # FITS image to open and JPEG counterpart
    fitsfile = argv[1]
    name, ext = os.path.splitext(fitsfile)
    jpegfile = "%s.jpg" % name

    # open as PIL object
    start = time.time()
    image = FitsImage(fitsfile)#.convert("RGB")
    stop = time.time()
    print "Converting to PIL object took %f sec" % (stop - start)
    image.show()
    # save as a jpeg
    start = time.time()
    image.save(jpegfile)
    stop = time.time()
    print "Saving to '%s' took %f sec" % (jpegfile, stop - start)

if __name__ == "__main__":
    main(sys.argv)
