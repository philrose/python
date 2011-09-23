#!/usr/bin/env python

# Copyright 2007 Google Inc.
# Author: Jeremy Brewer
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#      http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#
# Changelog:
#
# 9/25/07  Added call to fits_simple_verify() to verify input file is FITS.
#
# 9/19/07  Replaced the fitslib.py module with the pyfits library.
#          Modification made by Christopher Hanley, Space Telescope Science
#          Institute.

"""
Converts an image with WCS information to KML + warped image

Takes an image file (jpeg, tiff, etc.) and a FITS file containing WCS
information and warped the input image to be correctly display in Sky in
Google Earth.

Example usage:
wcs2kml.py --imagefile=<raster image> --fitsfile=<FITS file with WCS>

If the given imagefile is a FITS file, then only supply the --imagefile
flag (ignoring --fitsfile) and it will automatically be converted to raster
format using the fitsimage module.

Requirements:
- Python 2.3 or newer
- Python Imaging Library (Image)
- pyfits
- numpy
- Included FITS libraries: fitslib, fitsimage and wcslib
"""

import os
import sys
import math
import getopt
import time
import fitslib
import pyfits
import fitsimage
import wcslib
import Image

def round_float(value):
  """Rounds a float to the nearest integer value."""
  return int(value + 0.5)

class WrapAround(object):
  
  """
  Static class for determining if images cross the 0-360 boundary

  Nearly all of the source code that deals with spherical coordinates needs
  to deal properly with the 0-360 discontinuity.  This simple class exists
  so that it is possible to consistently deal with the discontinuity.

  Example usage:

  # List of ra values.
  ra = []

  # ...
  # Code for reading ra values omitted
  # ...

  ra_min = min(ra)
  ra_max = max(ra)

  # Does the image wrap around?
  is_wrapped = WrapAround.imageWrapsAround(ra_min, ra_max)

  # Modify the ra values so that they increase monotonically from
  # ra_min to ra_max.
  for i in xrange(len(ra)):
    ra[i] = WrapAround.makeRaMonotonic(ra[i])
  
  # Modify ra values so that they wrap from 360 back to 0.
  for i in xrange(len(ra)):
    ra[i] = WrapAround.restoreWrapAround(ra[i])
  """

  # Constant that defines when an image wraps around and which ra values
  # should be modified.
  MAX_DELTA_RA = 180.0
  THREE_SIXTY = 360.0
  
  def __init__(self):
    raise NotImplementedError("This class isn't intended to be instantiated")

  def imageWrapsAround(ra_min, ra_max):

    """
    Returns true when an image wraps around the 0-360 discontinuity.  This
    method should be called with the min and max ra values from a given
    image.
    """

    return abs(ra_min - ra_max) > WrapAround.MAX_DELTA_RA

  imageWrapsAround = staticmethod(imageWrapsAround)

  def makeRaMonotonic(ra):

    """
    Raises the input ra above 360 if the point wraps around.  Points with ra
    greater than MAX_DELTA_RA are increased by 360 because by definition
    no image can span from 0 to more than MAX_DELTA_RA.

    Note that this method should not be called for every projected point as
    it will simply move the discontinuity to MAX_DELTA_RA.  One should first
    determine if the image wraps around by finding the max and min ra
    values in an image and calling imageWrapsAround(), then apply this
    function to each point in the image.

    The function restoreWrapAround() is the inverse of this function.
    """

    if ra < WrapAround.MAX_DELTA_RA:
      return ra + WrapAround.THREE_SIXTY
    else:
      return ra

  makeRaMonotonic = staticmethod(makeRaMonotonic)

  def restoreWrapAround(ra):

    """
    Adjusts the input ra to lie within the proper 0-360 bounds.  This
    function is safe to use for all ra input and will work for for any
    value of ra, no matter how far away from the 0-360 limits it is.
    """

    ra_wrap = ra
    while ra_wrap > WrapAround.THREE_SIXTY:
      ra_wrap -= WrapAround.THREE_SIXTY
    while ra_wrap < 0.0:
      ra_wrap += WrapAround.THREE_SIXTY

    return ra_wrap

  restoreWrapAround = staticmethod(restoreWrapAround)

class Point(object):
  
  """
  Simple class for holding an pair of ra, dec and x, y coordinates

  This class holds a pair of spherical coordinates on the sky and their
  corresponding pixel coordinates in a raster image of this portion of the
  sky.  The right ascension (ra), declination (dec), x coordinate (x) and
  y coordinate (y) are the only fields.

  The functions for finding the bounding box of an image locates the four
  corners in projected space.  Because 16 coordinates are a lot to keep up
  with, the functions in this header return the coordinates as Points.
  """
  
  def __init__(self, ra=0.0, dec=0.0, x=0.0, y=0.0):
    self.setValues(ra, dec, x, y)

  def __str__(self):
    return "ra, dec: %.14f, %.14f => x, y: %.14f, %.14f" % (self.ra, self.dec,
                                                            self.x, self.y)

  def setValues(self, ra, dec, x, y):
    assert isinstance(ra, float)
    assert isinstance(dec, float)
    assert isinstance(x, float)
    assert isinstance(y, float)
    self.ra = ra
    self.dec = dec
    self.x = x
    self.y = y

  def distanceXY(self, point):
    
    """
    Computes the Euclidean distance between x and y of two points.
    """
    
    dx = point.x - self.x
    dy = point.y - self.y
    return math.sqrt(dx ** 2 + dy ** 2)

  def distanceRaDec(self, point):

    """
    Computes the Euclidean distance between ra and dec of two points.
    
    This distance is incorrect in real world coordinates because it does not
    account for the curvature of the sky.  However, it is correct in the lat-lon
    projection plane.
    """

    dr = point.ra - self.ra
    dd = point.dec - self.dec
    return math.sqrt(dr ** 2 + dd ** 2)

  def distanceRaDecExact(self, point):
    # Convert to radians from degrees.
    r1 = math.radians(self.ra)
    d1 = math.radians(self.dec)
    r2 = math.radians(point.ra)
    d2 = math.radians(point.dec)

    # Convert to x, y, z space on a sphere.
    x1 = math.cos(d1) * math.cos(r1)
    y1 = math.cos(d1) * math.sin(r1)
    z1 = math.sin(d1)

    x2 = math.cos(d2) * math.cos(r2)
    y2 = math.cos(d2) * math.sin(r2)
    z2 = math.sin(d2)

    dx = x1 - x2
    dy = y1 - y2
    dz = z1 - z2

    return math.sqrt(dx ** 2 + dy ** 2 + dz ** 2)

class BoundingBox(object):
  
  """
  Class for holding the bounding box of images in both x, y and ra, dec space

  Given an input image and a non-linear sperical projection, one can find
  the range in spherical coordinates that the image spans.  In general for
  a non-linear distortion these extrema will not occur at the corners.
  Similarly, there will not be symmetry between pairs of coordinates; i.e.,
  the declination for the maximum right ascension will not be the same as
  that for the minimum right ascension.  It is therefore necessary to keep
  up with both coordinates for each of the 4 corners.  This class is designed
  to make keeping up with all of these pairs of coordinates easier.

  Example usage:

  # This computes the bounding box and determines whether the image wraps
  # around the 0-360 discontinuity.
  box = BoundingBox(wcs, width, height)

  # Find out if the box wraps around 0-360.
  is_wrapped = box.isWrapped()

  # Access the corners of the box in both pixel and spherical coordinates.
  ra_min = box.ra_min
  ra = ra_min.ra
  dec = ra_min.dec
  x = ra_min.x 
  y = ra_min.y
  """
  
  # Constants for error checking.
  LARGE_BAD_VALUE = 999.0
  SMALL_BAD_VALUE = -999.0
  
  def __init__(self, wcs, width, height):
    self._is_wrapped = False
    self._ra_min = Point()
    self._ra_max = Point()
    self._dec_min = Point()
    self._dec_max = Point()
    self.findBoundingBox(wcs, width, height)

  def getRaMin(self):
    return self._ra_min

  ra_min = property(fget=getRaMin, doc="Access ra_min")

  def getRaMax(self):
    return self._ra_max

  ra_max = property(fget=getRaMax, doc="Access ra_max")

  def getDecMin(self):
    return self._dec_min

  dec_min = property(fget=getDecMin, doc="Access dec_min")

  def getDecMax(self):
    return self._dec_max

  dec_max = property(fget=getDecMax, doc="Access dec_max")

  def isWrapped(self):
    """Returns whether the bounding box wraps around the 0-360 limit."""
    return self._is_wrapped

  def _updateExtrema(self, wcs, x, y):
  
    """
    Checks the given point to determine if it is beyond the current extrema.
    The values of _ra_min, _ra_max, _dec_min, and _dec_max are updated if
    needed.
    """

    ra, dec = wcs.toRaDec(x, y)

    # Update ra to be monotonic across the image if the image wraps around
    # the 0-360 discontinuity.
    if self._is_wrapped:
      ra = WrapAround.makeRaMonotonic(ra)

    if ra > self._ra_max.ra:
      self._ra_max.setValues(ra, dec, x, y)
    elif ra < self._ra_min.ra:
      self._ra_min.setValues(ra, dec, x, y)

    # Check dec separately because the min or max ra, dec could occur at
    # the same point.
    if dec > self._dec_max.dec:
      self._dec_max.setValues(ra, dec, x, y)
    elif dec < self._dec_min.dec:
      self._dec_min.setValues(ra, dec, x, y)

  def _findBoundingBoxForKnownWrapped(self, wcs, width, height):

    """
    Internal method for determining the bounding box.  This method must be
    called twice to properly detemine the bounding box if the 0-360 boundary
    is crossed.
    """

    # Intentionally bad values.
    self._ra_min = Point(ra=BoundingBox.LARGE_BAD_VALUE,
                         dec=BoundingBox.LARGE_BAD_VALUE)
    self._ra_max = Point(ra=BoundingBox.SMALL_BAD_VALUE,
                         dec=BoundingBox.SMALL_BAD_VALUE)
    self._dec_min = Point(ra=BoundingBox.LARGE_BAD_VALUE,
                          dec=BoundingBox.LARGE_BAD_VALUE)
    self._dec_max = Point(ra=BoundingBox.SMALL_BAD_VALUE,
                          dec=BoundingBox.SMALL_BAD_VALUE)

    # The loops here sum integers instead of doubles for increased numerical
    # accuracy.  The added accuracy shouldn't be needed, but these loops aren't
    # performance critical, so we just keep the extra precision just in case.
    # NB: There are 2 duplicate corner calculations, wasting 4 calculations.

    # 1st edge
    y = 1.0
    for i in xrange(1, width + 1):
      x = float(i)
      self._updateExtrema(wcs, x, y)

    # 2nd edge
    y = float(height)
    for i in xrange(1, width + 1):
      x = float(i)
      self._updateExtrema(wcs, x, y)

    # 3rd edge
    x = 1.0
    for i in xrange(1, height + 1):
      y = float(i)
      self._updateExtrema(wcs, x, y)

    # 4th edge
    x = float(width)
    for i in xrange(1, height + 1):
      y = float(i)
      self._updateExtrema(wcs, x, y)

    # Sanity checks.
    assert self._ra_min.ra < BoundingBox.LARGE_BAD_VALUE
    assert self._ra_min.dec < BoundingBox.LARGE_BAD_VALUE
    
    assert self._ra_max.ra > BoundingBox.SMALL_BAD_VALUE
    assert self._ra_max.dec > BoundingBox.SMALL_BAD_VALUE

    assert self._dec_min.ra < BoundingBox.LARGE_BAD_VALUE
    assert self._dec_min.dec < BoundingBox.LARGE_BAD_VALUE
    
    assert self._dec_max.ra > BoundingBox.SMALL_BAD_VALUE
    assert self._dec_max.dec > BoundingBox.SMALL_BAD_VALUE

  def findBoundingBox(self, wcs, width, height):

    """
    Searches every edge pixel of an image to determine the spherical
    coordinate bounding box for the image.  This method also determines
    whether the image wraps around the 0-360 discontinuity.
    """

    assert isinstance(wcs, wcslib.BaseWcsProjection)
    assert isinstance(width, int)
    assert isinstance(height, int)

    self._is_wrapped = False
    self._findBoundingBoxForKnownWrapped(wcs, width, height)

    # If the image wraps around the 0-360 discontinuity, the max and min values
    # will be incorrect because ra is not monotonic across the image.  To fix
    # this we flag the image as wrapped and re-compute the max and min values.
    # findBoundingBox() adjusts ra to be monotonic across the image if the
    # is_wrapped_ flag is true.
    if WrapAround.imageWrapsAround(self._ra_min.ra, self._ra_max.ra):
      self._is_wrapped = True
      self._findBoundingBoxForKnownWrapped(wcs, width, height)

  def toKml(self, imagefile):
  
    """
    Generates a KML representation of the bounding box and returns it as
    a string.  The KML references the image whose filename is imagefile.
    
    The KML representation includes a <GroundOverlay> element which describes
    the bounding box of imagefile.  This method allows for the viewing of
    images in Google Earth, but imagefile must be warped to lat-lon projection.
    """
    
    # We need to ensure that ra increases monotonically from so that the
    # center calculation is correct.
    rmin = self._ra_min.ra
    rmax = self._ra_max.ra
    if rmax < rmin:
      rmax = WrapAround.makeRaMonotonic(rmax)

    ra_center = 0.5 * (rmin + rmax)
    ra_center = WrapAround.restoreWrapAround(ra_center)
    dec_center = 0.5 * (self._dec_min.dec + self._dec_max.dec)

    # Compute the bounding boxes.  Note that the ra values here should range
    # from 0 to 360 (see note below).
    east = WrapAround.restoreWrapAround(self._ra_max.ra)
    west = WrapAround.restoreWrapAround(self._ra_min.ra)
    north = self._dec_max.dec
    south = self._dec_min.dec

    assert 0.0 < east < 360.0
    assert 0.0 < west < 360.0
    assert 0.0 < ra_center < 360.0

    # Now we can convert to the KML -180 to 180 coordinate system.
    # NB: The ra values here must be between 0 and 360 for this conversion
    # to be correct (otherwise it's -360).
    ra_center -= 180.0
    east -= 180.0
    west -= 180.0

    # Google Earth assumes that east > west, otherwise it seems to switch
    # the coordinates (maybe it thinks you swapped them by mistake?).  For
    # some reason though, increasing east doesn't work, but decreasing west
    # does.
    if east < west:
      west -= 360.0

    kml = []
    kml.append('<?xml version="1.0" encoding="UTF-8"?>')
    kml.append('<kml xmlns="http://earth.google.com/kml/2.2" ' \
               'hint="target=sky">')
    kml.append("<Document>")
    kml.append("<GroundOverlay>")
    kml.append("  <name>Your uploaded image!</name>")
    kml.append("  <Icon>")
    kml.append("    <href>%s</href>" % imagefile)
    kml.append("  </Icon>")
    kml.append("  <LatLonBox>")
    kml.append("    <north>%.14f</north>" % north)
    kml.append("    <south>%.14f</south>" % south)
    kml.append("    <east>%.14f</east>" % east)
    kml.append("    <west>%.14f</west>" % west)
    kml.append("    <rotation>0.0</rotation>")
    kml.append("  </LatLonBox>")
    kml.append("  <LookAt>")
    kml.append("    <longitude>%.14f</longitude>" % ra_center)
    kml.append("    <latitude>%.14f</latitude>" % dec_center)
    kml.append("    <altitude>0.0</altitude>")
    kml.append("    <range>100000</range>")
    kml.append("    <heading>0.0</heading>")
    kml.append("    <altitudeMode>clampToGround</altitudeMode>")
    kml.append("  </LookAt>")
    kml.append("</GroundOverlay>")
    kml.append("</Document>")
    kml.append("</kml>")
    
    return "\n".join(kml)

class Matrix(object):
  
  """Class for accessing data in a flattened array."""
  
  def __init__(self, width, height, data=None, bg_color=(0, 0, 0, 0)):
    assert isinstance(width, int)
    assert isinstance(height, int)

    self.width = width
    self.height = height

    if not data:
      self.data = [bg_color] * (width * height)
    else:
      self.data = data

  def __getitem__(self, indeces):
    i = indeces[0]
    j = indeces[1]
    return self.data[i + j * self.width]

  def __setitem__(self, indeces, value):
    i = indeces[0]
    j = indeces[1]
    self.data[i + j * self.width] = value

class ImageOrigin(object):

  """
  Enumeration class for holding constants that specify where the origin in
  pixel space of an image is.

  Most image raster formats use an origin at the upper left corner of an
  image, but the FITS format assumes an origin in the lower left corner.
  It is therefore necessary to keep track of where the pixel origin is 
  when projection the image, otherwise it will be flipped about the input
  image y axis.

  Also note that because the conversion from FITS to raster format may or
  may not have corrected for this discrepancy already, so not all FITS
  images may need an origin of LOWER_LEFT.
  """

  UPPER_LEFT = 0
  LOWER_LEFT = 1

  def __init__(self):
    raise NotImplementedError("This class isn't intended to be instantiated")

  def validate(origin):
    assert origin in (ImageOrigin.UPPER_LEFT, ImageOrigin.LOWER_LEFT)

  validate = staticmethod(validate)

class SkyProjection(object):
  
  """
  Class that performs the image projection to lat-lon
  """
   
  def __init__(self, image, wcs):
    assert isinstance(image, Image.Image)
    assert isinstance(wcs, wcslib.BaseWcsProjection)

    # Convert input image to RGBA.
    self._image = image.convert("RGBA")
    self._wcs = wcs
    width, height = image.size

    # Bounding box for this image in lat-lon projection.
    self._bounding_box = BoundingBox(wcs, width, height)

    # Background color in RGBA for areas outside the input image.  The default
    # value is transparent.
    self._bg_color = (0, 0, 0, 0)
    
    # Where the origin in pixel space of the input image is.
    self._input_image_origin = ImageOrigin.UPPER_LEFT

    # Dimensions of the output projected image.  These must be set or
    # autmatically determined before the image can be projected.
    self._projected_width = 0
    self._projected_height = 0
    self._determineProjectedSize()

  def _determineProjectedSize(self):
    
    """
    Computes the size of the output dimensions by comparing pixel distances
    along east and west and updates _projected_width and _projected_height.
    """
    
    # Determine the side lengths along north and east.
    ra_min = self._bounding_box.ra_min
    ra_max = self._bounding_box.ra_max
    dec_min = self._bounding_box.dec_min

    east_side_len = round_float(ra_max.distanceXY(dec_min))
    north_side_len = round_float(ra_min.distanceXY(dec_min))

    # Compute the angle between east (x axis in new image) and the x axis in
    # the old image.
    cos_theta = (ra_max.ra - dec_min.ra) / ra_max.distanceRaDec(dec_min)
    theta = math.acos(cos_theta)

    # Known Issue:
    # When theta --> 0, my diagram of the geometry of the projected image
    # breaks down because the extrema move to the corners.  When this happens
    # there is ambiguity about which corner corresponds to the point I call
    # ra_min in my diagram: I assume it should be the LL corner, but it could
    # mathematically be the UL corner.  This is true for all of the corners --
    # they could be at the top or bottom of the projected image.
    #
    # This causes the image dimensions to be horribly wrong, but the
    # bounding box is still correct, and Earth seems to stretch the image OK.
    # Unfortunately I can't come up with a nice real fix for this, so I'm
    # just going to leave it until it breaks.

    # Determine the size of the output image.  This size is chosen so that the
    # original image remains the same dimensions within the projected image.
    self._projected_width = \
        round_float(math.cos(theta) * east_side_len + \
                    math.sin(theta) * north_side_len)
    self._projected_height = \
        round_float(math.sin(theta) * east_side_len + \
                    math.cos(theta) * north_side_len)

    assert self._projected_width > 0
    assert self._projected_height > 0

  def getBoundingBox(self):
    """Returns the BoundingBox object computed for this image."""
    return self._bounding_box

  boundingBox = property(fget=getBoundingBox, doc=getBoundingBox.__doc__)

  def getInputImageOrigin(self):
    return self._input_image_origin

  def setInputImageOrigin(self, origin):
    ImageOrigin.validate(origin)
    self._input_image_origin = origin

  inputImageOrigin = property(fget=getInputImageOrigin,
                              fset=setInputImageOrigin,
                              doc="Access the input image origin value")

  def setMaxSideLength(self, max_side_length):

    """
    Adjusts the output dimensions so that they do not exceed max_side_length
    pixels, maintaining the same aspect ratio.
    """

    assert isinstance(max_side_length, int)
    assert max_side_length > 0

    # Enforce a maximum side length that keeps the same aspect ratio.
    if self._projected_width > self._projected_height and \
       self._projected_width > max_side_length:
      old_width = float(self._projected_width)
      new_width = float(max_side_length)
      height = float(self._projected_height)
      self._projected_width = max_side_length
      self._projected_height = round_float(height * (new_width / old_width))
    elif self._projected_height > self._projected_width and \
         self._projected_height > max_side_length:
      old_height = float(self._projected_height)
      new_height = float(max_side_length)
      width = float(self._projected_width)
      self._projected_height = max_side_length
      self._projected_width = round_float(width * (new_height / old_height))

    assert self._projected_width > 0
    assert self._projected_height > 0

  maxSideLength = property(fset=setMaxSideLength, doc=setMaxSideLength.__doc__)

  def getProjectedWidth(self):
    return self._projected_width
    
  def setProjectedWidth(self, width):
    assert isinstance(width, int)
    self.projected_width_ = width

  projectedWidth = property(fget=getProjectedWidth, fset=setProjectedWidth,
                            doc="Access the width of the projected image")

  def getProjectedHeight(self):
    return self._projected_height
     
  def setProjectedHeight(self, height):
    assert isinstance(height, int)
    self.projected_height_ = height

  projectedHeight = property(fget=getProjectedHeight, fset=setProjectedHeight,
                             doc="Access the height of the projected image")

  def getBackgroundColor(self):
    return self._bg_color

  def setBackgroundColor(self, color):
    
    """
    Sets the background color that will be output for pixels outside of the
    input image when warping the image.  The input color must have 4
    channels.
    """

    assert len(color) == 4
    self._bg_color = tuple([x for x in color])

  backgroundColor = property(fget=getBackgroundColor, fset=setBackgroundColor,
                             doc="Access the background color of pixels " \
                                 "outside the original image")

  def warpImage(self):

    """
    Warps the underlying image to lat-lon projection.  The alpha channel of
    the input image is preserved if present.
    
    The warped image is returns as a Python Imaging Library object.
    """

    # Read the data from the original image in its original format.
    image_data = Matrix(self._image.size[0], self._image.size[1],
                        data=self._image.getdata())

    # Create the image for the projection.
    assert self._projected_width > 0
    assert self._projected_height > 0
    projected_image = Image.new("RGBA", (self._projected_width,
                                         self._projected_height),
                                self._bg_color)
    projected_image_data = Matrix(self._projected_width, 
                                  self._projected_height,
                                  bg_color=self._bg_color)

    ra_min = self._bounding_box.ra_min
    ra_max = self._bounding_box.ra_max
    dec_min = self._bounding_box.dec_min
    dec_max = self._bounding_box.dec_max

    # Scale factors for converting ra, dec to x, y in projected image.
    xscale = (ra_max.ra - ra_min.ra) / float(self._projected_width - 1)
    yscale = (dec_max.dec - dec_min.dec) / float(self._projected_height - 1)

    # For each pixel in the new image, compute ra, dec then x, y in the original
    # image and copy the pixel values.
    # NB: The loop procedes from (ra_max, dec_max) to (ra_min, dec_min), i.e.
    # from the upper left corner to the lower right corner of the projected
    # image.  Hence, i, j properly indexes the projected image in lat-lon
    # space.
    for i in xrange(self._projected_width):
      ra = ra_max.ra - i * xscale
      for j in xrange(self._projected_height):
        dec = dec_max.dec - j * yscale

        # Coordinates in original FITS image start at (1, 1) in the lower left
        # corner.  Here we convert them so that point (0, 0) is in the upper
        # left corner if needed.
        # NB: I tried doing this this by modifiying the WCS outside of the
        # inner loop, but the performance gains were minimal (5%).  The code
        # was much less clear though.
        x, y = self._wcs.toPixel(ra, dec, rnd=False)

        x -= 1.0
        if self._input_image_origin == ImageOrigin.LOWER_LEFT:
          y = self._image.size[1] - y
        else:
          y -= 1.0

        # Does this pixel lie inside the original image?
        inside = True
        if not 0 < x < self._image.size[0] or \
           not 0 < y < self._image.size[1]:
           inside = False

        if inside:
          # Copy the input pixel and add an alpha channel that is opaque.
          # Straight point sampling looks good because Earth applies its
          # own filtering.  We can ignore points outside of the input image
          # because the background color was set upon image creation.
          m = round_float(x)
          n = round_float(y)
          
          # Deal with points at the very outer edge.
          m = min(m, self._image.size[0] - 1)
          n = min(n, self._image.size[1] - 1)

          # This copy preserves the input alpha channel.
          projected_image_data[i, j] = image_data[m, n]
  
    # Update image data.
    projected_image.putdata(projected_image_data.data)
    return projected_image

def write_kmz(kmzfile, imgfile, kmlfile, other_files=None):
  
  """Compresses the given files into a KMZ archive.

  Args:
    kmzfile:     filename of output KMZ file
    imgilfe:     filename of input warped image file
    kmlfile:     filename of KML file (should be 'doc.kml')
    other_files: additional files to include in the KMZ archive
  """  
  
  if other_files is None:
    other_files = []

  if os.path.exists(kmzfile):
    os.remove(kmzfile)

  cmd = "zip %s %s %s %s" % (kmzfile, kmlfile, imgfile, " ".join(other_files))
  for out in os.popen(cmd):
    pass

def print_usage(argv):
  print
  print "Usage: %s [options] --imagefile=<image file> " \
        "--fitsfile=<fits file>" % argv[0]
  print
  print "Options:"
  print
  print "--input_image_origin_is_lower_left .. set origin of input image, "
  print "                                      default False"
  print "-h, --help .......................... displays this message"
  print "-o=, --outfile= ..................... name of output file, default"
  print "                                      <image prefix>.kmz"
  print "--max_side_length= .................. max side length in pixels of "
  print "                                      output image"
  print
  print "The image file can be a jpeg, gif, png, etc., and the FITS file"
  print "is only used when reading the WCS information."
  print

def main(argv):
  # Remove leading directory for error messages.
  argv[0] = os.path.basename(argv[0])
  
  # Default options.
  fitsfile = ""
  imagefile = ""
  input_image_origin_is_lower_left = False
  max_side_length = 3000
  outfile = ""

  # Parse command-line options.
  try:
    opts, files = getopt.getopt(argv[1:], "ho:", \
        ["fitsfile=", "help", "imagefile=", "input_image_origin_is_lower_left",
         "max_side_length=", "outfile="])
  except getopt.GetoptError, (msg, bad_opt):
    print "Error: %s" % msg
    print "Type '%s --help' for more information" % argv[0]
    sys.exit(2)

  # Process input flags.
  for opt, arg in opts:
    if opt == "--fitsfile":
      fitsfile = arg
    if opt == "--imagefile":
      imagefile = arg
    if opt == "--input_image_origin_is_lower_left":
      input_image_origin_is_lower_left = True
    if opt in ("-h", "--help"):
      print_usage(argv)
      sys.exit(0)
    if opt in ("-o", "--outfile"):
      outfile = arg
    if opt == "--max_side_length":
      try:
        max_side_length = int(arg)
      except ValueError:
        print "Error: invalid value '%s' for %s" % (arg, opt)
        sys.exit(2)

  # Check for an input image.
  if not imagefile:
    print "Error: no input imagefile was given"
    print "Type '%s --help' for more information" % argv[0]
    sys.exit(2)
  elif not os.path.exists(imagefile):
    print "Error: file '%s' doesn't exist" % imagefile
    sys.exit(2)

  # Determine if the input image is a FITS file.
  is_fits_file = True
  try:
    fitslib.fits_simple_verify(imagefile)
  except (IOError, ValueError, EOFError):
    is_fits_file = False
    if not fitsfile:
      print "Error: no FITS file input for --fitsfile"
      sys.exit(2)

  # Convert input FITS images to raster format if the --fitsfile option has
  # not been given.
  if is_fits_file and not fitsfile:
    sys.stdout.write("Converting input FITS image to raster format... ")
    sys.stdout.flush()
    image = fitsimage.FitsImage(imagefile)
    sys.stdout.write("done\n")
    fitsfile = imagefile    
    # This should always be True for fitsimage conversions.
    input_image_origin_is_lower_left = True
  else:
    image = Image.open(imagefile)

  # Check for input FITS file with WCS.
  if not os.path.exists(fitsfile):
    print "Error: file '%s' doesn't exist" % fitsfile
    sys.exit(2)

  # Read the WCS.
  sys.stdout.write("Reading WCS from %s... " % fitsfile)
  sys.stdout.flush()
  fitslib.fits_simple_verify(fitsfile)
  fits = pyfits.open(fitsfile)
  try:
    header = fits[0].header
    wcs = wcslib.WcsProjection(header)
  finally:
    fits.close()
  sys.stdout.write("done\n")
  
  # Set the various warping options and warp the input image.
  projection = SkyProjection(image, wcs)
  projection.backgroundColor = (0, 0, 0, 0)
  projection.maxSideLength = max_side_length
  if input_image_origin_is_lower_left:
    projection.inputImageOrigin = ImageOrigin.LOWER_LEFT
  else:
    projection.inputImageOrigin = ImageOrigin.UPPER_LEFT

  sys.stdout.write("Warping image... ")
  sys.stdout.flush()
  start = time.time()
  projected_image = projection.warpImage()
  stop = time.time()
  sys.stdout.write("done in %f sec\n" % (stop - start))

  # Write warped image to file.
  name, ext = os.path.splitext(imagefile)
  projected_imagefile = "%s_warped.png" % name

  sys.stdout.write("Writing warped image to %s... " % projected_imagefile)
  sys.stdout.flush()
  projected_image.save(projected_imagefile)
  sys.stdout.write("done\n")

  # Write bounding box to KML file.
  sys.stdout.write("Writing KML to doc.kml... ")
  sys.stdout.flush()
  bounding_box = projection.boundingBox
  f = open("doc.kml", "w")
  try:
    f.write(bounding_box.toKml(projected_imagefile))
  finally:
    f.close()
  sys.stdout.write("done\n")

  # Write KML and warped image to KMZ.
  if not outfile:
    outfile = "%s.kmz" % name
 
  sys.stdout.write("Writing KMZ to '%s'... " % outfile)
  sys.stdout.flush()
  write_kmz(outfile, projected_imagefile, "doc.kml")
  sys.stdout.write("done\n")
  
  # Clean up temporary files.
  sys.stdout.write("Cleaning up doc.kml and %s... " % projected_imagefile)
  sys.stdout.flush()
  os.remove("doc.kml")
  os.remove(projected_imagefile)
  sys.stdout.write("done\n")
  
  print "All done!"
  
if __name__ == "__main__":
  main(sys.argv)
