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
# 9/25/07   Added call to fits_simple_verify() to verify input file is FITS.
#
# 9/19/07  Modified to use the pyfits FITS I/O library and numpy. Modifications
#          made by Christopher Hanley, Space Telescope Science Institute.

"""
Takes a multi-extension FITS file and writes each extension to a separate
single extension FITS file
"""

import os
import sys
import fitslib
import pyfits

def main(argv):
  if len(argv) != 2:
    print "Usage: %s <FITS file>" % os.path.basename(argv[0])
    sys.exit(2)

  fitsfile = argv[1]
  if not os.path.exists(fitsfile):
    print "Error: file '%s' doesn't exist" % fitsfile
    sys.exit(2)

  name, ext = os.path.splitext(fitsfile)
  print "Extacting extensions from '%s'... " % fitsfile
  fitslib.fits_simple_verify(fitsfile)
  fits = pyfits.open(fitsfile)

  for i, hdu in enumerate(fits):
    new_fitsfile = "%s_ext_%d.fits" % (name, i)

    # Write out the new file
    pyfits.writeto(new_fitsfile, hdu.data, hdu.header)
    print "Wrote extension %d to %s" % (i, new_fitsfile)

  fits.close()

if __name__ == "__main__":
  main(sys.argv)
