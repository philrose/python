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
# 9/19/07  Replaced the fitslib module with the pyfits library.
#          Added the parse_card function to this script.  Modification made by
#          Christopher Hanley, Space Telescope Science Institute.

"""
Copies WCS information from an ascii file with FITS-like keywords into a
proper FITS header.
"""

import os
import sys
import fitslib
import pyfits

WCS_KEYWORDS = ["EPOCH", "EQUINOX", "CTYPE1", "CTYPE2", "CRPIX1", "CRPIX2",
                "CRVAL1", "CRVAL2", "CDELT1", "CDELT2", "CD1_1", "CD1_2",
                "CD2_1", "CD2_2", "WCSAXES", "CROTA", "CROTA2", "PC1_1",
                "PC1_2", "PC2_1", "PC2_2", "PC001001", "PC001002", "PC002001",
                "PC002002", "LATPOLE", "LONPOLE"]

def main(argv):
  if len(argv) < 2:
    print "Usage: %s <list of WCS ascii files>" % os.path.basename(argv[0])
    sys.exit(2)

  for filename in argv[1:]:
    name, ext = os.path.splitext(filename)
    fitsname = "%s.fits" % name
    print "Copying WCS from %s to FITS file %s..." % (filename, fitsname)

    header = {}
    f = open(filename)
    for line in f:
      name, value, comment = fitslib.parse_card(line)
      if name is not None and value is not None:
        header[name] = value
    f.close()

    # Create a new primary header with associated data.
    prihdu = pyfits.PrimaryHDU()

    # Create a new HDUList object.
    new_fits = pyfits.HDUList([prihdu])
    
    fits_header = new_fits[0].header
    for key in WCS_KEYWORDS:
      if key in header:
        fits_header.update(key, header[key])

    # Clean up WCS so only a CD matrix is included.
    if fits_header.has_key("CDELT1") or fits_header.has_key("CDELT2"):
      has_cd = False
      for key in ("CD1_1", "CD1_2", "CD2_1", "CD2_2"):
        if fits_header.has_key(key):
          has_cd = True

      if has_cd:
        for key in ("CDELT1", "CDELT2"):
          if fits_header.has_key(key):
            del fits_header[key]

        for key in ("CD1_1", "CD1_2", "CD2_1", "CD2_2"):
          if not fits_header.has_key(key):
            fits_header.update(key, 0.0)

    # Add a J2000 EPOCH if needed.
    if not fits_header.has_key("EPOCH") and not fits_header.has_key("EQUINOX"):
      print "Added a J2000 equinox"
      fits_header.update("EQUINOX", 2000.0)

    new_fits.writeto(fitsname)
    new_fits.close()

if __name__ == "__main__":
  main(sys.argv)
