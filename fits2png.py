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

"""
Converts all input FITS images to PNG using the zscale autocontrast algorithm.
"""

import os
import sys
import optparse
import fitsimage

def main(argv):
  usage = "%prog [options] <list of FITS files>"
  parser = optparse.OptionParser(usage)

  # Options controlling how to perform the autocontrast and scaling.
  parser.add_option("--contrast", "-c", default="zscale", type="string",
                    help="autocontrast algorithm to apply, either 'zscale', " \
                         "'percentile', or 'manual', default %default")
  parser.add_option("--scale", "-s", default="linear", type="string",
                    help="how to scale pixel values between min and max, " \
                         "either 'linear' or 'arcsinh', default %default")

  # Options controlling both autocontrast algorithms.
  parser.add_option("--num_points", default=5000, type="int",
                    help="# of points to subsample image with")
  parser.add_option("--num_per_row", default=250, type="int",
                    help="# of points to sample from each row")

  # Options controlling zscale autocontrast.
  parser.add_option("--contrast_value", default=0.25, type="float",
                    help="contrast option for zscale algorithm")

  # Options controlling percentile autocontrast.
  parser.add_option("--min_percent", default=3.0, type="float",
                    help="black level for percentile")
  parser.add_option("--max_percent", default=99.5, type="float",
                    help="white level for percentile")

  # Options controlling manual contrast.
  parser.add_option("--min", default=0.0, type="float",
                    help="black level for manual")
  parser.add_option("--max", default=32768.0, type="float",
                    help="white level for manual")

  # Options controlling arcsinh scaling.
  parser.add_option("--nonlinearity", default=3.0, type="float",
                    help="nonlinearity option for arcsinh scaling, " \
                         "default %default")

  flags, files = parser.parse_args(argv[1:])
  
  if len(files) < 1:
    print "usage: %s <list of FITS file>" % os.path.basename(argv[0])
    sys.exit(2)

  # Set up the options for the types of contrast algorithms.
  contrast_opts = {}
  contrast_opts["num_points"] = flags.num_points
  contrast_opts["num_per_row"] = flags.num_per_row

  if flags.contrast == "zscale":
    contrast_opts["contrast"] = flags.contrast_value
  elif flags.contrast == "percentile":
    contrast_opts["min_percent"] = flags.min_percent
    contrast_opts["max_percent"] = flags.max_percent
  elif flags.contrast == "manual":
    contrast_opts["min"] = flags.min
    contrast_opts["max"] = flags.max
  else:
    print "Error: invalid contrast algorithm '%s'" % flags.contrast
    sys.exit(2)

  if flags.scale not in ("linear", "arcsinh"):
    print "Error: invalid scale '%s'" % flags.scale
    sys.exit(2)

  scale_opts = {}
  scale_opts["nonlinearity"] = flags.nonlinearity
    
  for fitsfile in files:
    if not os.path.exists(fitsfile):
      print "Error: file '%s' doesn't exist" % fitsfile
      sys.exit(2)

    name, ext = os.path.splitext(fitsfile)
    pngfile = "%s.png" % name

    sys.stderr.write("Converting '%s'... " % fitsfile)
    image = fitsimage.FitsImage(fitsfile, contrast=flags.contrast,
                                contrast_opts=contrast_opts, scale=flags.scale,
                                scale_opts=scale_opts)
    image.save(pngfile)
    sys.stderr.write("done\n")

if __name__ == "__main__":
  main(sys.argv)
