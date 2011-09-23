#!/usr/bin/python

# Copyright (C) 2007 Google Inc.
# 
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

"""
Wrapper function for fitsregionator.py.  Use command line flags to specify
the parameters of the regionation work done by RegionateFITS.  If you plan
to make this data available over the web, use the "--url" flag to specify
the URL where the regionated kml in the output data will be stored.  This
means that others can access your data over the web by just downloading the
output.kml file (which is basically just a link to your data), rather than
downloading the entire data set in one go.
"""

import os
import sys
import getopt
import fits.fitsregionator

def usage(script_name):
  print """
    usage: %s [options] input.fit output.kml output_directory

    Basic Options:
    --orderby_field <string>        Choose a field in your FITs file to order
                                    your placemarks

    --header <integer>              Choose a FITS header to parse (default=1)

    --survey <string>               Specify a survey name to use when naming
                                    objects (default='MySurvey')

    --icon <stringURL>              URL to use for placemark icons.

    --url <stringURL>               URL where the regionated kml will be
                                    stored.


    Advanced Options:
    --min_lod <integer>             Minimum region size (in pixels) to use for
                                    loading deeper data (default=512)

    --max_lod <integer>             Maximum region size (in pixels) to use for
                                    loading deeper data (default=-1, if set to
                                    a positive value, this will make upper
                                    layer data disappear as you zoom in)

    --objects_per_region <integer>  Maximum number of objects in any given
                                    region.  When a region reaches this level
                                    any additional points are written to
                                    sub-regions.
    """ % script_name


def main():
  try:
    opts, args = getopt.getopt(sys.argv[1:], "o:h:s:i:r:u:l:n:",
                               ["orderby_field=",
                                "header=",
                                "survey=",
                                "icon=",
                                "url=",
                                "min_lod=",
                                "max_lod=",
                                "objects_per_node="])
  except getopt.GetoptError:
    # print help information and exit:
    print "Illegal option."

    print usage(os.path.basename(sys.argv[0]))

    sys.exit(2)

  orderbyField = ''
  hdu = 1
  surveyName = 'MySurvey'
  iconURL = 'http://mw1.google.com/mw-earth-vectordb/sky/sky1/pics/icon.png'
  minLod = 128
  maxLod = -1
  maxPerNode = 50
  url = ''

  for o, a in opts:
    if o == "--orderby_field":
      orderbyField = a
    if o == "--header":
      hdu = int(a)
    if o == "--survey":
      surveyName = a
    if o == "--icon":
      iconURL = a
    if o == "--url":
      url = a
    if o == "--min_lod":
      minLod = int(a)
    if o == "--max_lod":
      maxLod = int(a)
    if o == "--objects_per_node":
      maxPerNode = int(a)

  if len(args) != 3:
    print "Incorrect number of arguments."

    print usage(os.path.basename(sys.argv[0]))

    sys.exit(2)


  fitsfile = args[0]
  rootkml = args[1]
  outputDir = args[2]

  os.makedirs(outputDir)

  rtor = fits.fitsregionator.RegionateFITS(fitsfile,minLod,maxPerNode,rootkml,
                                           outputDir,surveyName,iconURL,
                                           orderbyField=orderbyField,
                                           hdu=hdu,maxLodPixels=maxLod,
                                           outputUrl=url)

  if not rtor:
    status = -1
  else:
    status = 0

  sys.exit(status)

if __name__ == "__main__":
    main()
