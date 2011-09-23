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
Convert a FITS data table into a collection of placemarks and write out to a
KML file.  The code will automatically look for two fields representing the
right ascension and declination for each object and exit if it is unable to
find any reasonable candidates.  Additionally, the user can specify a partcular
field to use when sorting the output (bearing in mind that the order will
affect objects' order of appearance if the resulting KML file is regionated).
If no field is specified, then the order will match that in the input FITS
file.
"""

import sys
import os
import getopt
import pyfits
import fits.fitsparse

def main():
  try:
    opts, args = getopt.getopt(sys.argv[1:], "o:h:s:i:", ["orderby-field=",
                                                          "header=",
                                                          "survey=",
                                                          "icon="])
  except getopt.GetoptError:
    # print help information and exit:
    print """
    usage: %s [options] input.fit output.kml

    Options:
    --orderby-field <string>  Choose a field in your FITs file to order
                              your placemarks

    --header <number>         Choose a FITS header to parse (default=1)

    --survey <string>         Specify a survey name to use when naming
                              objects (default='MySurvey')

    --icon <stringURL>        URL to use for placemark icons.
    """ % os.path.basename(sys.argv[0])
    sys.exit(2)

  orderbyField = ''
  header = 1
  surveyName = 'MySurvey'
  iconURL = 'http://mw1.google.com/mw-earth-vectordb/sky/sky1/pics/icon.png'

  for o, a in opts:
    if o == "--orderby-field":
      orderbyField = a
    if o == "--header":
      header = int(a)
    if o == "--survey":
      surveyName = a
    if o == "--icon":
      iconURL = a

  if len(args) != 2:
    print """
    usage: %s [options] input.fit output.kml

    Options:
    --orderby-field <string>  Choose a field in your FITs file to order
                              your placemarks

    --header <number>         Choose a FITS header to parse (default=1)

    --survey <string>         Specify a survey name to use when naming
                              objects (default='MySurvey')

    --icon <stringURL>        URL to use for placemark icons.
    """ % os.path.basename(sys.argv[0])
    sys.exit(2)

  fitsFile = args[0]
  kmlFile = args[1]


  fitsObjectList = fits.fitsparse.FITSParse(fitsFile,orderbyField,
                                            header,surveyName)


  print 'Making kml file...'

  kmlFileLink = open(kmlFile,'w')
  kmlFileLink.write('<?xml version="1.0" encoding="UTF-8"?>\n')
  kmlFileLink.write('<kml xmlns="http://earth.google.com/kml/2.2" '+
                    'hint="target=sky">\n')
  kmlFileLink.write('<Document>\n')
  kmlFileLink.write('<name>'+kmlFile.split('.')[0]+'</name>\n')
  kmlFileLink.write('  <Style id="FitsPoint">\n')
  kmlFileLink.write('    <IconStyle>\n')
  kmlFileLink.write('      <scale>1.0</scale>\n')
  kmlFileLink.write('      <Icon>\n')
  kmlFileLink.write('        <href>'+iconURL+'</href>\n')
  kmlFileLink.write('      </Icon>\n')
  kmlFileLink.write('    </IconStyle>\n')
  kmlFileLink.write('    <BalloonStyle>\n')
  kmlFileLink.write('      <text><![CDATA[<center>$[description]</center>]]></text>\n')
  kmlFileLink.write('      <color>ffffffff</color>\n')
  kmlFileLink.write('    </BalloonStyle>\n')
  kmlFileLink.write('  </Style>\n')

  for fitsPoint in fitsObjectList:
    kmlFileLink.write(fitsPoint.Placemark())

  kmlFileLink.write('</Document>\n')
  kmlFileLink.write('</kml>\n')

  kmlFileLink.close()

if __name__ == "__main__":
    main()
