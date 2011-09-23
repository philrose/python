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
Intelligently compresses tiles from a regionated image directory to JPEG.
"""

import os
import sys
import optparse
import Image

def main(argv):
  usage = "%prog [options] <regionated directory>"
  parser = optparse.OptionParser(usage)
  parser.add_option("--quality", "-q", default=60, type="int",
                    help="JPEG quality setting, default %default")
  flags, files = parser.parse_args(argv[1:])
  
  if len(files) != 1:
    print "usage: %s <regionated directory>" % os.path.basename(argv[0])
    sys.exit(2)

  path = files[0]

  if not os.path.exists(path):
    print "Error: directory '%s' doesn't exist"
    sys.exit(2)
  elif not os.path.isdir(path):
    print "Error: input '%s' is not a directory"
    sys.exit(2)

  all_files = os.listdir(path)
  os.chdir(path)
  print "Using JPEG quality %d" % flags.quality

  for kmlfile in all_files:
    name, ext = os.path.splitext(kmlfile)
    pngfile = "%s.png" % name
    
    if ext != ".kml":
      continue
    
    image = Image.open(pngfile)
    
    # Check for images without transparency.
    if image.mode in ("RGB", "L"):
      # Convert PNG to JPEG and remove PNG.
      jpgfile = "%s.jpg" % name
      print "Converting %s to %s" % (pngfile, jpgfile)
      image = Image.open(pngfile)
      image.save(jpgfile, quality=flags.quality)
      os.remove(pngfile)
      
      # Update the KML file to point to the proper file.
      f = open(kmlfile)
      kml = f.read()
      f.close()
      f = open(kmlfile, "w")
      f.write(kml.replace(".png", ".jpg"))
      f.close()

if __name__ == "__main__":
  main(sys.argv)
