#!/usr/bin/env python

# Library for reading the astronomy community's FITS format
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
# 10/19/07 Removed NAXIS syncing code when updating NAXIS1, NAXIS2, etc. since
#          it was buggy.  Improved code for writing EXTEND keyword.
#
# 10/17/07 Re-added Fits class to improve backwards compatibility.  The Fits
#          class is deprecated in favor of pyfits but we should try to maintain
#          some compatibility with the older codebase.  Fixed bug when adding
#          NAXISx keywords.  Migrated to numpy from Numeric.
#
# 9/25/07  Stripped out Fits class and related code for pyfits migration,
#          leaving the utility functions that supplement pyfits.  Wrote
#          fits_simple_verify() that should catch most invalid FITS file cases.
#
# 8/17/07  Added keyword arguments to Fits() that allow for specification of
#          how to handle various minor errors.  The first supported keyword is
#          allow_repeat_keywords which causes read_header() to ignore keywords
#          after they have been encountered previously instead of raising an
#          error.  Fixed bugs where special header keywords weren't being
#          updated properly.
#
# 7/24/07  Fixed bug in read_data() where the BITPIX value wasn't being updated
#          to the proper value.
#
# 6/17/07  Fixed bug in find_hdus() where header blocks consisting on a single
#          card with END were marked as empty.

"""
Library for reading data from the astronomy FITS file format

The API in this module is taken from pyfits and is quite simple.  It is best
illustrated with an example:

# opens foo.fits and parses to find the location of each hdu
fits = fitslib.Fits("foo.fits")

# returns first header, which acts just like a Python dictionary
hdr = fits[0].header

# examples of reading keywords and comments
naxis1 = hdr["NAXIS1"]
naxis1_comment = hdr.getComment("NAXIS1", "Default value")
naxis2 = hdr.get("NAXIS2, 0)

# examples of writing keywords
hdr["NAXIS1"] = 20
hdr.update("NAXIS2", 10, comment="Comment for NAXIS2")

# prints the header as FITS compatible format
print "Header 1:\n%s" % hdr

# Returns data in first hdu as a numpy array
data = fits[0].data

# write modifications to file
fits.writeTo("bar.fits")

# closes underlying file object
fits.close()

This module grew out of my work on a WCS correction service where I needed to
inspect odd FITS files that caused pyfits trouble.  It is currently deprecated
in favor of pyfits, but I left the original code in place to avoid having
multiple versions of this module around (since I still use it much of my
personal code).  The main function still used by this module is
fits_simple_verify(), which does 2 quick checks for FITS file validity designed
to catch use cases where users accidentally pass in a non-FITS file.

Limitations compared to pyfits:
- Header keyword parsing / writing probably not 100% correct, but close
- Comments and other non-value keywords are currently just ignored
- More FITS verification needed
- Code for creating FITS images is young and probably a little buggy
"""

__author__ = "Jeremy Brewer (jeremy.d.brewer@gmail.com)"
__copyright__ = "Copyright 2005, 2006, 2007 Jeremy Brewer"
__license__ = "BSD"
__version__ = "1.06"

import os
import sys
import numpy

# mappings from BITPIX values to array typecodes and back
TYPECODES = {8: "B", 16: "h", 32: "l", -32: "f", -64: "d"}
BITPIXES = {}
for key in TYPECODES:
    value = TYPECODES[key]
    BITPIXES[value] = key
del key
del value

# constants
FITS_BLOCK_SIZE = 2880
FITS_CARD_SIZE = 80
FITS_KEY_SIZE = 8
MAX_STR_LENGTH = 70

def fits_simple_verify(fitsfile):
    
    """
    Performs 2 simple checks on the input fitsfile, which is a string
    containing a path to a FITS file.  First, it checks that the first card is
    SIMPLE, and second it checks that the file 2880 byte aligned.
    
    This function is useful for performing quick verification of FITS files.
    
    Raises:
      ValueError:  if either of the 2 checks fails
      IOError:     if fitsfile doesn't exist
    """
    
    if not os.path.exists(fitsfile):
        raise IOError("file '%s' doesn't exist" % fitsfile)

    f = open(fitsfile)

    try:
        # check first card name
        card = f.read(len("SIMPLE"))
        if card != "SIMPLE":
            raise ValueError("input file is not a FITS file")

        # check file size
        stat_result = os.stat(fitsfile)
        file_size = stat_result.st_size
        if file_size % FITS_BLOCK_SIZE != 0:
            raise ValueError("FITS file is not 2880 byte aligned (corrupted?)")
    finally:
        f.close()

def convert_value(value):

    """
    Converts the input value (as a string) to its native type.
    """

    # strings are treated separately in parse_card()
    if value == "T":
        value = True
    elif value == "F":
        value = False
    elif "." in value:
        value = float(value)
    else:
        value = int(value)

    return value

def parse_card(card):

    """
    Parses an individual header card.  The value of the card is interpreted
    as string, boolean, int, or float and returned as the appropriate type of
    data.  The if card is blank or a comment (name = COMMENT), then the card
    is ignored and None is returned for each return value.
    """

    name = card[0:8].rstrip()
    value_indicator = card[8:10]

    # check if a value is present
    if value_indicator == "= ":
        field_str = card[10:]

        if field_str.startswith("'"):
            # strings are treated separately because they can contain the
            # comment character /
            i = field_str[1:].find("'")

            # split occurs after removing the string component
            fields = field_str[1:i + 1].split("/", 1)

            if len(fields) == 1:
                comment = None
            elif len(fields) == 2:
                comment = fields[1].strip()
            else:
                raise ValueError("invalid card '%s'" % card)

            # +1 because i is relative to index 1
            value = field_str[1:i + 1]

            # leading spaces are significant, trailing spaces are not, but we
            # need to deal with the special case where the string is blank
            if value.isspace():
                value = " "
            else:
                value = value.rstrip()

        else:
            # non-string cases
            fields = field_str.split("/", 1)

            if len(fields) == 1:
                value = fields[0].strip()
                comment = None
            elif len(fields) == 2:
                value = fields[0].strip()
                comment = fields[1].strip()
            else:
                raise ValueError("invalid card '%s'" % card)

            # convert the value field to its intrinsic type
            value = convert_value(value)
    else:
        # ignore comments, history, and blank lines
        name = None
        value = None
        comment = None

    return name, value, comment

def read_header(fitsfile, offset=0, **kwargs):

    """
    Reads a header of a FITs file and returns 2 dictionaries, one for the
    values and one for the comments, and a list of the keys (needed to
    preserve the original order).

    The input fitsfile can be any file object that supports read() and
    seek().
    """

    # normally an error is raised for repeat keywords
    allow_repeat_keywords = kwargs.get("allow_repeat_keywords", False)

    values = {}
    comments = {}
    keys = []
    fitsfile.seek(offset)

    # verify that the file is a FITs file  by reading the first keyword
    card = fitsfile.read(FITS_CARD_SIZE)
    if not card:
        raise IOError("cannot read FITs header")

    name, value, comment = parse_card(card)
    values[name] = value
    comments[name] = comment
    keys.append(name)

    # check that the header is a valid primary header or extension
    if offset == 0:
        if name != "SIMPLE":
            raise ValueError("invalid FITs file: 1st keyword is not SIMPLE")
    else:
        if name != "XTENSION":
            raise ValueError("invalid FITs file: 1st keyword is not XTENSION")

    # read the header one card at a time (note 2880 / 80 = 36)
    while True:
        card = fitsfile.read(FITS_CARD_SIZE)
        if not card or card.startswith("END"):
            break

        name, value, comment = parse_card(card)

        if name is None:
            continue
        else:
            if name not in values:
                values[name] = value
                comments[name] = comment
                keys.append(name)
            else:
                if allow_repeat_keywords:
                    pass
                else:
                    raise ValueError("keyword '%s' appears more than once" % \
                                     name)

    return values, comments, keys

def read_data(fitsfile, hdr, offset, transform=True):

    """
    Reads data from a FITS image at the given location (offset) with the
    corresponding header hdr, which must support dictionary like mapping of
    FITS header keywords and a get() method.
    
    If transform is True, then the data will be transformed as
    data = BSCALE * data + BZERO
    
    and the necessary keywords in hdr will be updated.
    """

    # start at the proper location
    fitsfile.seek(offset)

    # currently only images are supported
    ext_type = hdr.get("XTENSION", "IMAGE")
    if ext_type != "IMAGE":
        raise ValueError("only image extensions are supported (type = %s)" % \
                         ext_type)

    # read keywords needed for determinig type and size of data to read
    try:
        bitpix = hdr["BITPIX"]
        naxis = hdr["NAXIS"]
        # TODO: use gcount and pcount for non-image extensions
        #gcount = hdr.get("GCOUNT", 1)
        #pcount = hdr.get("PCOUNT", 0)
        bscale = hdr.get("BSCALE", 1.0)
        bzero = hdr.get("BZERO", 0.0)
        naxis_list = []
        for i in xrange(1, naxis + 1):
            key = "NAXIS%d" % i
            n = hdr[key]
            naxis_list.append(n)
    except KeyError, key:
        raise ValueError("input header has no key %s" % key)

    # determine the type of array
    typecode = TYPECODES.get(bitpix, None)
    if typecode is None:
        raise ValueError("unsupported BITPIX value %s" % bitpix)

    # total number of elements to read
    num_elems = 1
    for n in naxis_list:
        num_elems *= n

    # read data
    data = numpy.fromfile(fitsfile, count=num_elems,
                          dtype=numpy.dtype(typecode))

    # FITS data is big endian
    if sys.byteorder == "little":
        data = data.byteswap(True)

    if transform:
        # apply linear transform (often used so 16 bits is enough)
        data = data * bscale + bzero

        # update header keywords (transform changes array type)
        hdr["BITPIX"] = BITPIXES[data.dtype.char]
        hdr["BSCALE"] = 1.0
        hdr["BZERO"] = 0.0

    # apply proper shape
    data.shape = naxis_list

    return data

def find_hdus(fitsfile):

    """
    Parses a FITs file to find the list of the byte locations of all HDUs.
    """

    # the file should be opened in mode 'rb'
    if fitsfile.mode != "rb":
        raise IOError("input FITs file must be opened in mode 'rb'")

    # hdr gives the location of the last header found, and all FITs files
    # must have a primary header
    hdr = 0
    hdus = []

    # current location in FITs file and whether a header or data block is
    # being read currently
    location = 0
    reading_hdr = True

    while True:
        # look at the next block
        location += FITS_BLOCK_SIZE

        if reading_hdr:
            # read current block in reverse, card by card, until the first
            # non-space card is found
            card_location = location - FITS_CARD_SIZE
            block_start = location - FITS_BLOCK_SIZE

            # read until one non-space card is found
            while True:
                # read one card
                fitsfile.seek(card_location)
                card = fitsfile.read(FITS_CARD_SIZE)
                if not card:
                    raise EOFError("EOF encountered before end of header")

                # remove padding
                card = card.rstrip()
                if card:
                    # found a card
                    break

                # skip back by one card
                card_location -= FITS_CARD_SIZE
                if card_location < block_start:
                    raise IOError("block at %d is all space" % block_start)

            if card == "END":
                # header has ended, determine if the next block is another
                # header or data block
                fitsfile.seek(location)
                key = fitsfile.read(FITS_KEY_SIZE)

                if not key:
                    # EOF reached
                    hdus.append((hdr, None))
                    break
                elif key == "XTENSION":
                    # found header
                    hdus.append((hdr, None))
                    hdr = location
                else:
                    # found data
                    hdus.append((hdr, location))
                    reading_hdr = False

        else:
            # read only a single byte
            fitsfile.seek(location - 1)
            byte = fitsfile.read(1)
            if not byte:
                raise EOFError("EOF encountered before end of data")

            if byte == "\0":
                # data may have ended
                key = fitsfile.read(FITS_KEY_SIZE)
                if not key:
                    # EOF reached
                    break
                elif key == "XTENSION":
                    # found header
                    hdr = location
                    reading_hdr = True

    return hdus

def pad_length(size):

    """
    Returns the number of bytes needed to pad the input integer to be a
    multiple of the FITS block size FITS_BLOCK_SIZE.
    """

    # outer modulus needed for when size is a multiple of FITS_BLOCK_SIZE
    return (FITS_BLOCK_SIZE - size % FITS_BLOCK_SIZE) % FITS_BLOCK_SIZE

def pad_str(s, pad):

    """
    Pads the input string using the supplied pad character to be a multiple
    of FITS_BLOCK_SIZE.  If the supplied pad string contains more than one
    character, only the first is used.
    """

    pad_len = pad_length(len(s))
    padding = pad[0] * pad_len

    return s + padding

def format_value(value):

    """
    Formats the value based on its type to be a FITS header compatible string.
    """

    if isinstance(value, str):
        s = "'%-8s'" % value
        val = "%-20s" % s
    elif isinstance(value, bool):
        if value:
            val = "T"
        else:
            val = "F"
    elif isinstance(value, int):
        val = "%20d" % value
    elif isinstance(value, float):
        if value >= 0.0:
            val = "%#20.14E" % value
        else:
            val = "%#20.13E" % value
    else:
        raise TypeError("invalid value type %s" % type(value).__name__)

    return val

def format_card(name, value, comment=None):

    """
    Formats a FITS header card given the name, value, and comment
    """

    val = format_value(value)

    if comment:
        # some strings are longer than 20, but have a maximum of 68 chars + 2
        # for the quotes
        if len(val) > 20:
            spaces = " " * (MAX_STR_LENGTH - len(val))
            card = "%-8s= %s%s" % (name, val, spaces)
        else:
            card = "%-8s= %20s / %-47s" % (name, val, comment)
    else:
        # long string case
        if len(val) > 20:
            spaces = " " * (MAX_STR_LENGTH - len(val))
            card = "%-8s= %s%s" % (name, val, spaces)
        else:
            card = "%-8s= %20s %-49s" % (name, val, " ")

    return card

def format_header(values, comments, keys, pad=False):

    """
    Formats a header to be FITS compatible, optionally padding the returned
    string to be a multiple of FITS_BLOCK_SIZE.  The headers are written in
    the order specified by keys.
    """

    # print keys in order
    h = []
    for key in keys:
        value = values[key]
        comment = comments.get(key, None)
        card = format_card(key, value, comment)
        h.append(card)

    h.append("%-80s" % "END")
    h = "".join(h)

    # add padding
    if pad:
        h = pad_str(h, " ")

    return h

def write_data(fitsfile, data, offset=None):

    """
    Writes array data to file with padding.  If offset is None, the data is
    written at the current position.  Otherwise it is written at the value
    of offset.
    """

    # seek to the proper position
    if offset is not None:
        fitsfile.seek(offset)

    # write the data as big endian
    if sys.byteorder == "little":
        copy = data.byteswap(False)
        data_str = copy.tostring()
    else:
        data_str = data.tostring()

    # TODO: write data in chunks so that huge arrays aren't copied into
    #       memory twice
    fitsfile.write(data_str)

    # pad the rest of the block with ascii nulls
    num_bytes = data.itemsize * numpy.size(data)
    padding = "\0" * pad_length(num_bytes)
    fitsfile.write(padding)

class Header(object):

    """
    Class that represents FITS headers, which consist of an ordered set of
    key value pairs.  This object behaves like an ordered dictionary where
    keywords are retained in the order they are given.

    Currently comment, history, and blank keywords are ignored and cannot be
    represented.
    """

    def __init__(self, values={}, comments={}, keys=[], name="",
                 extension=False):

        """
        Creates a new Header object, optionally building the object out of
        a dictionary of FITS keywords with values[name] = value,
        comments[name] = comment and a list of keyword names specifying the
        output order.  Both values and keys must be present and
        consistent (i.e. each key must be in values) to create a
        non-default Header object.
        """

        # check args
        if not isinstance(values, dict):
            raise TypeError("values is not of type dict")
        if not isinstance(comments, dict):
            raise TypeError("comments is not of type dict")
        if not isinstance (keys, (list, tuple)):
            raise TypeError("keys is not of type list or tuple")

        if values and keys:
            # check for compatibility between values, comments, and keys
            for key in keys:
                if key not in values:
                    raise KeyError("key '%s' not present in values" % key)

            for key in comments:
                if key not in values:
                    raise KeyError("comment '%s' not present in values" % key)

            self._values = values
            self._comments = comments
            self._keys = keys
            name = values.get("EXTNAME", "PRIMARY")
            self._name = name.upper()
        else:
            # default keywords
            if not extension:           
                self._values = {"SIMPLE": True, "BITPIX": 8, "NAXIS": 0,
                                "EXTEND": True}
                self._keys = ["SIMPLE", "BITPIX", "NAXIS", "EXTEND"]
            else:
                self._values = {"XTENSION": "Image", "BITPIX": 8, "NAXIS": 0}
                self._keys = ["XTENSION", "BITPIX", "NAXIS"]
            
            self._comments = {}
            
            # TODO: currently the extension name isn't used for newly
            #       created headers, although it can be set through EXTNAME
            if name:
                self._name = name
            else:
                self._name = "PRIMARY"

    def getName(self):
        
        """
        Returns the name of this Header.
        """
        
        return self._name
    
    name = property(fget=getName, doc=getName.__doc__)

    def __len__(self):

        """
        Returns the number of keywords contained in the header, not counting
        comment, history, and blank keywords.
        """

        return len(self._keys)

    def __getitem__(self, key):

        """
        Allows for dictionary like access of keyword values.  Only the value
        of each keyword is returned.  Raises a KeyError if they supplied
        keyword isn't found.
        """

        key = str(key).upper()
        value = self._values.get(key, None)
        if value is None:
            raise KeyError("keyword '%s' not present" % key)
        return value

    def __setitem__(self, key, value):

        """
        Sets the value for a keyword without changing the comment.  The
        update() method allows setting both the value and the comment for a
        keyword.
        """

        key = str(key).upper()

        if not isinstance(value, (bool, int, float, str)):
            raise TypeError("invalid type for keyword value %s" % \
                            type(value).__name__)

        old_value = self._values.get(key, None)

        # determine where to insert the new value, making special exceptions
        # for the important keywords
        if key in ("SIMPLE", "XTENSION"):
            # only needed so users can change extensions
            self._keys[0] = key
        elif key == "NAXIS":
            # only needed for special case below
            pass
        elif key.startswith("NAXIS"):
            # find the previous NAXIS keyword
            try:
                i = int(key[5:])
            except ValueError:
                raise KeyError("invalid NAXIS keyword '%s'" % key)

            if i == 1:
                # special case for NAXIS1, the first keyword
                j = 2
            else:           
                prev_key = "NAXIS%d" % (i - 1)
                try:
                    j = self._keys.index(prev_key)
                except ValueError:
                    raise KeyError("previous key '%s' not found" % prev_key)
            
            # only add the NAXISx keyword if not present
            if old_value is None:
                self._keys.insert(j + 1, key)
        elif key == "EXTEND":
            # find the last NAXIS keyword
            for i in xrange(999):
                naxis_key = "NAXIS%d" % i
                if naxis_key not in self._values:
                    break
            if old_value is None:
                self._keys.insert(i + 1, key)
        else:
            # for normal keywords, only change position if the keyword isn't
            # present
            if old_value is None:
                self._keys.append(key)

        # update the dictionary after any errors that could have occured have
        # been raised
        self._values[key] = value

    def __delitem__(self, key):

        """
        Removes a keyword from the header.  Raises a KeyError if the given
        keyword isn't present.
        """

        key = str(key).upper()

        if key not in self._values:
            raise KeyError("keyword '%s' not present" % key)

        # TODO: Insert checks to keep NAXIS up to date and to ensure that
        # NAXIS1 can't be deleted when NAXIS2 is present.  This might be
        # better done in a verify() method that checks that the user input
        # things in a sane way.  That way, it would be possible to clean up
        # incorrectly formatted headers.

        if key in ("SIMPLE", "EXTENSION", "NAXIS"):
            raise KeyError("cannot delete keyword '%s'" % key)

        del self._values[key]
        i = self._keys.index(key)
        del self._keys[i]

        if key in self._comments:
            del self._comments[key]

    def __iter__(self):

        """
        Allows iteration over the keys in a header in their specified order.
        """

        for key in self._keys:
            yield key

    def iterkeys(self):

        """
        Allows iteration over the keys in a header in their specified order
        (identical to __iter__()).
        """

        for key in self._keys:
            yield key

    def itervalues(self):
        
        """
        Allows iteration over the values in a header in the specified order
        of the keywords.
        """

        for key in self._keys:
            value = self._keys[key]
            yield value

    def iteritems(self):
        
        """
        Allows iteration over the keys and values in a header in the
        specified order of the keywords.
        """
        
        for key in self._keys:
            value = self._keys[key]
            yield key, value

    def keys(self):
        
        """
        Returns a list of all the keywords in a header.
        """
        
        keylist = []
        
        for key in self.iterkeys():
            keylist.append(key)

        return keylist

    def values(self):
        
        """
        Returns a list of all the values in a header.
        """
        
        valuelist = []
        
        for value in self.itervalues():
            valuelist.append(value)
        
        return valuelist

    def items(self):
        
        """
        Returns a list of keyword / value pairs for the header.
        """
        
        itemlist = []
        
        for key, value in self.iteritems():
            itemlist.append((key, value))

        return itemlist
       
    def get(self, key, default=None):

        """
        Returns the specified keyword if present, otherwise the supplied
        default, which is None by default.  Analagous to the get() method of
        dictionary objects.
        """

        key = str(key).upper()
        return self._values.get(key, default)

    def getComment(self, key, default=None):

        """
        Like get(), but returns the comment instead of the value using the
        supplied default if the keyword isn't present.
        """

        key = str(key).upper()
        return self._comments.get(key, default)

    def update(self, key, value, comment=None):

        """
        Updates a keyword if it is present or adds a new keyword if not
        present.  This method also allows for updating comments.
        """

        key = str(key).upper()
        self.__setitem__(key, value)
        self._comments[key] = comment

    def clearNaxisKeys(self):
        
        """
        Sets NAXIS = 0 and removes all previous NAXISx keywords.
        """
        
        self._values["NAXIS"] = 0
        for i in xrange(999):
            naxis_key = "NAXIS%d" % i
            if naxis_key in self._values:
                self.__delitem__(naxis_key)
            else:
                break

    def has_key(self, key):

        """
        Returns whether the header contains the given keyword.
        """

        key = str(key).upper()
        return key in self._values

    def __contains__(self, key):

        """
        Returns whether the header contains the given keyword (identical to
        the has_key() method).
        """

        key = str(key).upper()
        return key in self._values

    def __str__(self):

        """
        Returns a non-padded string representation of a FITS header.  This
        method is identical to calling the format_header() function with
        pad=False.
        """

        return format_header(self._values, self._comments, self._keys,
                             pad=False)

    def format(self):

        """
        Returns a padded string representation of a FITS header.  This
        method is identical to calling the format_header() function with
        pad=True.
        """

        return format_header(self._values, self._comments, self._keys,
                             pad=True)

class HDU(object):

    """
    Class for representing HDUs (Header Data Units) in FITS files.
    """

    def __init__(self, fits=None, hdr_pos=None, data_pos=None,
                 extension=False):

        self._fits = fits
        self._hdr_pos = hdr_pos
        self._data_pos = data_pos
        self._hdr = None
        self._data = None

        if fits is not None and hdr_pos is not None:
            # read the header immediately
            self._readHeader()

            # the name of the HDU is given by the name in the header
            self._name = self._hdr.name
            
            # check whether this is an extension
            if hdr_pos == 0:
                self._is_extension = False
            else:
                self._is_extension = True
        else:
            # create an empty header
            self._hdr = Header(extension=extension)
            self._is_extension = extension

    def isPrimary(self):
        
        """
        Returns whether this is a primary extension.
        """
        
        return not self._is_extension

    def isExtension(self):

        """
        Returns whether this is a a non-primary extension.
        """

        return self._is_extension

    def hasData(self):

        """
        Returns whether the HDU has data in addition to a header.
        """

        return self._data_pos is not None

    def hasReadData(self):
        
        """
        Returns whether the HDU has read its associated data from file.
        """
        
        return self._data is not None

    def _readHeader(self):

        """
        Reads in the header for the current HDU and creates a Header object.
        There is no need to call this method separately because it is
        automatically invoked.
        """

        values, comments, keys = \
            read_header(self._fits._fitsfile, self._hdr_pos,
                        **self._fits._kwargs)
        self._hdr = Header(values, comments, keys)

    def readData(self, transform=True):

        """
        Reads in the data for the current HDU if present.  If there is no
        associated data, a ValueError is raised.
        
        By default, data from an HDU is read into memory only when requested.
        This method can be used to override this behavior.
        """

        if not self.hasData():
            raise ValueError("HDU does not contain a data section")

        self._data = read_data(self._fits._fitsfile, self._hdr,
                               self._data_pos, transform=transform)

    def getName(self):
        
        """
        Returns the name of the HDU, which is the same as the name of the
        header contained in this HDU.
        """
        
        return self._name

    name = property(fget=getName, doc=getName.__doc__)

    def getHeader(self):

        """
        Returns the header in this HDU.
        """

        return self._hdr

    def setHeader(self, hdr):
        
        """
        Sets the header in this HDU.
        """
        
        if not isinstance(hdr, Header):
            raise TypeError("trying to set the header to a non-Header object")
        
        self._hdr = hdr

    header = property(fget=getHeader, fset=setHeader,
                      doc="Access the header in this HDU")

    def getData(self, transform=True):

        """
        Returns the data in this HDU.  Raises a ValueError if there is no
        data present in the HDU.
        
        If transform is True, the data is transformed according to
        data = BSCALE * data + BZERO
        
        and the necessary keywords in the header are updated.
        """

        if not self.hasReadData():
            self.readData(transform=transform)

        return self._data

    def setData(self, data):
        
        """
        Sets the data in this HDU.
        """
        
        # TODO: the type of object numpy returns is unclear, so I can't
        #       test that here

        # check that the type is supported
        bitpix = BITPIXES.get(data.dtype.char, None)
        if bitpix is None:
            raise ValueError("unsupported numpy typecode %s" % data.dtype.char)

        # update the keywords in the header to be consistent        
        self._hdr["BITPIX"] = bitpix
        self._hdr["BSCALE"] = 1.0
        self._hdr["BZERO"] = 0.0
        self._hdr.clearNaxisKeys()
        self._hdr["NAXIS"] = len(data.shape)
        
        for i, n in enumerate(data.shape):
            key = "NAXIS%d" % (i + 1)
            self._hdr[key] = n

        # update internal variables to recognize that there is data
        self._data_pos = True
        self._data = data

    data = property(fget=getData, fset=setData,
                    doc="Access the data in this HDU")

class Fits(object):

    """
    Class for representing FITS files.
    """

    def __init__(self, fitsfile=None, **kwargs):

        """
        Opens or creates a FITS file object.  Newly created Fits objects
        contain a single HDU instance with a basic header.  Additional
        HDUs can be added through the append() method.

        Input:  fitsfile -- can be a filename or a file object that supports
                            mode, seek, read, and close, or None if the Fits
                            object is to be newly created
                kwargs   -- additional keywords arguments used for controlling
                            how to handle minor errors encountered when
                            reading FITS files (see below)
        
        
        Keyword arguments:

        allow_repeat_keywords -- don't raise an error if the same keyword is
                                 found multiple times in a header (the first
                                 value is the one used)
        """

        self._kwargs = kwargs

        if fitsfile is None:
            # create a new Fits object with a single HDU
            hdu = HDU(extension=False)
            self._hdus = [hdu]
            self._fitsfile = None
        else:
            # open and read the FITS file
            if isinstance(fitsfile, str):
                # fitsfile is a filename
                if not os.path.exists(fitsfile):
                    raise IOError("file '%s' doesn't exist" % fitsfile)

                self._fitsfile = open(fitsfile, "rb")
            else:
                # fitsfile is a file object
                if not hasattr(fitsfile, "close"):
                    raise TypeError("fitsfile lacks close method")
                elif not hasattr(fitsfile, "read"):
                    raise TypeError("fitsfile lacks read method")
                elif not hasattr(fitsfile, "seek"):
                    raise TypeError("fitsfile lacks seek method")
                elif not hasattr(fitsfile, "mode"):
                    raise TypeError("fitsfile lacks mode attribute")

                self._fitsfile = fitsfile

            # parse the file to find the locations of all HDUs
            hdu_positions = find_hdus(self._fitsfile)

            # create a list of HDU objects for each location (the headers
            # are read in, but the data is not)
            self._hdus = []
            for hdr_pos, data_pos in hdu_positions:
                hdu = HDU(self, hdr_pos, data_pos)
                self._hdus.append(hdu)

    def __del__(self):

        """
        Closes the internal file referenced by the Fits object.
        """

        if self._fitsfile:
            self._fitsfile.close()

    def close(self):

        """
        Closes the internal file object referenced by the Fits object.
        """

        if self._fitsfile:
            self._fitsfile.close()

    def __len__(self):

        """
        Returns the number of HDUs in the FITS file.
        """

        return len(self._hdus)

    def __iter__(self):

        """
        Allows for iteration over all HDUs.
        """

        for h in self._hdus:
            yield h

    def __getitem__(self, i):

        """
        Returns the ith HDU.
        """

        return self._hdus[i]

    def __setitem__(self, i, hdu):
        
        """
        Sets the ith HDU.
        """
        
        if not isinstance(hdu, HDU):
            raise TypeError("trying to set an HDU to non HDU class")
        
        self._hdus[i] = hdu

    def __delitem__(self, i):
        
        """
        Deletes the ith HDU.
        """
        
        del self._hdus[i]

    def getHDU(self, i):

        """
        Returns the ith HDU (identical to __getitem__()).
        """

        return self.__getitem__(i)

    def setHDU(self, i, hdu):
        
        """
        Sets the ith HDU (identical to __setitem__()).
        """
        
        self.__setitem__(i, hdu)

    def append(self, hdu):
        
        """
        Adds an additional HDU.
        """
        
        if not isinstance(hdu, HDU):
            raise TypeError("trying to set an HDU to non HDU class")

        self._hdus.append(hdu)

    def readAllData(self, transform=False):
        
        """
        Forces the entire FITS file to be read into memory.  Portions that
        have already been read will not be read again.
        """
        
        # read all data from the file before overwriting
        for hdu in self._hdus:
            if hdu.hasData() and not hdu.hasReadData():
                hdu.readData(transform=transform)

    def writeTo(self, filename, overwrite=False):

        """
        Writes a Fits object to file.  In addition, the remaining parts of
        the FITS file will be read into memory using readAllData() and the
        underlying file object will be closed using close() to ensure proper
        behavior in the case of overwriting the internal file referenced.
        """

        if os.path.exists(filename) and not overwrite:
            raise IOError("file '%s' already exists" % filename)

        # force reading of the entire file into memory
        self.readAllData(transform=False)

        # close internal file in case the file is the same
        self.close()

        f = open(filename, "wb")
        try:
            for hdu in self._hdus:
                # write header
                hdr = hdu.header
                hdr_str = hdr.format()
                f.write(hdr_str)

                # write data if present
                if hdu.hasData():
                    data = hdu.data
                    write_data(f, data)
        finally:
            f.close()

def main(argv):
    import time

    # check args
    if len(argv) != 2:
        print "Usage: %s <fits-file>" % os.path.basename(argv[0])
        sys.exit(2)

    if not os.path.exists(argv[1]):
        print "Error: file '%s' doesn't exist" % argv[1]
        sys.exit(1)

    # create a Fits object (parses file to find structure)
    fits = Fits(argv[1])

    # print out each header
    for i, hdu in enumerate(fits):
        hdr = hdu.header
        print "Header %d '%s':\n%s\n" % (i, hdu.name, hdr)

    print
    print "Testing reading of data:"
    start = time.time()
    data = fits[0].data
    stop = time.time()
    print "Read %d objects in %f sec" % (numpy.size(data), stop - start)
    print "Data type is %s" % data.dtype.char
    shape_str = map(str, data.shape)
    print "Data shape is (%s)" % ", ".join(shape_str)
    print "Altered keywords:"
    print "BITPIX = %s" % hdr["BITPIX"]
    print "BSCALE = %s" % hdr["BSCALE"]
    print "BZERO = %s" % hdr["BZERO"]

    # write the fits file to a copy
    fits.writeTo("test.fits", overwrite=True)
    print "Wrote copy to file 'test.fits'"

    fits.close()
    
    # now create a new Fits file from scratch
    fits = Fits()
    data = 255 * numpy.ones((200, 200), dtype=numpy.dtype("B"))
    hdr = Header()
    hdr["test"] = "A test"
    fits[0].header = hdr
    fits[0].data = data
    print "NAXIS = %d" % hdr["NAXIS"]
    fits.writeTo("test2.fits", overwrite=True)
    print "Wrote white 200 x 200 image to 'test2.fits'"

if __name__ == "__main__":
    main(sys.argv)

