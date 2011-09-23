#!/usr/bin/env python

# Class for lines
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
# 2/7/06    Updated documentation.

"""
Line class module

The Line class is a simple class for holding the slope and intercept of a
line.  Line objects are callable -- when called, they evaluate at the given
x value, e.g. y = line(x) gives the value of line at x.

Example Usage.:

l = line.Line(2.0, 3.0)
x1, x2, x3 = 3.14, 0.0, 1.0

print "The line is %s" % l
print "f(%f) = %f" % (x1, l(x1))
print "f(%f) = %f" % (x2, l(x2))
print "f(%f) = %f" % (x3, l(x3))

lp = l.perpToAtX(0.0)
print "Line perpendicular to original at x = 0 is %s" % lp
lp = l.perpToAtY(3.0)
print "Line perpendicular to original at y = 3 is %s" % lp

flip = l.flipXY()
print "Line with flipped x, y is %s" % flip

fit = line.LinearFit(2.0, 3.0, 0.987)
print "Linear fit is %s" % fit
"""

__author__ = "Jeremy Brewer (jeremy.d.brewer@gmail.com)"
__copyright__ = "Copyright 2005, 2006, 2007 Jeremy Brewer"
__license__ = "BSD"
__version__ = "1.0"

class Line(object):

    """Class for describing lines"""

    def __init__(self, slope=0.0, intercept=0.0):

        """
        Initializes a line to have the given slope and intercept.
        
        Input:  slope     -- slope of the line
                intercept -- intercept of the line
        """

        try:
            self.slope = float(slope)
        except ValueError:
            raise TypeError("invalid slope value '%s'" % slope)
        try:
            self.intercept = float(intercept)
        except ValueError:
            raise TypeError("invalid intercept value '%s'" % intercept)

    def __str__(self):
        """Returns a string representation of the line"""
        return "y = %f * x + %f" % (self.slope, self.intercept)

    def __call__(self, x):
        """Evaluates a line at a given position x"""
        assert isinstance(x, float)
        return self.slope * x + self.intercept

    def perpToAtX(self, x):
    
        """Returns a line perpendicular to this line at the given x position"""
        
        assert isinstance(x, float)
        perpSlope = -1.0 / self.slope
        perpIntercept = x * (self.slope + 1.0 / self.slope) + self.intercept
        return Line(perpSlope, perpIntercept)

    def perpToAtY(self, y):

        """Returns a line perpendicular to this line at the given x position"""

        assert isinstance(y, float)
        x = (y - self.intercept) / self.slope
        return self.perpToAtX(x)

    def flipXY(self):
        
        """Creates a line where x and y have been flipped"""
        
        if self.slope == 0.0:
            raise ZeroDivisionError("cannot flip line with slope = 0")

        newSlope = 1.0 / self.slope
        newIntercept = -self.intercept / self.slope
        return Line(newSlope, newIntercept)

class LinearFit(Line):

    """Class for describing linear fits"""
    
    def __init__(self, slope=0.0, intercept=0.0, r2=0.0):
    
        """
        Initializes a linear fit to have the given slope, intercept, and
        correlation coefficient.
        
        Input:  slope     -- slope of the line
                intercept -- intercept of the line
                r2        -- correlation coefficient of the line
        """

        Line.__init__(self, slope, intercept)

        try:
            self.r2 = float(r2)
        except ValueError:
            raise TypeError("invalid r2 value '%s'" % r2)

    def __str__(self):
        """Returns a string representation of the linear fit"""
        return "y = %f * x + %f, r^2 = %f" % (self.slope, self.intercept,
                                              self.r2)

# testing code
if __name__ == "__main__":
    l = Line(2.0, 3.0)
    x1, x2, x3 = 3.14, 0.0, 1.0

    print "The line is %s" % l
    print "f(%f) = %f" % (x1, l(x1))
    print "f(%f) = %f" % (x2, l(x2))
    print "f(%f) = %f" % (x3, l(x3))

    lp = l.perpToAtX(0.0)
    print "Line perpendicular to original at x = 0 is %s" % lp
    lp = l.perpToAtY(3.0)
    print "Line perpendicular to original at y = 3 is %s" % lp

    flip = l.flipXY()
    print "Line with flipped x, y is %s" % flip
    
    fit = LinearFit(2.0, 3.0, 0.987)
    print "Linear fit is %s" % fit
