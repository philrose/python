#!/usr/bin/env python

# Class for performing least squares fits
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
# 3/12/07   Fixed bug in sigmaIterate() where the perpendicular offset distance
#           wasn't being used when deciding which points to reject when
#           perp_offset=True.  Wrote functions distance() and perp_distance()
#           to fix this problem.  This has a small (10^-3) effect on the final
#           fit slope.
#
# 3/9/07    Fixed bug in sigmaIterate() where the max_sigma parameter was
#           being ignored (it was always 3).  Fixed bugs in how sigmaIterate()
#           counted the number of rejections.  Now the loop properly terminates
#           when all or none of the points are rejected.  Added an optional
#           initial_fit parameter to sigmaIterate() that allows for
#           specification of the first fit line to use when initially
#           rejecting points.

"""
The PointArray class is designed for iterative least squares fitting routines
where points must be rejected on each pass.  It makes it easy to perform
these fits and to report on which points were actually used in the fitting.

Example Usage:

num_points = 100
x = []
y = []
for i in xrange(num_points - 5):
    x.append(float(i))
    y.append(float(i + 1))

# create 5 bad points
x.append(15.2)
y.append(65.8)

x.append(56.7)
y.append(14.6)

x.append(23.5)
y.append(67.8)

x.append(12.1)
y.append(30.0)

x.append(4.0)
y.append(50.0)

# create a PointArray
p = pointarray.PointArray(x, y, min_err=1.0e-4)
print "Points:\n%s" % p

fit = p.sigmaIterate()
print "Iterative least squares fit: %s" % fit   
print "Rejected points:"
for pt in p.rejectedPoints():
    print pt
"""

__author__ = "Jeremy Brewer (jeremy.d.brewer@gmail.com)"
__copyright__ = "Copyright 2005, 2006, 2007 Jeremy Brewer"
__license__ = "BSD"
__version__ = "1.0"

import math
import line

def distance(point, line):
    
    """
    Returns the delta y offset distance between a given point (x, y) and a
    line evaluated at the x value of the point.
    """
    
    return point.y - line(point.x)
    
def perp_distance(point, line):
    
    """
    Returns the perpendicular offset distance between a given point (x, y)
    and a line evaluated at the x value of the point.
    """
    
    return (point.y - line(point.x)) / math.sqrt(1.0 + line.slope ** 2)

class Point(object):

    """Class for representing points"""
    
    def __init__(self, x, y, y_err=0.0):
        self.x = x
        self.y = y
        self.y_err = y_err
        self.isRejected = False

    def __str__(self):
        if self.isRejected:
            s = "(%f, %f +- %f) R" % (self.x, self.y, self.y_err)
        else:
            s = "(%f, %f +- %f)" % (self.x, self.y, self.y_err)
        return s

    def allow(self):
        self.isRejected = False

    def reject(self):
        self.isRejected = True

class PointArray(object):

    """Class for arrays of points"""
    
    def __init__(self, x, y, y_err=None, min_err=1.0e-14):

        # check args
        assert len(x) == len(y)
        assert min_err > 0.0

        if y_err:
            assert len(x) == len(y_err)

            # build a list of y errors that is filtered for the minimum err
            err = []
            for e in y_err:
                if e >= min_err:
                    err.append(e)
                else:
                    err.append(min_err)
        else:
            err = [min_err] * len(x)

        # create a list of points
        self.points = map(Point, x, y, err)

    def __str__(self):
        return "\n".join(map(str, self.points))

    def __len__(self):
        return len(self.points)

    def __iter__(self):
        """Allows iteration in for statements"""
        for point in self.points:
            yield point

    def allowedPoints(self):
        """Allows iteration over only allowed points"""
        for pt in self.points:
            if not pt.isRejected:
                yield pt

    def rejectedPoints(self):
        """Allows iteration over only rejected points"""
        for pt in self.points:
            if pt.isRejected:
                yield pt

    def __getitem__(self, i):
        return self.points[i]

    def __setitem__(self, i, value):
        if not isinstance(value, Point):
            raise TypeError("object to set is not a Point")
        self.points[i] = value

    def allowAll(self):
    
        """
        Resets all points to be allowed.
        """

        for pt in self.rejectedPoints():
            pt.isRejected = False

    def leastSquaresFit(self):
        
        """
        Performs a least squares fit on the input data.
        """
        
        S = 0.0
        Sx = 0.0
        Sy = 0.0
        Sxy = 0.0
        Sxx = 0.0
        Syy = 0.0
        
        # compute sums
        for pt in self.allowedPoints():
            variance = pt.y_err ** 2
            S += 1.0 / variance
            Sx += pt.x / variance
            Sy += pt.y / variance
            Sxx += (pt.x ** 2) / variance
            Syy += (pt.y ** 2) / variance
            Sxy += (pt.x * pt.y) / variance

        # check for all points rejected
        if S == 0.0:
            return line.LinearFit()

        # compute the slope using a technique to minimize roundoff (see
        # Numerical Recipes for details)
        Stt = 0.0
        slope = 0.0

        for pt in self.allowedPoints():
            t = (pt.x - Sx / S) / pt.y_err
            Stt += t ** 2
            slope += (t * pt.y) / pt.y_err 

        slope /= Stt
        intercept = (Sy - Sx * slope) / S
        r2 = (Sxy * S - Sx * Sy) / \
             math.sqrt((S * Sxx - Sx * Sx) * (S * Syy - Sy * Sy))
              
        return line.LinearFit(slope, intercept, r2)

    def perpLeastSquaresFit(self):

        """
        Performs a perpendicular offset least squares fit on the input data.
        """

        S = 0.0
        Sx = 0.0
        Sy = 0.0
        Sxy = 0.0
        Sxx = 0.0
        Syy = 0.0
        
        # compute sums
        for pt in self.allowedPoints():
            variance = pt.y_err ** 2
            S += 1.0 / variance
            Sx += pt.x / variance
            Sy += pt.y / variance
            Sxx += (pt.x ** 2) / variance
            Syy += (pt.y ** 2) / variance
            Sxy += (pt.x * pt.y) / variance

        # check for all points rejected
        if S == 0.0:
            return line.LinearFit()

        B = ((S * Syy - Sy * Sy) - (S * Sxx - Sx * Sx)) / \
            (2.0 * (Sx * Sy - Sxy * S))

        # there are two solutions for the slope
        m = -B + math.sqrt(B * B + 1.0)
        m2 = -B - math.sqrt(B * B + 1.0)
        r2 = (Sxy * S - Sx * Sy) / \
                math.sqrt((S * Sxx - Sx * Sx) * (S * Syy - Sy * Sy))

        # slope for regular least squares fitting
        delta = S * Sxx - Sx * Sx
        ls_slope = (S * Sxy - Sx * Sy) / delta

        # choose the slope that is closest to normal least squares fitting
        diff = abs(ls_slope - m)
        diff2 = abs(ls_slope - m2)

        if diff <= diff2:
            slope = m
        else:
            slope = m2
            
        intercept = (Sy - slope * Sx) / S
        return line.LinearFit(slope, intercept, r2)

    def stddev(self, fit):
    
        """
        Returns the standard deviation of the difference from the input fit
        line for allowed points.  Returns -1 if there are too few allowed
        points to compute the standard deviation.
        """

        count = 0
        variance = 0.0

        for pt in self.allowedPoints():
            # take <y> as the value of the fit point at each point x
            delta = pt.y - fit(pt.x)
            variance += delta ** 2
            count += 1

        if count <= 1:
            return -1.0
        else:
            # the variance is bias corrected
            variance /= float(count - 1)
            sigma = math.sqrt(variance)
            return sigma

    def perpStddev(self, fit):

        """
        Returns the standard deviation of the perpendicular offset difference
        from the input fit line for allowed points.  Returns -1 if there are
        too few allowed points to compute the standard deviation.
        """

        count = 0
        variance = 0.0

        for pt in self.allowedPoints():
            # take <y> as the value of the fit point at each point x
            delta = (pt.y - fit(pt.x)) / math.sqrt(1.0 + fit.slope ** 2)
            variance += delta ** 2
            count += 1

        if count <= 1:
            return -1.0
        else:
            # the variance is bias corrected
            variance /= float(count - 1)
            sigma = math.sqrt(variance)
            return sigma

    def sigmaIterate(self, max_iter=5, max_sigma=3.0, perp_offset=False,
                     initial_fit=None):
        
        """
        Performs an iterative sigma clipping fit on the input data set.  A
        least squares line is fit to the data, then points further than
        max_sigma stddev away will be thrown out and the process repeated
        until no more points are rejected, all points are rejected, or the
        maximum number of iterations is passed.
        
        The sigma clipping algorithm uses either standard or perpendicular
        offset least squares fitting depending on the perp_offset flag.  By
        default normal least squares fitting is used.
        
        An optional initial fit to use can be supplied via initial_fit.  This
        is useful for noisy data where the "true" fit is approximately known
        beforehand.
        """

        # determine fit algorithm to use
        if perp_offset:
            lsFit = self.perpLeastSquaresFit
            stddev = self.perpStddev
            dist = perp_distance
        else:
            lsFit = self.leastSquaresFit
            stddev = self.stddev
            dist = distance
        
        # total number of rejected points
        total = 0
        
        for pt in self.points:
            if pt.isRejected:
                total += 1
        
        # initial fit
        if initial_fit is None:
            fit = lsFit()
        else:
            fit = initial_fit

        for i in xrange(max_iter):
            # standard deviation from fit line
            three_sigma = max_sigma * stddev(fit)
            
            # number of newly rejected points on each pass
            count = 0

            # throw away outliers
            for pt in self.allowedPoints():
                diff = abs(dist(pt, fit))
                if diff > three_sigma:
                    pt.isRejected = True
                    count += 1

            total += count

            # exit if none or all of the points were rejected
            if count == 0:
                break
            elif total == len(self.points):
                raise ValueError("all points were rejected")

            # update the fit
            fit = lsFit()

        return fit

if __name__ == "__main__":
    num_points = 100
    x = []
    y = []
    for i in xrange(num_points - 5):
        x.append(float(i))
        y.append(float(i + 1))

    # create 5 bad points
    x.append(15.2)
    y.append(65.8)

    x.append(56.7)
    y.append(14.6)

    x.append(23.5)
    y.append(67.8)

    x.append(12.1)
    y.append(30.0)

    x.append(4.0)
    y.append(50.0)

    # create a PointArray
    p = PointArray(x, y, min_err=1.0e-4)
    #print "Points:\n%s" % p

    fit = p.leastSquaresFit()
    print "Least squares fit: %s" % fit

    fit = p.perpLeastSquaresFit()
    print "Perp least squares fit: %s" % fit
    
    fit = p.sigmaIterate()
    print "Iterative least squares fit: %s" % fit   
    print "Rejected points:"
    for pt in p.rejectedPoints():
        print pt

    p.allowAll()
    
    fit = p.sigmaIterate(perp_offset=True)
    print "Iterative perp least squares fit: %s" % fit

    p.allowAll()
    
    fit = p.sigmaIterate(perp_offset=True, initial_fit=line.Line(1.0, 1.0))
    print "Iterative perp least squares fit: %s" % fit
