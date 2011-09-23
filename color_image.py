import pyfits
import sys
import numpy as num
import sigclip
import numarray.linear_algebra.mlab as mlab
import Image

# becker@astro.washington.edu
VERSION = 1.0

# globals
scale        = 2.5
scales       = [0.005 * scale, 0.0065 * scale, 0.0075 * scale] 
nonlinearity = 2.0  # 0 for linear, Inf for logrithmic.  Try 2.
skysub       = False

if __name__ == '__main__':
    rfits = sys.argv[1]
    gfits = sys.argv[2]
    bfits = sys.argv[3]
    outim = sys.argv[4]
    if len(sys.argv) == 6:
        skysub = True

    print 'Reading in data'
    rdata = pyfits.open(rfits)[0].data
    gdata = pyfits.open(gfits)[0].data
    bdata = pyfits.open(bfits)[0].data

    if skysub == True:
    
        print 'Clipping out bright/faint points'
        rdata_clipped = sigclip.sigclip(rdata, med=1)[0]
        gdata_clipped = sigclip.sigclip(gdata, med=1)[0]
        bdata_clipped = sigclip.sigclip(bdata, med=1)[0]
    
        print 'Estimating sky (mode)'
        rmode         = 2.5 * mlab.median(rdata_clipped) - 1.5 * mlab.mean(rdata_clipped)
        gmode         = 2.5 * mlab.median(gdata_clipped) - 1.5 * mlab.mean(gdata_clipped)
        bmode         = 2.5 * mlab.median(bdata_clipped) - 1.5 * mlab.mean(bdata_clipped)
        print '   R = %.2f; G = %.2f; B = %2f' % (rmode, gmode, bmode)

        print 'Subtracting sky'
        rdata = rdata - rmode
        gdata = gdata - gmode
        bdata = bdata - bmode


    print 'Scaling data (multiplicative)'
    rdata = rdata * scales[0]
    gdata = gdata * scales[1]
    bdata = bdata * scales[2]

    print 'Scaling data (arcsinh)'
    total = rdata + gdata + bdata
    if nonlinearity == 0:
        scaling = total
    else:
        scaling = num.arcsinh(total*nonlinearity) / nonlinearity
    rdata = rdata * scaling / total
    gdata = gdata * scaling / total
    bdata = bdata * scaling / total

    print 'Making sure colors do not saturate'
    # get the maximum of all 3 data values at each pixel
    mdata = num.maximum(num.maximum(rdata, gdata), bdata)
    # if mdata < 1, make it equal to 1
    # if mdata > 1, keep it
    mdata = num.maximum(mdata, 1)
    # divide by max
    rdata = rdata / mdata
    gdata = gdata / mdata
    bdata = bdata / mdata

    print 'Scaling to 8 bit image'
    rdata = num.minimum(num.floor(num.maximum(rdata * 256., 0)), 255)
    gdata = num.minimum(num.floor(num.maximum(gdata * 256., 0)), 255)
    bdata = num.minimum(num.floor(num.maximum(bdata * 256., 0)), 255)

    print 'Writing image', outim
    R = Image.frombuffer('L', (rdata.shape[1], rdata.shape[0]), rdata.astype(num.Int8))
    G = Image.frombuffer('L', (gdata.shape[1], gdata.shape[0]), gdata.astype(num.Int8))
    B = Image.frombuffer('L', (bdata.shape[1], bdata.shape[0]), bdata.astype(num.Int8))

    im = Image.merge('RGB', (R, G, B))
    im.save(outim)

