# im2rgbfits CL0024.png -over -header det.fits
# WILL HONOR WCS FROM headerfile

# im2rgbfits.py
# ~/ACS/CL0024/color/production/color.py
# ALSO SEE pyfits.pdf (Pyfits manual)

from coetools import *
from PIL import Image
import pyfits

def im2rgbfits(infile, rgbfile='', overwrite=False, headerfile=None, flip=False):
    if rgbfile == '':
        rgbfile = decapfile(infile) + '_RGB.fits'
        
    if exists(rgbfile):
        if overwrite:
            delfile(rgbfile)
        else:
            print rgbfile, 'EXISTS'
            sys.exit(1)
    
    #im = Image.open(infile)
    #print 'Loading data...'
    #data = array(im.getdata())
    #nxc, nyc = im.size
    #data.shape = (nyc,nxc,3)
    #data = transpose(data, (2,0,1))
    data = loadrgb(infile)
    
    #hdu = pyfits.PrimaryHDU()
    header = headerfile and pyfits.getheader(headerfile)
    hdu = pyfits.PrimaryHDU(None, header)
    hdulist = pyfits.HDUList([hdu])
    hdulist.writeto(rgbfile)

    try:  # If there's a 'SCI' extension, then that's where the WCS is
        header = pyfits.getheader(headerfile, 'SCI')
    except:
        pass
    
    if header <> None:
        if 'EXTNAME' in header.keys():
            del(header['EXTNAME'])
    
    for i in range(3):
        print 'RGB'[i]
        data1 = data[i]
        if flip:
            data1 = flipud(data1)
        pyfits.append(rgbfile, data1, header)
        
    print rgbfile, 'NOW READY FOR "Open RGB Fits Image" in IRAF'


if __name__ == '__main__':
    infile = sys.argv[1]
    
    outfile = ''
    if len(sys.argv) > 2:
        file2 = sys.argv[2]
        if file2[0] <> '-':
            outfile = file2
    
    params = params_cl()
    overwrite = 'over' in params.keys()
    headerfile = params.get('header', None)
    
    im2rgbfits(infile, outfile, overwrite=overwrite, headerfile=headerfile)


#hdulist = pyfits.open(rgbfile)
#hdulist.info()
