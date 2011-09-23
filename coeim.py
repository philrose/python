from numpy import *
#import Image
from PIL import Image

def minmax(x, range=None):
    if range:
        lo, hi = range
        good = between(lo, x, hi)
        x = compress(good, x)
    return min(x), max(x)

def scale255minmax(data):
    lo, hi = minmax(ravel(data))
    scaled = (data - lo) / float(hi - lo) * 255
    return scaled.astype(int)

def savergb1(rgb, outfile):  # SLOW!
    #R, G, B = rgb
    #data = 256 * array([R, G, B])
    ny, nx = rgb[0].shape
    data = array(rgb).astype(int)
    data = transpose(data, (1,2,0))
    data.shape = (ny*nx,3)
    
    datal = data.tolist()
    datat = tuple([ tuple(sublist) for sublist in datal ])
    
    im = Image.new('RGB', (nx, ny))
    im.putdata(datat)
    im = im.transpose(Image.FLIP_TOP_BOTTOM)
    #im.show()
    im.save(outfile)

# r, g, b = data  (3, ny, nx)
def rgb2im(rgb):
    rgb = array(rgb).astype(uint8)
    rgb = transpose(rgb, (1,2,0))
    im = Image.fromarray(rgb, 'RGB')
    im = im.transpose(Image.FLIP_TOP_BOTTOM)
    return im

# r, g, b = data  (3, ny, nx)
def savergb(rgb, outfile):
    im = rgb2im(rgb)
    im.save(outfile)
    return im

def savegray(data, outfile, scale=False):
    ny, nx = data.shape
    im = Image.new('L', (nx,ny))
    if scale:
        data = scale255minmax(data)
    
    im.putdata(data.ravel())
    im = im.transpose(Image.FLIP_TOP_BOTTOM)
    im.save(outfile)
    return im

def loadrgb(infile):
    im = Image.open(infile)
    im = im.transpose(Image.FLIP_TOP_BOTTOM)
    # rgb = array(im.getdata())
    rgb = asarray(im)  # numpy
    print rgb.shape
    #nx, ny = im.size
    #rgb.shape = (ny,nx,3)
    rgb = transpose(rgb, (2,0,1))
    rgb = rgb[:3]  # in case there's an alpha channel on the end
    rgb.flags.writeable = True  # DEFAULT IS CAN'T EDIT IT!
    return rgb

def loadgray(infile):
    """Load grayscale image"""
    im = Image.open(infile)
    im = im.transpose(Image.FLIP_TOP_BOTTOM)
    data = asarray(im)  # numpy
    data.flags.writeable = True  # DEFAULT IS CAN'T EDIT IT!
    return data
