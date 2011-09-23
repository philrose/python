import pyfits
import numpy as np
import re,os,subprocess,sys
import matplotlib.pyplot as plt
import matplotlib.nxutils as nxutils
from matplotlib import cm
from PadovaTracksUtils import *

def smooth(x,window_len=11,window='hanning'):
    """
    taken from http://www.scipy.org/Cookbook/SignalSmooth
    smooth the data using a window with requested size.
    
    This method is based on the convolution of a scaled window with the signal.
    The signal is prepared by introducing reflected copies of the signal 
    (with the window size) in both ends so that transient parts are minimized
    in the begining and end part of the output signal.
    
    input:
        x: the input signal 
        window_len: the dimension of the smoothing window; should be an odd integer
        window: the type of window from 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'
            flat window will produce a moving average smoothing.
    
    output:
        the smoothed signal
        
    example:
    
    t=linspace(-2,2,0.1)
    x=sin(t)+randn(len(t))*0.1
    y=smooth(x)
    
    see also: 
    
    np.hanning, np.hamming, np.bartlett, np.blackman, np.convolve
    scipy.signal.lfilter
 
    TODO: the window parameter could be the window itself if an array instead of a string   
    """
    
    if x.ndim != 1:
        raise ValueError, "smooth only accepts 1 dimension arrays."
    if x.size < window_len:
        raise ValueError, "Input vector needs to be bigger than window size."
    
    if window_len<3:
        return x
    
    if not window in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
        raise ValueError, "Window is on of 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'"
    
    s=np.r_[x[window_len-1:0:-1],x,x[-1:-window_len:-1]]
    #print(len(s))
    if window == 'flat': #moving average
        w=np.ones(window_len,'d')
    else:
        w=eval('np.'+window+'(window_len)')
    
    y=np.convolve(w/w.sum(),s,mode='valid')
    return y

def check_mkdir(outdir):
    if not os.path.isdir(outdir): 
        os.mkdir(outdir)
        print 'made %s'%outdir
    return os.path.abspath(outdir)

def read_array(filename, commentchar='#',separator=' '):
    """ Read a file with an arbitrary number of columns.
        The type of data in each column is arbitrary
        It will be cast to the given dtype at runtime
        mydescr = np.dtype([('column1', 'int32'), ('column2Name', 'uint32'), ('col3', 'uint64'), ('c4', 'float32')])
        myrecarray = read_array('file.csv', mydescr)
    """
    file = open(filename, 'r')
    firstline = file.readline()
    secondline = file.readline()
    file.close()
    col_names = firstline.replace('commentchar','').split()
    types = [type(is_numeric(s)) for s in secondline.split()]
    str_types = [str(t).split()[-1].split('>')[0].replace('\'','') for t in types]
    dtype = np.dtype([(col_name,str_type) for col_name,str_type in zip(col_names,str_types)])
    
    cast = np.cast
    data = [[] for dummy in xrange(len(dtype))]
    for line in open(filename, 'r'):
        fields = line.strip().split()
        for i, number in enumerate(fields):
            data[i].append(number)
    for i in xrange(len(dtype)):
        data[i] = cast[dtype[i]](data[i])
    return np.rec.array(data, dtype=dtype)

def call_polyfit(x,y,order):
    '''
    ex:
    fit = call_polyfit(color,mag2,order)
    plt.plot(color,fit,'o',color='red')
    '''
    z,residuals, rank, singular_values, rcond = np.polyfit(x,y,order,full=True)
    p = np.poly1d(z)
    print 'roots:',p.r
    print 'coeffs [0]x**n [1]x**n-1 etc:', p.c
    return p(x)

def normalize_to_slope(lowlim,highlim,arrOrig,arrtoNorm):
    NOrig = len(np.nonzero( (lowlim>arrOrig) & (highlim<arrOrig) )[0])
    NtoNorm = len(np.nonzero( (lowlim>arrtoNorm) & (highlim<arrtoNorm) )[0])
    print NOrig,NtoNorm
    return float(NOrig)/float(NtoNorm)

def continuous_hist(arr,bins,offset):
    '''
    bins is a np array.
    Must be equally sized bins.
    makes a histogram of arr with bins=bins then repeats with bins+offset
    returns a combined histogram and bins+offsets
    '''
    hist = []
    allbins = []
    for k in range(int(1./offset)):
        xbins = bins+offset*float(k)
        allbins.append(xbins[:-1])
        hist.append(np.histogram(arr,bins=xbins)[0])
        
    sort_bins = np.concatenate(np.transpose(allbins))
    sort_hist = np.concatenate(np.transpose(hist))
    return sort_hist,sort_bins

def mag2polygon(magmax,magmin=99,colmin=-99,colmax=99):
    '''
    makes a parallelogram out of four cmd points. 
    upper left -> upper right -> lower right -> lower left
    '''
    verts = np.array([(colmin,magmax),(colmax,magmax),(colmax,magmin),(colmin,magmin)])
    return verts
    
def mags2polygon(mag1,mag2):
    return np.column_stack((mag1-mag2,mag2))

def manual_hist(arr,bins):
    '''
    in: arr,bins
    out: hist,bins,inds
    '''
    inds,hist = [],[]
    for i in range(len(bins)-1):
        within = np.nonzero((arr < bins[i+1]) & (arr >= bins[i]))[0]
        # in case nothing there
        if len(within) == 0: 
            inds.append(np.ndarray([1])*0-1.)
            hist.append(0)
        else:
            inds.append(within)
            hist.append(len(within))
    return hist,bins,inds

def continuous_hist_manual(arr,bins,offset,flux,frac=False):
    '''
    the same as continuous_hist, except it returns the indices of array arr for each bin.
    flux is np array of the mags in question, eg, use vegamag2flux
    offset is float bin spacing for each histogram
    bins is a np array.
    Must be equally spaced bins.
    '''
    int_flxs = []
    frac_flxs = []
    allbins = []
    hist = []
    for k in range(int(1./offset)):
        xbins = bins+offset*float(k)
        allbins.append(xbins[:-1])
        hist.append(manual_hist(arr,xbins)[0])
        inds = manual_hist(arr,xbins)[2]
        int_flx=[]
        for j in range(len(inds)):
            if len(inds[j]) == 1: 
                int_flx.append(0)
            else: 
                int_flx.append(np.sum(flux[inds[j]]))
                
        frac_flux = [np.sum(flux[k:]) for k in range(len(inds))]
        frac_flxs.append(frac_flux)
        int_flxs.append(int_flx)
        
    sort_bins = np.concatenate(np.transpose(allbins))
    sort_hist = np.concatenate(np.transpose(hist))
    if frac==True: sort_int_flxs = np.concatenate(np.transpose(frac_flxs))
    else: sort_int_flxs = np.concatenate(np.transpose(int_flxs))
    
    return sort_hist,sort_bins,sort_int_flxs

# griddata.py - 2010-07-11 ccampo
def griddata(x, y, z, binsize=0.01, retbin=True, retloc=True):
    """
    take from 
    http://www.scipy.org/Cookbook/Matplotlib/Gridding_irregularly_spaced_data
    on 02/08/11
    Place unevenly spaced 2D data on a grid by 2D binning (nearest
    neighbor interpolation).
    
    Parameters
    ----------
    x : ndarray (1D)
        The idependent data x-axis of the grid.
    y : ndarray (1D)
        The idependent data y-axis of the grid.
    z : ndarray (1D)
        The dependent data in the form z = f(x,y).
    binsize : scalar, optional
        The full width and height of each bin on the grid.  If each
        bin is a cube, then this is the x and y dimension.  This is
        the step in both directions, x and y. Defaults to 0.01.
    retbin : boolean, optional
        Function returns `bins` variable (see below for description)
        if set to True.  Defaults to True.
    retloc : boolean, optional
        Function returns `wherebins` variable (see below for description)
        if set to True.  Defaults to True.
   
    Returns
    -------
    grid : ndarray (2D)
        The evenly gridded data.  The value of each cell is the median
        value of the contents of the bin.
    bins : ndarray (2D)
        A grid the same shape as `grid`, except the value of each cell
        is the number of points in that bin.  Returns only if
        `retbin` is set to True.
    wherebin : list (2D)
        A 2D list the same shape as `grid` and `bins` where each cell
        contains the indicies of `z` which contain the values stored
        in the particular bin.

    Revisions
    ---------
    2010-07-11  ccampo  Initial version
    """
    # get extrema values.
    xmin, xmax = x.min(), x.max()
    ymin, ymax = y.min(), y.max()

    # make coordinate arrays.
    xi      = np.arange(xmin, xmax+binsize, binsize)
    yi      = np.arange(ymin, ymax+binsize, binsize)
    xi, yi = np.meshgrid(xi,yi)

    # make the grid.
    grid           = np.zeros(xi.shape, dtype=x.dtype)
    nrow, ncol = grid.shape
    if retbin: bins = np.copy(grid)

    # create list in same shape as grid to store indices
    if retloc:
        wherebin = np.copy(grid)
        wherebin = wherebin.tolist()

    # fill in the grid.
    for row in range(nrow):
        for col in range(ncol):
            xc = xi[row, col]    # x coordinate.
            yc = yi[row, col]    # y coordinate.

            # find the position that xc and yc correspond to.
            posx = np.abs(x - xc)
            posy = np.abs(y - yc)
            ibin = np.logical_and(posx < binsize/2., posy < binsize/2.)
            ind  = np.where(wbin == True)[0]

            # fill the bin.
            bin = z[ibin]
            if retloc: wherebin[row][col] = ind
            if retbin: bins[row, col] = bin.size
            if bin.size != 0:
                binval         = np.median(bin)
                grid[row, col] = binval
            else:
                grid[row, col] = np.nan   # fill empty bins with nans.

    # return the grid
    if retbin:
        if retloc:
            return grid, bins, wherebin
        else:
            return grid, bins
    else:
        if retloc:
            return grid, wherebin
        else:
            return grid

def lh_bins(xlo,xhi,y):
    # Stephanie's xbins.pro slightly pythonic.
    ny=len(y)
    xbins,ybins=np.zeros(ny*2+2),np.zeros(ny*2+2)
    j=0
    for i in range(0,len(ybins)-2,2):
        xbins[i:(i+2)]= xlo[j]
        ybins[(i+1):(i+3)]=y[j]
        j=j+1
    
    xbins[-2:]=xhi[j-1]
    return xbins,ybins

def get_afile(src,search_string):
    import glob
    try:
        files = glob.glob1(src,search_string)
    except IndexError:
        print 'Can''t find',search_string,'in',src
        sys.exit(2)
    return [src+f for f in files]

def closest_match(num,somearray):
    index = 0
    somearray = np.nan_to_num(somearray)
    difference = abs(num-somearray[0])
    for i in range(len(somearray)):
        if difference > abs(num-somearray[i]):
            difference = abs(num-somearray[i])
            index = i
    return index

def files_with_extension(filelist,ext,extra=None):
    files = []
    for file in filelist:
        if file.endswith(ext):
            if extra != None:
                if re.search(extra,file):
                    files.append(file)
            else:
                files.append(file)
    if len(files) == 0: print 'no files found with ',ext,extra
    return files

def uniqify(seq, idfun=None): 
    # order preserving
    if idfun is None:
        def idfun(x): return x
    seen = {}
    result = []
    for item in seq:
        marker = idfun(item)
        # in old Python versions:
        # if seen.has_key(marker)
        # but in new ones:
        if marker in seen: continue
        seen[marker] = 1
        result.append(item)
    return result

def area(p):
    # finds area of ordered closed polygon using greens thm:
    # int (-y*dx + x*dy) = 2 Area
    # abs is because order direction (clock- or counterclock-wise) could
    # give neg value for area.
    return 0.5 * abs(sum(x*dy - y*dx for ((x, y), (dx, dy)) in segments(p)))

def segments(p):
    # this is used in function "area"
    dp = p[1:] - p[0:-1]
    dp = np.append(dp,[[0.,0.]],axis=0)
    return zip(p, dp)

def discrete_colors(Ncolors,colormap='gist_rainbow'):
    colors = []
    cmap = cm.get_cmap(colormap)
    for i in range(Ncolors):
        colors.append( cmap(1.*i/Ncolors) ) # color will now be an RGBA tuple
    return colors
    
def age_strings(ages):
    str_ages= []
    for age in ages:
        str_ages.append(str('%.0e' % age).replace('e+06',' Myr').replace('e+07','0 Myr').replace('e+08','00 Myr').replace('e+09',' Gyr').replace('e+10','0 Gyr'))
    return str_ages

def mass_strings(masses):
    str_masses = []
    for mass in masses:
        str_masses.append(str(mass)+' Msun')
    return str_masses

def bfits2acii(fitsfiles):
    sys.path.append('/Users/Phil/Dropbox/research/python/')
    from Astronomy import get_tab5_trgb_Av_dmod
    for fitsfile in fitsfiles:
        f = pyfits.open(base+fitsfile)
        data = f[1].data 
        cols = f[1].columns
        col_keys = cols.names
        camera = f[0].header['CAMERA']
        ra = data.field('RA')
        dec = data.field('DEC')
        mag1 = data.field('MAG1_'+camera)
        mag1_err =data.field('MAG1_ERR')
        mag1_std =data.field('MAG1_STD')
        mag2 = data.field('MAG2_'+camera)
        mag2_err =data.field('MAG2_ERR')
        mag2_std =data.field('MAG2_STD')
        Target = fitsfile.split('/')[-1].split('_')[1]
        trgb,av,dmod = get_tab5_trgb_Av_dmod(Target)
        o = open(base+fitsfile.replace('fits','dat'),'w')
        o.write('# Target=%s dmod=%s Av=%s Camera=%s \n' % (Target,dmod,av,camera))
        o.write('# RA DEC MAG1 MAG1_ERR MAG1_STD MAG2 MAG2_ERR MAG2_STD \n')
        for i in range(len(ra)):
            o.write('%f %f %f %f %f %f %f %f \n' % (ra[i],dec[i],mag1[i],mag1_err[i],mag1_std[i],mag2[i],mag2_err[i],mag2_std[i]))
        o.close()
        print 'wrote '+Target+': '+base+fitsfile.replace('fits','dat')

class Table(object):
    '''
    use with read_table(filename)
    self.data_array
    self.key_dict
    self.name
    '''
    def __init__(self, data_array, col_keys,name):
        self.data_array = data_array
        self.key_dict = dict(zip(col_keys,range(len(col_keys))))
        self.name = name
        
    def get_row(self,i):
        return self.data_array[i,:]
            
    def get_col(self,key):
        return self.data_array[:,self.key_dict[key]]
    
def get_numeric_data(filename):
    f=open(filename,'r')
    col_keys = f.readline().replace('#','').strip().split()
    f.close()
    print 'get_numeric_data uses np.loadtxt. Use \n from GenUtils impot read_table'
    data = np.loadtxt(filename)
    
    return Table(data,col_keys,filename)

def fortransucks(line):
    # with fortrans fixed column sizes and the need to tell it to make it a number
    # an exponential, if the number is really frickin small, like E-100 it will
    # not write the letter E, just a -100 and that lends for a confusion when
    # entering the number again. I'm just calling that type of number 0.
    tmp2 = []
    tmp = line.split()
    for item in tmp:
        try:
            float(item)
        except ValueError:
            fortransucks_balls = item.split('-')
            if len(fortransucks_balls[-1]) == 3:
                item = ' 0. '
            else:
                weirdness = item.replace('..','.')
                try:
                    float(weirdness)
                    item = weirdness
                except ValueError:
                    print 'Something is fucked, dude.'
        
        tmp2.append(item)
    return tmp2

def read_table(filename,comment_char='#'):
    f=open(filename,'r')
    lines = f.readlines()
    f.close()
    col_keys = lines[0].strip().replace(comment_char,'').split()
    Ncols = len(col_keys)
    Nrows = len([l for l in lines if not l.startswith(comment_char)])
    data = np.ndarray(shape=(Nrows,Ncols), dtype=float)
    row = 0
    for i,line in enumerate(lines):
        if i==0:
            col_keys = lines[0].strip().replace(comment_char,'').split()
            
        if line.startswith(comment_char): continue
        try:
            data[row] = line.strip().split()
        except ValueError:
            data[row] = map(float,fortransucks(line))
        row +=1
    tab = Table(data,col_keys,filename)
    return tab

class BinaryFitsTable(object):
    def __init__(self, data_array, col_keys):
        self.data_array = data_array
        self.key_dict = dict(zip(col_keys,range(len(col_keys))))
        
    def get_row(self,i):
        return self.data_array[i,:]
        
    def get_col(self,key):
        return self.data_array[:,self.key_dict[key]]

def readBinaryFitsTable(filename):
    fits = pyfits.open(filename)
    bdata = fits[1].data
    cols = fits[1].columns
    col_keys = cols.names
    
    Ncols = len(col_keys)
    Nrows = len(bdata.field(0))
    data = np.column_stack(([field for field in bdata.field(field)]))
    # this doesn't work yet... how to fill a nxn data table or is there a better way?
    #BinTab = BinaryFitsTable(data,col_keys)
    return BinTab

def closest_match(num,somearray):
    index = -1
    try:
        difference = abs(num-somearray[0])
    except IndexError:
        print 'closest_match(num,somearray)'
        return 'bust'
    for i in range(len(somearray)):
        if difference > abs(num-somearray[i]):
            difference = abs(num-somearray[i])
            index = i
    return index

def inside(x,y,u,v,verts=False):
    """
    returns the indices of u,v that are within the boundries of x,y.
    """
    if verts != True:
        verts = get_verts(x,y,nbinsx=60,nbinsy=60,smooth=1)
    else:
        verts =  np.column_stack((x,y))
    points = np.column_stack((u,v))
    mask = nxutils.points_inside_poly(points, verts)
    ind = np.nonzero(mask)[0]
    return ind
    
def get_verts(x,y,dx=None,dy=None,nbinsx=None,nbinsy=None,smooth=None):
    if smooth == None:
        smooth = 0
    else:
        smooth = 1
    ymin = y.min()
    ymax = y.max()
    xmin = x.min()
    xmax = x.max()
    
    if dx == None and dy == None:
        dx = (xmax-xmin)/nbinsx
        dy = (ymax-ymin)/nbinsy
        
    if nbinsx == None and nbinsy == None:
        nbinsy = (ymax-ymin)/dy
        nbinsx = (xmax-xmin)/dx
    
    ymid = []
    min_x = [] 
    max_x = []
    for j in range(nbinsy):
        yinner = ymin+j*dy
        youter = ymin+(j+1)*dy
        # counter intuitive because I'm dealing with mags...
        ind = np.nonzero((y > yinner) & (y < youter))[0]  
        if len(ind) > 0:
            if smooth == 1:
                min_x.append(np.average(x[ind])-3.*np.std(x[ind]))
                max_x.append(np.average(x[ind])+3.*np.std(x[ind]))
                ymid.append((yinner+youter)/2.)
            else:
                min_x.append(min(x[ind]))
                max_x.append(max(x[ind]))
                ymid.append((yinner+youter)/2.)
    
    max_x.reverse() 
    ymidr = ymid[:]
    ymidr.reverse()
    
    # close polygon
    max_x.append(min_x[0])
    
    # close polygon
    ymidr.append(ymid[0])
    
    # get verticies of polygon
    xs = np.concatenate((min_x,max_x))
    ys = np.concatenate((ymid,ymidr))
    verts = np.column_stack((xs,ys))
    
    return verts

def savetext(file,x,y,fmtstring,ask='yes'):
    go = 'y'
    if ask == 'yes':
        go = raw_input('Write to '+file+' (y/n)? ')
    if go != 'n':
        file = filecheck(file)
        np.savetxt(file,np.transpose((x,y)),fmt=fmtstring)
        print 'Saved file: '+file
        return 1
    else:
        print 'No file written.'
        return 0

def filecheck(filename):
    """
    Usage   
    vertfile = sys.argv[1]+'.verts'
    vertfile = filecheck(vertfile)
    """
    if os.path.isfile(filename):
        fc1 = raw_input(filename + ' exists. Overwrite? (y/n) ')
        if fc1 == 'n': 
            fc2 = raw_input('Name of new file (enter quits): ')
            if not fc2:
                print 'Fine, have it your way.'
                sys.exit()
        else:
            fc2 = filename
    else:
        fc2 = filename
    return fc2

def bin_up(x,y,nbinsx=None,nbinsy=None,dx=None,dy=None):
    """
    Adapted from contour_plus writen by Oliver Fraser. 
    Commented out a means to bin up by setting dx,dy. I found it more useful
    to choose the number of bins in the x and y direction (nbinsx,nbinsy)
    
    # xrange = np.arange(xmin, xmax, dx)
    # yrange = np.arange(ymin, ymax, dy)
    # nbinsx = xrange.shape[0]
    # nbinsy = yrange.shape[0]
    
    Returns Z, xrange,yrange
    """
    npts = len(x)
    xmin  = float( x.min() ) # cast these for matplotlib's w/o np support
    xmax  = float( x.max() )
    ymin  = float( y.min() )
    ymax  = float( y.max() )
    if nbinsx == None:
        nbinsx = (xmax - xmin) / dx
    if nbinsy == None:
        nbinsy = (ymax - ymin) / dy
    if dx == None:
        dx = (xmax - xmin) / nbinsx
    if dy == None:
        dy = (ymax - ymin) / nbinsy
    
    xrange = np.linspace(xmin, xmax, nbinsx)
    yrange = np.linspace(ymin, ymax, nbinsy)
    Z    = np.zeros((nbinsy, nbinsx))
    xbin = np.zeros(npts)
    ybin = np.zeros(npts)
    # assign each point to a bin
    for i in range(npts):
        xbin[i] = min(int((x[i] - xmin) / dx), nbinsx-1)
        ybin[i] = min(int((y[i] - ymin) / dy), nbinsy-1)
        # and count how many are in each bin
        Z[ ybin[i] ][ xbin[i] ] += 1
    
    return Z,xrange,yrange

class LineBuilder:
    """
    Usage
    line, = ax.plot([0],[0],color='green',lw=3)
    
    linebuilder = LineBuilder(line)
    plt.show()
    
    verts = np.column_stack((linebuilder.xs,linebuilder.ys))
    np.savetxt(filename,np.transpose((linebuilder.xs,linebuilder.ys)),fmt='%7.4f %7.4f')
    """
    def __init__(self, line):
        self.line = line
        self.xs = list(line.get_xdata())
        self.ys = list(line.get_ydata())
        self.cid = line.figure.canvas.mpl_connect('button_press_event', self)
        
    def __call__(self, event):
        print 'Select region with left mouse button. To quit press right mouse button.'
        print 'clicked button=%d, xdata=%f, ydata=%f'%(event.button, event.xdata, event.ydata)
        if event.button == 3: 
            plt.close('all')
            self.xs.pop(0)
            self.ys.pop(0)
            self.xs.append(self.xs[0])
            self.ys.append(self.ys[0])
        else:
            if event.inaxes!=self.line.axes: return
            self.xs.append(event.xdata)
            self.ys.append(event.ydata)
            self.line.set_data(self.xs, self.ys)
            self.line.figure.canvas.draw()  

def galnames(gal):
    ngal = gal.replace('GC','')
    ngal = ngal.replace('SO','')
    ngal = ngal.replace('C','')
    if re.search('HALO',gal): ngal = ngal.split('-HALO')[0]
    if re.search('WIDE',gal): ngal = ngal.split('-WIDE')[0]
    if re.search('DEEP',gal): ngal = ngal.split('-DEEP')[0]
    if re.search('FIELD',gal): ngal = ngal.split('-FIELD')[0]
    if re.search('SGS',gal): ngal = ngal.split('-SGS')[0]
    if re.search('M81K61',gal): ngal = 'KDG61'
    if re.search('M81K64',gal): ngal = 'KDG64'
    if re.search('DDO71',gal): ngal = 'KDG63'
    if re.search('U-5139',gal): ngal = 'HoI'
    if re.search('M81F12D1',gal): ngal = 'KK77'
    if re.search('MESSIER-081-DWARF-A',gal): ngal = 'KDG52'
    if re.search('ANTLIA',gal): ngal = gal.title()
    if re.search('U-04459',gal): ngal = 'DDO53'
    if re.search('U8651',gal): ngal = 'DDO181'
    if re.search('U8760',gal): ngal = 'DDO183'
    if re.search('U9128',gal): ngal = 'DDO187'
    if re.search('M81F6D1',gal): ngal = 'FM1'
    if len(ngal.split('-')) > 2: 
        ngal = ngal.split('-')[0]+ngal.split('-')[1]
    return ngal
    
def tablefile2dict(filename,top_line_comment):
    f = open(filename,'r')
    lines = f.readlines()
    f.close()
    i = 0
    for line in lines:
        if line.startswith(top_line_comment):
            line= line.replace(top_line_comment,'')
            col_keys = line.split()
            i += 1
    datalines_string = lines[i:]
    d = fill_dict(col_keys,datalines_string)
    return d

def fill_dict(keys,datalines_string):
    d = {}
    for key in keys:
        d[key] = []
    for data_line in datalines_string:
        row = data_line.split()
        try:
            for key,i in zip(keys,range(len(keys))):
                d[key].append(is_numeric(row[i]))
        except IndexError:
            print 'missing data!'
            d[key].append('Nan')
    return d

def get_file(ext):
    """
    finds files in current directory and returns a list. If you know there
    is only one file and don't want a list do:
    file = get_file(ext)[0]

    Uncomment the print line if you want to know what was found.
    """
    ls = os.listdir('.')
    file = []
    for l in ls:
        if l.endswith(ext):
            file.append(l)
    # print 'Found ',len(file),' files ending with ',ext, ' in ',os.getcwd()
    if len(file) < 1 :
        print 'No file found with ',ext, 'in ',os.getcwd()

    return file

def get_dirlist(dilimiter=None):
    dirs = []
    list = os.listdir('.')
    for item in list:
        if os.path.isdir(item):
            if dilimiter != None:
                if re.search(dilimiter,item):
                    dirs.append(item)
            else:
                dirs.append(item)
    return dirs

def is_numeric(lit):
    """
    Return value of numeric literal string or ValueError exception
    From http://rosettacode.org/wiki/Determine_if_a_string_is_numeric#Python
    """
    # Empty String
    if len(lit) <= 0:
        return lit    
    # Handle '0'
    if lit == '0': return 0
    # Hex/Binary
    if len(lit) > 1: # sometimes just '-' means no data...
        litneg = lit[1:] if lit[0] == '-' else lit
        if litneg[0] == '0':
            if litneg[1] in 'xX':
                return int(lit,16)
            elif litneg[1] in 'bB':
                return int(lit,2)
            else:           
                try:
                    return int(lit,8)
                except ValueError:
                    pass
    # Int/Float/Complex
    try:
        return int(lit)
    except ValueError:
        pass
    try:
        return float(lit)
    except ValueError:
        pass
    try:
        return complex(lit)
    except ValueError:
        pass
    return lit

def merge(dict1,dict2,key1,key2):
    """
    returns dict3 with the keys from dict1 and dict2 and the rows where
    dict1[key1] == dict2[key2]
    """
    dict3 = {}
    for key,item in dict2.iteritems():
        dict3[key] = []
    for key,item in dict1.iteritems():
        dict3[key] = []
    for i in range(len(dict1[key1])):
        for j in range(len(dict2[key2])):
            if dict1[key1][i] == dict2[key2][j]:
                for key in dict2.keys():
                    dict3[key].append(dict2[key][j])
                for key in dict1.keys():
                    dict3[key].append(dict1[key][i])
    return dict3

def above_below(x,y,m,b):
    above = np.where(np.greater(y, m*x+b),1,0)
    aindx = above.ravel().nonzero()
    
    below = np.where(np.less(y, m*x+b),1,0)
    bindx = below.ravel().nonzero()
    return aindx,bindx
