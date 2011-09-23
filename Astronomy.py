import numpy as np
import sys
from params import TABLE_DIR

#def reddening_vector_patch(R,x,y,dy=0.5,**kwargs):
#    '''
#    example: kwargs={color='red',linewidth=1,length_includes_head=False,head_width=.1}
#    '''
#    dx = np.sqrt(R**2-dy**2)
#    from matplotlib.patches import FancyArrow
#    arr = FancyArrow(x,y,dx,dy,**kwargs)
#    return arr


def reddening_vector_patch(R,mag2_lim = (19,20),color_start=0.4,**kwargs):
    xs = np.linspace(color_start,10,10000)
    ys = xs*R+mag2_lim[0]-xs[0]*R
    plt_inds = np.nonzero(ys<mag2_lim[1])[0]
    
    x = xs[plt_inds][0]
    y = ys[plt_inds][0]
    dx = xs[plt_inds][-1]-x
    dy = ys[plt_inds][-1]-y
    from matplotlib.patches import FancyArrow
    arr = FancyArrow(x,y,dx,dy,**kwargs)
    return arr
    
def reddening_vector(Filter1,Filter2,camera,Rv,Teff=15000):
    '''
    wfc3uvis data from Leo: TABLE_DIR+'extinction_coeff.dat'
    '''
    R = [0]
    if camera == 'wfc3uvis':
        try:
            ec = read_table(TABLE_DIR+'tables/extinction_coeff.dat')
        except NameError:
            from GenUtils import read_table
            ec = read_table(TABLE_DIR+'extinction_coeff.dat')
        Teffs = ec.get_col('Teff')
        Rvs = ec.get_col('Rv')
        A1s = ec.get_col('A('+Filter1+')')
        A2s = ec.get_col('A('+Filter2+')')
        pick = np.nonzero((Rvs == Rv) & (Teffs == Teff))[0]
        A1 = A1s[pick]
        A2 = A2s[pick]
        R = A2/(A1-A2)
    else:
        print 'A_i/A_V information only for wfc3uvis.'
    return R[0]

def hms2dec(h,m,s):
    hours = float(h)+float(m)/60.+float(s)/3600.
    dec = hours*180./12.
    return dec

def dms2deg(d,m,s):
    deg = float(d)+float(m)/60.+float(s)/3600.
    return deg

def deg2hms(deg):
     x = deg*(24./360.)
     h = int(x)
     y = x % h * 60.
     m = int(y)
     s = y % m * 60
     return h,m,s

def deg2dms(deg):
     d = int(deg)
     y = deg % d * 60.
     m = int(y)
     s = y % m * 60
     return d,m,s

def get_cmdlimits(color,mag):
    xmin = color.min()
    xmax = color.max()
    ymin = mag.max()
    ymax = mag.min()
    
    if xmin < -1.5:
        xmin = -1.5
    if xmax > 4.:
        xmax = 4.
    if abs(ymax - ymin) > 8:
        ymax = ymin -8.
    return [xmin,xmax,ymin,ymax]

def get_tab4_completeness_mag(Target,filter):
    #filename needs PropID_Target.
    from LatexUtils import ReadLatexTable
    import re
    tab4 = ReadLatexTable('/Users/Phil/research/Papers/ANGST/Dalcanton09src.tar/tab4.tex')
    j = 0
    for i in range(len(tab4['Target name'])):
        if Target == tab4['Target name'][i]:
            if filter == tab4['Filter'][i]:
                if j>1: print 'more than one match found: ',Target,filter,fifty,camera
                fifty = tab4['50% Completeness mag'][i]
                camera = tab4['Instrument'][i]
                j+=1
    
    #if j>1: print 'more than one match found...'
    try:
        fifty = float(fifty)
    except NameError:
        print Target,filter,'combination not found in tab4.'
        camera = 0
        fifty = 0
    return fifty,camera

def get_tab5_trgb_Av_dmod(Target):
    from LatexUtils import ReadLatexTable
    import re,os
    #tab5 = ReadLatexTable('/astro/users/philrose/python/tab5.tex')
    tab5 = ReadLatexTable(os.environ['PYTHON_CODES']+'/tables/tab5.tex')
    trgb = 0
    Av = 0
    dmod = 0
    for i in range(len(tab5['Target Name'])):
        if Target == tab5['Target Name'][i]:
            trgb = tab5['m_TRGB'][i]
            Av = tab5['A_V'][i]
            dmod = tab5['m-M_0'][i]
    if trgb == 0:
        for i in range(len(tab5['Target Name'])): 
            if re.search(Target[:-1],tab5['Target Name'][i]):
                print 'Table 5 search: Using ',tab5['Target Name'][i],' instead of ',Target
                trgb = tab5['m_TRGB'][i]
                Av = tab5['A_V'][i]
                dmod = tab5['m-M_0'][i]
    if trgb == 0:
        for i in range(len(tab5['Target Name'])): 
            if re.search(Target.split('-')[0],tab5['Target Name'][i]):
                trgb = tab5['m_TRGB'][i]
                Av = tab5['A_V'][i]
                dmod = tab5['m-M_0'][i]
    if trgb == 0:
        print 'Can''t find listing for '+Target+' in ANGST Table 5: /Users/Phil/research/python/tables/tab5.tex.'
        go =raw_input('Would you like to enter the correct Target (y/N)? ')
        if go == 'y' or go == 'Y':
            Target = raw_input('Target: ')
            for i in range(len(tab5['Target Name'])):
                if Target == tab5['Target Name'][i]:
                    trgb = tab5['m_TRGB'][i]
                    Av = tab5['A_V'][i]
                    dmod = tab5['m-M_0'][i]
            print 'passing ',trgb,Av,dmod
        else:
            trgb, Av,dmod = 0.,0.,0.
            print 'passing ',trgb,Av,dmod
    return trgb,Av,dmod
    
def HSTmag2Mag(mag,Target,filter):
    if Target == 'M31':
        trgb = ''
        Av = 0.0
        dmod = 24.47
        print 'No Av yet!!'
    else:
        trgb,Av,dmod = get_tab5_trgb_Av_dmod(Target)
        print 'TRGB, Av, dmod for '+Target+'=',trgb,Av,dmod
        if (Av == 0. and dmod == 0.): 
            A = 0.
            print 'Returning apparent mag for ' + Target+', '+filter
    # m-M = 5 log d - 5 + Av
    # M = m - dmod - A_F*Av
    if (filter == 'F435W'):
        A = 1.30 * Av
    if (filter == 'F475W'):
        A = 1.15 * Av 
    if (filter == 'F555W'):
        A = 1.00 * Av
    if (filter == 'F606W'):
        A = 0.87 * Av
    if (filter == 'F814W'):
        A = 0.57 * Av
    if (filter == 'F110W'):
        A = 0.33926 * Av
    if (filter == 'F160W'):
        A = 0.20443 * Av
    if (filter == 'F275W' or filter == 'F275W1'):
        A = 1.94436 * Av
    if (filter == 'F336W' or filter == 'F336W1'):
        A = 1.65798 * Av
    return mag-dmod-A

def HSTMag2mag(Mag,Target,filter):
    if Target == 'M31':
        trgb = ''
        Av = 0.206
        dmod = 24.47
    else:
        trgb,Av,dmod = get_tab5_trgb_Av_dmod(Target)
        print 'TRGB, Av, dmod for '+Target+'=',trgb,Av,dmod
        if (Av == 0. and dmod == 0.): 
            A = 0.
            print 'Returning apparent mag for ' + Target+', '+filter
    # m-M = 5 log d - 5 + Av
    # M = m - dmod - A_F*Av
    if (filter == 'F435W'):
        A = 1.30 * Av
    if (filter == 'F475W'):
        A = 1.15 * Av 
    if (filter == 'F555W'):
        A = 1.00 * Av
    if (filter == 'F606W'):
        A = 0.87 * Av
    if (filter == 'F814W'):
        A = 0.57 * Av
    if (filter == 'F110W'):
        A = 0.33926 * Av
    if (filter == 'F160W'):
        A = 0.20443 * Av
    if (filter == 'F275W' or filter == 'F275W1'):
        A = 1.94436 * Av
    if (filter == 'F336W' or filter == 'F336W1'):
        A = 1.65798 * Av
    return Mag+dmod+A

def mag2Mag(mag,dmod,a):
    return mag-dmod-a

def Mag2mag(Mag,dmod=24.47,a=0.0):
    return Mag+dmod+a

def Mag2magUVIS(Mag,Filter,Av=0.206):
    ''' 
    This is just a short hand. 
    All variables except Filter and Mag are globals.
    '''
    if Filter == 'F275W':A = Av*af275w
    if Filter == 'F336W':A = Av*af336w
    return Mag+dmod+A

def flux2vegamag(flux,filter,camera):
    if camera == 'acs':
        zeropt_table=TABLE_DIR+'acs_zeropts.dat'
        f_vega = 3.46e-9 # erg cm-2 s-1
    elif camera == 'wfc3uvis':
        zeropt_table=TABLE_DIR+'wfcuvis.tab'
        if filter == 'F275W': f_vega = 3.727e-9 # erg cm-2 s-1 A-1
        if filter == 'F336W': f_vega = 3.245e-9 # erg cm-2 s-1 A-1
    elif camera == 'wfc3ir':
        zeropt_table=TABLE_DIR+'wfcir.tab'
        f_vega = 3.46e-9 # erg cm-2 s-1
    else: 
        print 'camera %s not found. choose acs, wfc3ir, wfc3uvis'
        sys.exit()
    
    lines = open(zeropt_table).readlines()[1:]
    filters = [line.split()[0] for line in lines]
    vegamag0 = [float(line.split()[5]) for line in lines]
    
    try:
        ind = filters.index(filter)
    except ValueError:
        print 'vegamag2flux: filter %s not found.'%filter
        print 'available filters: ',filters
        sys.exit()
        
    return -2.5*np.log10(flux/f_vega)#+vegamag0[ind]


def vegamag2flux(mag,filter,camera):
    '''
    Returns flux in units of  erg cm^{-2} s^{-1} A-1
    F_vega at 5556 A: http://arxiv.org/abs/astro-ph/0403712
    
    wfc3 zeropoints from http://www.stsci.edu/hst/wfc3/phot_zp_lbn# 
    f_vega= $3.46 x10^{-9} erg cm^{-2} s^{-1} A-1 $
    ACS/WFC Zeropoints at -81 C
    Note: For all data obtained after July 4, 2006
    Accessed http://www.stsci.edu/hst/acs/analysis/zeropoints 22/06/11
    if camera == 'acs':
        zeropt_table='/Users/Phil/research/python/tables/acs_zeropts.dat'
    elif camera == 'wfc3uvis':
        zeropt_table='/Users/Phil/research/PHAT/Filters/wfcuvis.tab'
    elif camera == 'wfc3ir':
        zeropt_table='/Users/Phil/research/PHAT/Filters/wfcir.tab'
    '''
    if camera == 'acs':
        zeropt_table=TABLE_DIR+'acs_zeropts.dat'
        f_vega = 3.46e-9 # erg cm-2 s-1
    elif camera == 'wfc3uvis':
        zeropt_table=TABLE_DIR+'wfcuvis.tab'
        if filter == 'F275W': f_vega = 3.727e-9 # erg cm-2 s-1 A-1
        if filter == 'F336W': f_vega = 3.245e-9 # erg cm-2 s-1 A-1
    elif camera == 'wfc3ir':
        zeropt_table=TABLE_DIR+'wfcir.tab'
        f_vega = 3.46e-9 # erg cm-2 s-1
    else: 
        print 'camera %s not found. choose acs, wfc3ir, wfc3uvis'
        sys.exit()
    
    lines = open(zeropt_table).readlines()[1:]
    filters = [line.split()[0] for line in lines]
    vegamag0 = [float(line.split()[5]) for line in lines]
    
    try:
        ind = filters.index(filter)
    except ValueError:
        print 'vegamag2flux: filter %s not found.'%filter
        print 'available filters: ',filters
        sys.exit()
    return 10**(-0.4*(mag))*f_vega#-vegamag0[ind]))*f_vega
'''
old!
def vegamag2flux(mag,filter):

    Zeropoints from http://www.stsci.edu/hst/wfc3/phot_zp_lbn
    F_vega at 5556 A: http://arxiv.org/abs/astro-ph/0403712
    returns flux in erg cm-2 s-1

    if filter == 'F275W' or filter == 'F275W1' : zero = 22.65
    elif filter == 'F336W' or filter == 'F336W1' : zero = 23.46
    else: 
        print filter,'not found! update astronomy.py \n http://www.stsci.edu/hst/wfc3/phot_zp_lbn'
        return 0
    f_vega = 3.46e-9 # erg cm-2 s-1
    return 10**(-0.4*(mag-zero))*f_vega
'''
def plot_cmd(fitsfile,M2,Abs=False,oplot=False,color='black'):
    from MatchUtils import readbintab
    import matplotlib.pyplot as plt
    fits = readbintab(fitsfile)
    filter1 = fits['names'][2]
    filter2 = fits['names'][3]
    Target = fits['names'][1]
    title = fits['names'][0]+' '+fits['names'][1]
    if Abs == True:
        mag1 = HSTmag2Mag(fits['Mag1'],Target,filter1)
        mag2 = HSTmag2Mag(fits['Mag2'],Target,filter2)
    else:
        mag1 = fits['Mag1']
        mag2 = fits['Mag2']
    color = mag1-mag2
    
    if oplot == False: fig = plt.figure()
    if M2 == 'I':
        magy = mag2
        plt.ylabel(r'$\rm{%s}$'%filter2)
    if M2 == 'V':
        magy = mag1
        plt.ylabel(r'$\rm{%s}$'%filter1)
    plt.xlabel(r'$\rm{%s-%s}$'%(filter1,filter2))
    plt.title('$\rm{%s}$'%title.replace(' ','\ '))
    plt.scatter(color,magy,s=1,marker='o',color=color,edgecolor=None)
    plt.axis(get_cmdlimits(color,mag2))
    return mag1,mag2

def rotate(x,y,r,theta):
    x_rot = x*np.cos(theta)-y*np.sin(theta)
    y_rot = x*np.sin(theta)+y*np.cos(theta)
    #x_rot = r*np.cos(theta)+x
    #y_rot = r*np.sin(theta)+y
    return x_rot,y_rot

def rotate_box(xc,yc,w,h,angle):
    # Up Down Left Right
    theta = angle*np.pi/180.
    print theta,angle
    xl = -w/2.
    xr = w/2.
    yu = h/2
    yd = -h/2
    print yu,yd,yc
    r = np.sqrt((w/2.)**2+(h/2)**2)
    xul,yul = rotate(xl,yu,r,theta)
    xur,yur = rotate(xr,yu,r,theta)
    xdr,ydr = rotate(xr,yd,r,theta)
    xdl,ydl = rotate(xl,yd,r,theta)
    
    xs = [xul+xc,xur+xc,xc+xdr,xc+xdl,xul+xc]
    ys = [yul+yc,yur+yc,yc+ydr,yc+ydl,yul+yc]
    return np.column_stack((xs,ys))

def color_mag_errs(fitstable,yfilt='I',binsize=0.1,camera='UVIS'):
    '''    
    p2 = plt.errorbar(errcol+3, errmag, xerr=errcolerr, yerr=errmagerr,ecolor='white',lw=3,capsize=0,fmt=None)
    p2 = plt.errorbar(errcol+3, errmag, xerr=errcolerr, yerr=errmagerr,ecolor='black',lw=2,capsize=0,fmt=None)
    '''
    import pyfits
    import numpy as np
    fits = pyfits.open(fitstable)
    data = fits[1].data
    try:
        mag1 = data.field('MAG1_'+fits[0].header['CAMERA'])
        mag2 = data.field('MAG2_'+fits[0].header['CAMERA'])
    except KeyError:
        mag1 = data.field('MAG1_'+camera)
        mag2 = data.field('MAG2_'+camera)
    
    mag1err = data.field('MAG1_ERR')
    mag2err = data.field('MAG2_ERR')
    
    x = mag1-mag2 # color
    # v or i for y-axis?
    if yfilt == 'I':
        y = mag2
        yerr = mag2err
    else:
        y = mag1
        yerr = mag1err
    
    nbins = (y.max()-y.min())/binsize
    errclr = -1.
    errmag = np.arange(int(nbins/5)-1)*1.
    errcol = np.arange(int(nbins/5)-1)*1.
    errmagerr = np.arange(int(nbins/5)-1)*1.
    errcolerr = np.arange(int(nbins/5)-1)*1.
    
    for q in range(len(errmag)-1):
        test = y.min()+5.*(q+2)*binsize+2.5*binsize
        test2 = np.nonzero( (mag2 > test-2.5*binsize) & (mag2 <= test+2.5*binsize) & (mag1-mag2>-0.5) & (mag1-mag2 < 2.5))[0]
        if len(test2) < 5 : continue
        errmag[q] = y.min()+5.*(q+2)*binsize+2.5*binsize
        errcol[q] = errclr
        errmagerr[q] = np.sqrt(np.mean(yerr[np.nonzero((y > errmag[q]-2.5*binsize) & (y < errmag[q]+2.5*binsize))[0]]**2))
        errcolerr[q] = np.sqrt(np.mean(
                            mag1err[np.nonzero((y>errmag[q]-2.5*binsize) & (y < errmag[q]+2.5*binsize) & (x > -0.5) & (x<2.5))[0]]**2 +
                            mag2err[np.nonzero((y>errmag[q]-2.5*binsize) & (y < errmag[q]+2.5*binsize) & (x > -0.5) & (x<2.5))[0]]**2))
    return errcol[:-1],errmag[:-1],errcolerr[:-1],errmagerr[:-1]

def mean_proj_radius(ra,dec,center_ra=10.68463,center_dec=41.26693,pos_angle=20.,eccentricity=0.804):
    '''
    eqn 3.2 pg 45 of Andrew West PhDTh 
    pos_angle is in degrees
    
    19.4.11: Changed default position angle to 20. Plotting showed it was off,
    angle/2. guess worked well. Still may be off...
    '''
    theta = pos_angle*np.pi/180.
    dra = ra-center_ra
    ddec = dec-center_dec
    
    term1 = dra*np.cos(theta)+ddec*np.sin(theta)
    term2 = dra*np.sin(theta)-ddec*np.cos(theta)
    term3 = 1.-eccentricity**2
    return np.sqrt(term1**2+term2**2/term3)

def dmod2kpc(dmod):
    return 10**(0.2*dmod+1.)/1000.

def arcsec2radian(arcsec):
    return arcsec/206264.806247

def arcsec2kpc(arcsec,dmod=24.47):
    print 'using dmod =',dmod
    r = dmod2kpc(dmod)
    theta = arcsec2radian(arcsec)
    return r*theta


def solar_consts():
    lsun = 3.845e33
    mbol = -26.83
    Mbol = 4.74
    BC = -0.08
    
    
