import sys,os,re
import numpy as np

# For polygon
import Polygon
import matplotlib.nxutils as nxutils

# For parsing fits files and reg files
import pyfits
import pyregion
from pyregion.mpl_helper import properties_func_default
import pywcsgrid2
from kapteyn import wcs

from PadovaTracksUtils import *
from GenUtils import get_afile,is_numeric,discrete_colors,area,bin_up,tablefile2dict
from Astronomy import Mag2mag,color_mag_errs,vegamag2flux

m31 = '/Users/Phil/research/PHAT/m31_fullfield.fits'

# MAC
data_src = '/Users/Phil/research/PHAT/Data/'
base_fits = '/Users/Phil/research/PHAT/Data/12058_M31-B01-F10-UVIS_F336W_drz.chip1.fits'

# astro network
data_src = '/astro/net/angst2/philrose/PHAT/FullBricks/'
base_fits = '/astro/net/angst2/philrose/PHAT/Fields/drzs/12058_M31-B01-F10-UVIS_F336W_drz.chip1.fits'

print 'fits are centered with',base_fits
f = pyfits.open(base_fits)
f_shape = (f[0].header["NAXIS1"], f[0].header["NAXIS2"])
f_data = f[0].data
f_proj = wcs.Projection(f[0].header)
f_header = f[0].header
f.close()
def read_basefits():
    return f


def radec_to_reg(ra,dec,outfile='radec.reg',shape='circle',size=0.5,header='default'):
    if header == 'default':
        head = 'global color=green dashlist=8 3 width=1 font=\"helvetica 10 normal\" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1 \n fk5'
    else:
        head = header
    
    f = open(outfile,'w')
    f.write('%s'%head)
    [f.write('%s( %f, %f,%f")\n'%(shape,ra[i],dec[i],size)) for i in range(len(ra))]
    print 'radec_to_reg wrote',outfile
    f.close()

def reg_to_array(reg_file,commentchar='#',tag=None):
    import re
    ra,dec,tags = np.array([]),np.array([]),np.array([])
    for r in open(reg_file,'r').readlines():
        r = r.strip()
        if r.startswith(commentchar): continue
        if r.startswith('circle'):
            ra = np.append(ra,float(r.split(',')[0].split('(')[-1]))
            dec = np.append(dec,float(r.split(',')[1]))
            tags = np.append(tag,str(r.split('{')[-1].replace('}','')))
        if r.startswith('polygon'):
            if tag==None:
                data = map(float,r.split('(')[-1].split(')')[0].split(','))
                ra = np.append(ra,data[::2])
                dec = np.append(dec,data[1::2])
                tags = np.append(tag,str(r.split('{')[-1].replace('}','')))
            else:
                if re.search(tag,r):
                    data = map(float,r.split('(')[-1].split(')')[0].split(','))
                    ra = data[::2]
                    dec = data[1::2]
    return ra,dec,tags

def read_6_filters(file='/Users/Phil/research/PHAT/Data/sixfilt-brick01-v4st.fits'):
    fitsobj=pyfits.open(file)
    data = fitsobj[1].data
    Data= {}
    keys = ['acs_ra','ir_ra','uv_ra','acs_dec','ir_dec','uv_dec','acs_nmatch','ir_nmatch','uv_nmatch','acs_mag1','ir_mag1','uv_mag1','acs_mag2','ir_mag2','uv_mag2','acs_mag1err','ir_mag1err','uv_mag1err','acs_mag2err','ir_mag2err','uv_mag2err','acs_mag1samvar','ir_mag1samvar','uv_mag1samvar','acs_mag2samvar','ir_mag2samvar','uv_mag2samvar','ra','dec']
    for key in keys:
        Data[key]=data.field(key)
    fitsobj.close()
    return Data

def closest_match(num,somearray):
    index = -1
    somearray = np.nan_to_num(somearray)
    difference = abs(num-somearray[0])
    for i in range(len(somearray)):
        if difference > abs(num-somearray[i]):
            difference = abs(num-somearray[i])
            index = i
    return index

def read_fullbrick(table):
    t = pyfits.open(table)
    data = t[1].data
    return data

def area_within_reg(reg_name,base_fits=m31):
    '''
    returns deg^s 
    deg^2 to arcsec^2: areas*3600.*3600.
    deg^2 to arcmin^2: areas*3600.
    '''
    # finds area of concentric reg shapes
    # use GenUtils.area for the area of a polygon.
    tot_areas,areas = [],[]
    proj = wcs.Projection(f_header)
    r = pyregion.open(reg_name).as_imagecoord(header=f_header)
    patch_list, artist_list = r.get_mpl_patches_texts()
    for p in patch_list[::-1]:
        verts_xy = p.get_xy()
        pixels = (verts_xy[:,0],verts_xy[:,1])
        ra_v,dec_v = proj.toworld(pixels)
        verts = np.column_stack((ra_v,dec_v))
        tot_areas.append(area(verts))
        
    for i in range(len(tot_areas)):
        if i == len(tot_areas)-1:
            areas.append(tot_areas[i])
        else:
            areas.append(tot_areas[i]-tot_areas[i+1])

    # deg^2 to arcsec^2: areas*3600.*3600.
    # deg^2 to arcmin^2: areas*3600.
    return np.array(areas)

def read_all_fake(allfakes,mag1_cut=999,mag2_cut=999):
    '''
    col_num: 0 ... 19
    ? chip X       Y      counts magin
    0 1  197.04   11.05   3021.4 20.751   2947.7 20.751   1860.4 20.751   1833.6 20.751   8742.9 20.606   8898.6 20.606   6888.9 20.606   6899.2 20.606 

    20 21  22     23       24     25      26    27     28    29    30   
    ? chip X       Y      Chi     SNR   sharp  round majax crowd type 
    0 1  197.03   11.05   0.96    80.5  0.008  0.222 110 0.000 1    

    31         32    33   34     35      36     37      38    39     40    41     42     43     44     45    46  
    counts     bg   mag           mag_err chi    snr     sharp round  crowd fwhm   ell    psfa   psfb   psfc   err
    9481.6    11.74 20.778 99.999 0.029   0.78    37.1  0.004  0.025 0.000  1.000  0.000  1.000  1.000  0.000  3   

    47        48     49      50     51     52      53    54
    32043.0   134.65 20.578 99.999 0.015   1.09    71.4  0.011  0.384 0.000  1.000  0.000  1.000  1.000  0.000  3

    2947.0     3.65 20.778 99.999 0.029   0.78    37.1  0.004  0.025 0.000  1.000  0.000  1.000  1.000  0.000  3       
    0.0     0.00 99.999 99.999 0.000   0.00     0.0  0.000  0.000 0.000  0.000  0.000  0.000  0.000  0.000  9       
    0.0     0.00 99.999 99.999 0.000   0.00     0.0  0.000  0.000 0.000  0.000  0.000  0.000  0.000  0.000  9       
    0.0     0.00 99.999 99.999 0.000   0.00     0.0  0.000  0.000 0.000  0.000  0.000  0.000  0.000  0.000  9    
    8972.0    37.70 20.578 99.999 0.015   1.09    71.4  0.011  0.384 0.000  1.000  0.000  1.000  1.000  0.000  3       
    0.0     0.00 99.999 99.999 0.000   0.00     0.0  0.000  0.000 0.000  0.000  0.000  0.000  0.000  0.000  9       
    0.0     0.00 99.999 99.999 0.000   0.00     0.0  0.000  0.000 0.000  0.000  0.000  0.000  0.000  0.000  9       
    0.0     0.00 99.999 99.999 0.000   0.00     0.0  0.000  0.000 0.000  0.000  0.000  0.000  0.000  0.000  9
    '''
    if type(allfakes) is not list: 
        x,y,mag1in,mag2in,mag1out,snr1,sharp1,mag2out,snr2,sharp2 = np.loadtxt(allfakes,unpack=True,usecols=(2,3,5,15,33,37,38,49,53,54))
    else:
        for j,allfake in enumerate(allfakes):
            print 'reading',allfake.split('/')[-1]
            # see huge comment or *.phot.columns in MC directory for column headings.
            data = np.loadtxt(allfake,unpack=True,usecols=(2,3,5,15,33,37,38,49,53,54))
            if j == 0: 
                x,y,mag1in,mag2in,mag1out,snr1,sharp1,mag2out,snr2,sharp2 = data
                continue
            else: 
                x = np.append(x,data[0])
                y = np.append(y,data[1])
                mag1in = np.append(mag1in,data[2])
                mag2in = np.append(mag2in,data[3])
                mag1out = np.append(mag1out,data[4])
                snr1 = np.append(snr1,data[5])
                sharp1 = np.append(sharp1,data[6])
                mag2out = np.append(mag2out,data[7])
                snr2 = np.append(snr2,data[8])
                sharp2 = np.append(sharp2,data[9])
    # do quality cut
    ind1 = quality_cut(sharp1,snr1,mag1in,mag_cut=mag1_cut)
    ind2 = quality_cut(sharp2,snr2,mag2in,mag_cut=mag2_cut)       
    return x,y,mag1in,mag2in,mag1out,mag2out,ind1,ind2

def read_phot(photfile):
    # see *.phot.columns in MC directory for column headings.
    x,y,mag1,snr1,sharp1,mag2,snr2,sharp2 = np.loadtxt(photfile,unpack=True,usecols=(2,3,13,17,18,29,33,34))
    # indices of each star that pass quality cuts (in each filter)
    ind1 = quality_cut(sharp1,snr1,mag1)
    ind2 = quality_cut(sharp2,snr2,mag2)
    ind12 = list(set(ind1) & set(ind2))
    
    # single detections are defined by stars in filt A that pass quality cut
    # and their filt B couterpart is not resolved (i.e. mag B = 99)
    single1 = np.nonzero(mag2[ind1] > 90)[0]
    single2 = np.nonzero(mag1[ind2] > 90)[0]
    
    return x,y,mag1,mag2,ind1,ind2,ind12

def XY2radec(pyfits_obj,x,y):
    # returns ra,dec based on (an already opened) fits image.
    proj = wcs.Projection(pyfits_obj[0].header)
    pixels = (x,y)
    return proj.toworld(pixels)
    
def points_inside_reg(pyfits_obj,reg_name,ra_points,dec_points,plot=False):
    '''
    Takes an open fits file (pyfits_obj), a reg file (reg_name), ra and dec
    array, and returns the indices of the ra dec array inside each shape in
    the region file. Will also plot...
    inds - cumulative
    '''
    if plot == True:
        ax=pywcsgrid2.subplot(111, header=pyfits_obj[0].header)
        ax.grid()
        ax.imshow(pyfits_obj[0].data, cmap=cm.gray, origin="lower")
        
    # load projection, reg shapes, and get them in python
    proj = wcs.Projection(pyfits_obj[0].header)
    r = pyregion.open(reg_name).as_imagecoord(header=pyfits_obj[0].header)
    patch_list, artist_list = r.get_mpl_patches_texts()
    
    radec = np.column_stack((ra_points,dec_points))
    # get masks of ra,dec in each shape
    masks = []
    for p in patch_list[::-1]:
        if plot == True: ax.add_patch(p)
        verts_xy = p.get_xy()
        pixels = (verts_xy[:,0],verts_xy[:,1])
        ra_v,dec_v = proj.toworld(pixels)
        verts = np.column_stack((ra_v,dec_v))
        masks.append( nxutils.points_inside_poly(radec, verts))
        
    # subtract masks from outside to inside to not over count.
    inds = []
    for i in range(len(masks)-1):
        isolated_mask = (masks[i]-masks[i+1])
        if len(isolated_mask)>0:
            within = np.nonzero(isolated_mask)[0]
            inds.append(within)
            if plot == True: 
                ax["fk5"].plot(ra_points[within],dec_points[within],'.',alpha=0.2,label='%i'%i,color=cols[i])
                for t in artist_list:
                    ax.add_artist(t)
                ax.add_compass(loc=1)
                ax.legend()
    return inds

def fake_fraction_recovered(Fra1,Fdec1,Fra2,Fdec2,Mag1in_all,Mag2in_all,fra1_qc,fdec1_qc,fra2_qc,fdec2_qc,Mag1in,Mag2in,reg_name,base_fits=m31,labs=None,continuous=True,fraction=True,output=False):
    '''
    Takes a list of the *all.fake files, their associated drizzled fits image
    (for x,y to ra, dec conversions) (be sure the lists are in the same order)
    a reg file to check which stars are within - first written with polynomial
    reg shapes as a function of radius. The regfile and the fake stars are 
    alligned with the base_fits file, which is set by default as m31 full field,
    a global variable.

    Returns the fraction of stars recovered per reg shape in each filter of 
    all.fake.
    See points_inside_reg for what gets defined as star recovered.

    Before I wrote remove_overlaps:
    for i in range(len(allfakes)):
        d = pyfits.open(drzs[i])
        print 'reading',allfakes[i].split('/')[-1]
        # reading in all fakes into big ol np arrays. 
        # NB ind1f,ind2f are the indices that survive the quality cut.
        xf,yf,mag1in,mag2in,mag1out,mag2out,ind1f,ind2f = read_all_fake(allfakes[i])
        fra,fdec = XY2radec(d,xf,yf)
        if i == 0:
            Mag1in_all = mag1in
            Mag2in_all = mag2in
            Mag1in = mag1in[ind1f]
            Mag2in = mag2in[ind2f]
            Mag1out = mag1out[ind1f]
            Mag2out = mag2out[ind2f]
            fra1_qc = fra[ind1f]
            fdec1_qc = fdec[ind1f]
            fra2_qc = fra[ind2f]
            fdec2_qc = fdec[ind2f]
            Fra = fra[:]
            Fdec = fdec[:]
            continue
        else:
            Mag1in = np.append(Mag1in,mag1in[ind1f])
            Mag2in = np.append(Mag2in,mag2in[ind2f])
            Mag1out = np.append(Mag1out,mag1out[ind1f])
            Mag2out = np.append(Mag2out,mag2out[ind2f])
            # asts that survive the quality cuts
            fra1_qc = np.append(fra1_qc,fra[ind1f])
            fdec1_qc = np.append(fdec1_qc,fdec[ind1f])
            fra2_qc = np.append(fra2_qc,fra[ind2f])
            fdec2_qc = np.append(fdec2_qc,fdec[ind2f])
            # all asts
            Fra = np.append(Fra,fra)
            Fdec = np.append(Fdec,fdec)
            Mag1in_all = np.append(Mag1in_all,mag1in)
            Mag2in_all = np.append(Mag2in_all,mag2in)
    '''

    if output == True:
        out1 = open('M31_B01_Completeness_F275W.dat','w')
        out2 = open('M31_B01_Completeness_F336W.dat','w')
        out1.write('# radial bin, mag bin, number %s \n')
        out2.write('# radial bin, mag bin, number %s \n')

    # asts that survive the quality cut
    inds1_qc = points_inside_reg(f,reg_name,fra1_qc,fdec1_qc)
    inds2_qc = points_inside_reg(f,reg_name,fra2_qc,fdec2_qc)
    # all asts
    inds1 = points_inside_reg(f,reg_name,Fra1,Fdec1)
    inds2 = points_inside_reg(f,reg_name,Fra2,Fdec2)    

    #cols=discrete_colors(len(inds1))

    # plots
    fig = plt.figure(3, figsize=(6,6))
    left, width = 0.15, 0.35
    bottom, height = 0.1, 0.8
    left2 = left+width + 0.04
    
    axis1 = [left, bottom, width, height]
    axis2 = [left2, bottom, width, height]
    
    ax1 = plt.axes(axis1)
    ax2 = plt.axes(axis2)
    
    ax1.set_xlabel(r'$F275W$',fontsize = ftsize)
    ax2.set_xlabel(r'$F336W$',fontsize = ftsize)
    
    ax1.set_xlim( (26, 16) ) # set all axes limits here
    ax2.set_xlim( ax1.get_xlim() )
    
    if fraction==True:
        ax1.set_ylim( (0.98, 1.001) )
        ax1.set_ylabel(r'$\rm{Fraction\ of\ Artificial\ Stars\ Recovered}$',fontsize = ftsize)
        ax2.set_ylim( ax1.get_ylim() )
    else:
        ax1.set_ylabel(r'$\rm{\#\ of\ Artificial\ Stars\ Recovered}$',fontsize = ftsize)
    
    nullfmt   = NullFormatter() # no labels
    ax2.yaxis.set_major_formatter(nullfmt)
    #nbins = np.array([16,18,20,22,24,26])
    nbins = np.array([16,17,18,19,20,21,22,23,24,25,26])
    offset = 0.1
    mag1_cuts,mag2_cuts = [], []
    for i in range(len(inds1_qc)):
        labsA1 = 'Full bin %i F275W: %i/%i = %.2f'%(i,len(inds1_qc[i]),len(inds1[i]),float(len(inds1_qc[i]))/float(len(inds1[i])))
        labsA2 = 'Full bin %i F336W: %i/%i = %.2f'%(i,len(inds2_qc[i]),len(inds2[i]),float(len(inds2_qc[i]))/float(len(inds2[i])))
        if labs == None:
            labs1 = labsA1
            labs2 = labsA2
        else:
            labs1 = labs[i]
            labs2 = labs[i]
        if continuous == True:
            hist01, bins = continuous_hist(Mag1in_all[inds1[i]],nbins,offset) # all asts
            hist1, bins = continuous_hist(Mag1in[inds1_qc[i]],nbins,offset) # quality cuts
            hist02, bins = continuous_hist(Mag2in_all[inds2[i]],nbins,offset)
            hist2, bins = continuous_hist(Mag2in[inds2_qc[i]],nbins,offset)
        else:
            hist01,bins = np.histogram(Mag1in_all[inds1[i]],bins=nbins)
            hist1,bins = np.histogram(Mag1in[inds1_qc[i]],bins=nbins)
            hist02,bins = np.histogram(Mag2in_all[inds2[i]],bins=nbins)
            hist2,bins = np.histogram(Mag2in[inds2_qc[i]],bins=nbins)
            bins = bins[1:]
        
        f_recovered1 = (hist1*1.)/(hist01*1.) # normalized (and made into float)
        f_recovered2 = (hist2*1.)/(hist02*1.)        
        comp_limit1 = closest_match(fr,f_recovered1)
        comp_limit2 = closest_match(fr,f_recovered2)
        if comp_limit1 == -1:
            print 'mag 1 comp not found?',f_recovered1
        if comp_limit2 == -1:
            print 'mag 2 comp not found?',f_recovered2

        if output==True:
            for h,b in zip(f_recovered1,bins):
                #out1.write('%i %.1f %.3e \n'%(i,b,h/areas[i]))
                out1.write('%i %.1f %.3f \n'%(i,b,h))
                
        if output==True:
            for h,b in zip(f_recovered2,bins):
                #out2.write('%i %.1f %.3e \n'%(i,b,h/areas[i]))
                out2.write('%i %.1f %.3f \n'%(i,b,h))
    
        print labsA1,'%.2f complete:'%fr,bins[comp_limit1]
        print labsA2,'%.2f complete:'%fr,bins[comp_limit2]
        mag1_cuts.append(bins[comp_limit1])
        mag2_cuts.append(bins[comp_limit2])
        
        pos_err1= (np.sqrt(hist1*1.*float(len(inds1_qc[i]))/float(len(inds1[i]))))/100.
        if fraction == True:
            ax1.plot(bins,f_recovered1,drawstyle='steps',linewidth=2,color='white')
            ax1.plot(bins,f_recovered1,drawstyle='steps',linewidth=1.5,color=cols[i],label=labs1)
            ax2.plot(bins,f_recovered2,drawstyle='steps',zorder=1,linewidth=2,color='white')
            ax2.plot(bins,f_recovered2,drawstyle='steps',color=cols[i],zorder=1.5,linewidth=2,label=labs2)
        else:
            ax1.semilogy(bins[1:],hist1*1.,drawstyle='steps',color=cols[i],zorder=1)
            ax1.semilogy(bins[1:],hist01*1.,'--',drawstyle='steps',color=cols[i],zorder=1)
            ax2.semilogy(bins[1:],hist2*1.,drawstyle='steps',color=cols[i],zorder=1)
            ax2.semilogy(bins[1:],hist02*1.,'--',drawstyle='steps',color=cols[i],zorder=1)

    if output==True:
        out1.close()
        out2.close()
        print 'wrote M31_B01_Completeness_F275W.dat'
        print 'wrote M31_B01_Completeness_F336W.dat'
    ax1.legend(loc=3)
    #ax2.legend(loc=3)
    plt.savefig('completeness.png')
    print 'wrote completeness.png'
    plt.close()
    print 'found mag1_cut=%f and mag2_cut=%f'%(np.min(mag1_cuts),np.min(mag2_cuts))
    return np.min(mag1_cuts),np.min(mag2_cuts)

def read_table_mixed_type(file,delimiter=','):
    # writen to read a csv file with mixed float/str
    # add funcionality for comment lines?
    lines = open(file,'r').readlines()
    col_keys = lines[0].split(delimiter)
    Ncols = len(col_keys)
    Nrows = len(lines)-1
    data = {}
    for key in col_keys:
        data[key] = []
    for line in lines[1::]:
        item = line.strip().split(delimiter)
        for j,key in enumerate(col_keys):
            data[key].append(is_numeric(item[j]))
    return data

def exp_tab2reg(file):
    # I already did a grep of UVIS and then F336W, you know, fuck re.
    out = open('12058_exposure_coverage_UVIS_F336W.reg','w')
    out.write('# Region file format: DS9 version 4.1 \n')
    out.write('global color=green dashlist=8 3 width=1 font="helvetica 10 normal" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1\n')
    out.write('fk5\n')
    lines = open(file,'r').readlines()
    for line in lines:
        item = line.strip().split(',')
        tag = '_'.join((item[14],item[8],str(int(float(item[10])))))
        x,poly1,poly2 = item[-1].strip().split('Polygon J2000 ')
        out.write('polygon(%s) # tag=%s-1\n' % (poly1.strip().replace(' ',','),tag))
        out.write('polygon(%s) # tag=%s-2\n' % (poly2.strip().replace(' ',','),tag))
    out.close()

def fixed_color(shape, saved_attrs):
    attr_list, attr_dict = saved_attrs
    attr_dict["color"] = "blue"
    kwargs = properties_func_default(shape, (attr_list, attr_dict))
    return kwargs

def mpl_patch_XY2radec(patch,proj=f_proj):
    # returns ra,dec
    verts = patch.get_xy()
    pixels = (verts[:,0],verts[:,1])
    return proj.toworld(pixels)

def ndarray2tuple(ndarray):
    return tuple([tuple(l) for l in ndarray])

def poly_over_under(ra,dec,patch_left,patch_right):
    '''
    Short cut to Polygon.Polygon
    (p and q are polygons)
    p & q: intersection: a polygon with the area that is covered by both p and q
    p | q: union: a polygon containing the area that is covered by p or q or both
    p - q: difference: a polygon with the area of p that is not covered by q
    p + q: sum: same as union
    p ^ q: xor: a polygon with the area that is covered by exactly one of p and q

    len(p):
    p[i]:

    number of contours
    contour with index i, the same as p.contour(i), slicing is not yet supported
    '''
    # returns indices of ra, dec that are not in the overlap region 
    # of the two patches.
    # use ra,dec.
    lr,ld = mpl_patch_XY2radec(patch_left)
    rr,rd = mpl_patch_XY2radec(patch_right)
    # Polygon needs tuples...
    left = Polygon.Polygon(ndarray2tuple(np.column_stack((lr,ld))))
    right = Polygon.Polygon(ndarray2tuple(np.column_stack((rr,rd))))
    # p & q: intersection: a polygon with the area that is covered by both p and q
    overlap = right & left
    # vertices of overlapped region
    if len(overlap[0]) > 0:
        verts = np.transpose(np.array(overlap[0]))
        radec = np.column_stack((ra,dec))
        mask = nxutils.points_inside_poly(radec, np.array(overlap[0]))
    # sloppy? flip T/F: subtract 1 from bool array and then get rid of neg sign
        not_overlapped = np.nonzero(abs(mask-1))[0]
    else: not_overlapped = []
    return not_overlapped

def in_chips(ra,dec,patch_chip1,patch_chip2,proj=f_proj):
    # this is meant to exlude stars that fall in the chip gap
    # TO DO: make each region edge 2 arcsec smaller...
    within = []
    for patch in [patch_chip1,patch_chip2]:
        verts_xy = patch.get_xy()
        pixels = (verts_xy[:,0],verts_xy[:,1])
        ra_v,dec_v = proj.toworld(pixels)
        verts = np.column_stack((ra_v,dec_v))
        radec = np.column_stack((ra,dec))
        mask = nxutils.points_inside_poly(radec, verts)
        within.extend(np.nonzero(mask)[0])
    return within
    
def match_2inds(inds1,inds2):
    # takes two sets of indices and returns a new matching set
    return list(set(inds1) & set(inds2))

def match_3inds(inds1,inds2,inds3):
    return list(set(inds1) & set(inds2) & set(inds3))
    
def remove_field_overlaps(photfiles,drzs,footprints,base_fits=m31,plot=False,fake=False,matched_only=False,qc=True):
    '''
    This only works for field 3,4,5,9,10,11,15,16,17. It can be adapted for 
    all fields though, 3 is now the top left field, 9 and 15 are the left most
    fields on the next rows.

    footprints is a region file from exposure_coverage.csv
    pre-processing:
    1. grep UVIS exposure_coverage.csv > poop
    2. grep F336W poop > poop1
    3. grep 550 poop1 > poop2
    4. Then some hard coding in the file:
      I decided to update the exposure footprints file since full bricks are 
      on their way...
      the problem is that chip 1 vs chip 2 flip their orientation on fields 1,2,3
      from those of 4,5,6 such that chip 1 of field 3 overlaps with chip 1 of 
      field 4 while chip 2 of field 4 and 5 overlap with chip 1 of 5 and 6.
      I swapped the chip values of 1,2,3,7,8,9,13,14,15 so it's consistent with
      4,5,6,10,11,12,16,17,18. 
    5. Then I cut out all fields BUT 3,4,5,9,10,11,15,16,17

    Hacked adaptability for all.fake files. If .fake files are coming in, just
    say they are photfiles and set the fake kwarg to True.

    For all.fake: need to run twice (ug). Once with qc=True and once with qc=False.

    adding more fields? Check the if statements and the [i-*]s.
    '''
    clrs = discrete_colors(len(photfiles))
    if plot == True:
        ax=pywcsgrid2.subplot(111, header=f_header)
        ax.grid()
        ax.imshow(f_data, cmap=cm.gray, origin="lower")
    r = pyregion.open(footprints).as_imagecoord(header=f_header)
    # a is the left chip of the field (see # 4 in comment)
    r1a = pyregion.ShapeList([rr for rr in r if re.search('550-1',rr.attr[1].get("tag"))])
    # b is the right chip of the field (see # 4 in comment)
    r1b = pyregion.ShapeList([rr for rr in r if re.search('550-2',rr.attr[1].get("tag"))])

    patch_list1a, artist_list1a = r1a.get_mpl_patches_texts()
    patch_list1b, artist_list1b = r1b.get_mpl_patches_texts(fixed_color)
    if plot == True:
        for p in patch_list1a + patch_list1b:
            ax.add_patch(p)
    for i in range(len(photfiles)):
        inds1,inds2 = '',''
        print 'reading',photfiles[i]
        if fake==True: 
            if qc == False:
                x,y,mag1,mag2,xxx,yyy,qc_ind1,qc_ind2 = read_all_fake(photfiles[i])
            else:
                x,y,mag1,mag2,xxx,yyy,qc_ind1,qc_ind2 = read_all_fake(photfiles[i],mag1_cut=999.,mag2_cut=999.)
        else: x,y,mag1,mag2,qc_ind1,qc_ind2,matches = read_phot(photfiles[i])
        # HACK!!! This is to only use the phot stars found in both filters (with qc)
        if matched_only == True:
            qc_ind1 = matches
            qc_ind2 = matches
        d = pyfits.open(drzs[i])
        # get out of image coords
        ra,dec = XY2radec(d,x,y)
        # currect for chip gap and edges:
        within = in_chips(ra,dec,patch_list1a[i],patch_list1b[i])
        # field 3 is the chosen one to sit on top:
        if re.search('F03',photfiles[i]):
            inds1 = match_2inds(within,qc_ind1)
            inds2 = match_2inds(within,qc_ind2)
            if qc == False:
                inds1 = within
                inds2 = within
            Mag1 = mag1[inds1]
            Mag2 = mag2[inds2]
            Ra1 = ra[inds1]
            Dec1 = dec[inds1]
            Ra2 = ra[inds2]
            Dec2 = dec[inds2]
        # top row will only overlap with one to the right
        elif re.search('F0',photfiles[i]) or re.search('F05',photfiles[i]) or re.search('F06',photfiles[i]):
            not_overlapped = poly_over_under(ra,dec,patch_list1b[i-1],patch_list1a[i])
            inds1 = match_3inds(within,qc_ind1,not_overlapped)
            inds2 = match_3inds(within,qc_ind2,not_overlapped)
            if qc == False:
                inds1 = match_2inds(within,not_overlapped)
                inds2 = match_2inds(within,not_overlapped)
        # next two rows will overlap above and to the right
        else:
            # correct for overlap by field in above row
            #    by left chip
            top_a_a = poly_over_under(ra,dec,patch_list1a[i-4],patch_list1a[i])
            top_a_b = poly_over_under(ra,dec,patch_list1a[i-4],patch_list1b[i])
            top_b_b = poly_over_under(ra,dec,patch_list1b[i-4],patch_list1b[i])
            tops = list(set(top_a_a) & set(top_a_b) & set(top_b_b))
            # correct for overlap by field to the left, F09 is the left-most.
            if re.search('F09',photfiles[i]) or re.search('F15',photfiles[i]):
                inds1 = list(set(qc_ind1) & set(within) & set(tops))
                inds2 = list(set(qc_ind2) & set(within) & set(tops))
            else:
                top_left_b_a = poly_over_under(ra,dec,patch_list1b[i-5],patch_list1a[i])
                the_left = poly_over_under(ra,dec,patch_list1b[i-1],patch_list1a[i])
                inds1 = list(set(qc_ind1) & set(within) & set(tops) & set(top_left_b_a) & set(the_left))
                inds2 = list(set(qc_ind2) & set(within) & set(tops) & set(top_left_b_a) & set(the_left))
                if qc == False:
                    inds1 = list(set(within) & set(tops) & set(top_left_b_a) & set(the_left))
                    inds2 = list(set(within) & set(tops) & set(top_left_b_a) & set(the_left))
                    
        if plot==True: ax['fk5'].plot(ra[inds1],dec[inds1],'o',ms=4,mec=clrs[i],color=clrs[i])
        Mag1 = np.append(Mag1,mag1[inds1])
        Mag2 = np.append(Mag2,mag2[inds2])
        Ra1 = np.append(Ra1,ra[inds1])
        Dec1 = np.append(Dec1,dec[inds1])
        Ra2 = np.append(Ra2,ra[inds2])
        Dec2 = np.append(Dec2,dec[inds2])
        
    return Ra1,Dec1,Ra2,Dec2,Mag1,Mag2

def quality_cut(sharp,snr,mag,mag_cut=99.):
    # set the quality cut here, returns quality indices
    print 'using quality cut from PHATDataUtils.py'
    sharpsq = sharp*sharp
    return np.nonzero((sharpsq < 0.1) & (snr > 4.) & (mag < mag_cut))[0]

def quality_cut_multiple_detection(fitstable_data,sharp_cut=0.075,crowd_cut=1.0,snr_cut=4.,mag1_cut=99.,mag2_cut=99.,color_cut=[-99.,99.]):
    # dan's: crowd_cut = 1.0 -- default
    #        (sharp1 + sharp2)^2 < 0.075
    # ben's: crowd1 + crowd2 < 0.1
    #        (sharp1 + sharp2)^2 < 0.075
    sharp_crit = np.nonzero( ((fitstable_data.field('SHARP1') + fitstable_data.field('SHARP2')))**2. < sharp_cut)[0]
    crowd_crit = np.nonzero(fitstable_data.field('CROWD1') + fitstable_data.field('CROWD2') < crowd_cut)[0]
    snr_crit = np.nonzero((fitstable_data.field('SNR1')>snr_cut) & (fitstable_data.field('SNR2') > snr_cut))[0]
    mag_crit = np.nonzero( (fitstable_data.field('MAG1_UVIS') < mag1_cut) & (fitstable_data.field('MAG2_UVIS') < mag2_cut))[0]
    color = fitstable_data.field('MAG1_UVIS')-fitstable_data.field('MAG2_UVIS')
    color_crit = np.nonzero( (color >color_cut[0]) & (color <color_cut[1]) )[0]
    return list(set(sharp_crit) & set(crowd_crit) & set(snr_crit) & set(mag_crit) & set(color_crit) )

'''
def mag_errs(errmag,)

binsize = 0.1
mag1
m2 = mag2

clr = m1-m2

errclr = xmin+0.5
errmag = findgen(nint(nmbins/5)-1)
errcol = findgen(nint(nmbins/5)-1)
errmagerr = findgen(nint(nmbins/5)-1)
errcolerr = findgen(nint(nmbins/5)-1)

for q=0,n_elements(errmag)-1 do begin
    test = mmin + 5*(q+2)*binsize + 2.5*binsize
    test2 = where(mag2 gt test-2.5*binsize and mag2 le test+2.5*binsize and mag1
-mag2 gt -0.5 and mag1-mag2 lt 2.5)
    if test2[0] eq -1 then goto, jump1
    if n_elements(test2) lt 5 then goto, jump1
    errmag[q] = mmin + 5*(q+2)*binsize + 2.5*binsize
    errcol[q] = errclr 
   
     errmagerr[q] =  mean(mag2err(where(mag2 gt errmag[q]-2.5*binsize and mag2 le
  errmag[q]+2.5*binsize))^2)^0.5
  
  
  
  
    errcolerr[q] =  mean(mag1err(where(mag2 gt errmag[q]-2.5*binsize and mag2 le
  errmag[q]+2.5*binsize and mag1-mag2 gt -0.5 and mag1-mag2 lt 2.5))^2 + mag2err
(where(mag2 gt errmag[q]-2.5*binsize and mag2 le errmag[q]+2.5*binsize and mag1-
mag2 gt -0.5 and mag1-mag2 lt 2.5))^2)^0.5
  jump1:
endfor



'''
