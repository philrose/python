#
#  plotUVrealquick.py
#  
#
#  Created by Philip Rosenfield on 11/5/10.
#  Copyright (c) 2010 __MyCompanyName__. All rights reserved.
#
import sys,re,pyfits
sys.path.append('/Users/Phil/Dropbox/research/python')
from GenUtils import inside
def finite_mags_6filts(data,cam):
    # either takes cam like 'acs' and gives mag1,mag2
    # or takes filter like F275W and returns mag            
    if len(cam) <= 3:
        mag1 = data.field(cam+'_mag1')
        mag2 = data.field(cam+'_mag2')
        # wanted to keep em all the same length, 
        # so just figured 999 would be out of the 
        # plotting regime...
        nans = np.isnan(mag2)
        mag2[nans] = 999.
        mag1[nans] = 999.
        return mag1,mag2
    else:
        cam = filt_to_mag(cam)
        mag = data.field(cam)
        nans = np.isnan(mag)
        mag[nans] = 999.
        return mag
def filt_to_mag(filt):
    if filt == 'F275W': cam = 'uv_mag1'
    if filt == 'F336W': cam = 'uv_mag2'
    if filt == 'F475W': cam = 'acs_mag1'
    if filt == 'F814W': cam = 'acs_mag2'
    if filt == 'F110W': cam = 'ir_mag1'
    if filt == 'F160W': cam = 'ir_mag2'
    try:
        x =str(cam)
    except ValueError:
        print 'Filt ',filt,' not found.'
        cam = 0
    return cam

def load_verts(filename):
    return np.loadtext(filename,unpack=True)
 


fits = pyfits.open('/Users/Phil/Dropbox/research/PHAT/dustin_corr3_all.fits')
data = fits[1].data

filtset1 = 'F475W-F814W'
filtset2 = 'F814W-F110W'
filt1 = filtset1.split('-')[0]
filt2 = filtset1.split('-')[1]
filt3 = filtset2.split('-')[0]
filt4 = filtset2.split('-')[1]
Mag1 = finite_mags_6filts(data,filt1)
Mag2 = finite_mags_6filts(data,filt2)
Mag3 = finite_mags_6filts(data,filt3)
Mag4 = finite_mags_6filts(data,filt4)
Color1 = Mag1-Mag2
Color2 = Mag3-Mag4
ColorAxis = [min(Color1[(np.nonzero(Color1>-100))]),
             max(Color1[(np.nonzero(Color1>-100))]),
             min(Color2[(np.nonzero(Color2>-100))]),
             max(Color2[(np.nonzero(Color2>-100))])]

p = plt.scatter(Color1,Color2,marker='o',s=1.5,edgecolor='grey',facecolor='grey',zorder=1)

points =  np.column_stack((Color1,Color2))
xs,ys = load_verts('/Users/Phil/Dropbox/research/PHAT/PHAT6_F475W-F814W_F814W-F110W_1.verts')
sel_ind = inside(xs,ys,Color1,Color2,verts=True)
p = plt.scatter(Color1[sel_ind],Color2[sel_ind],s=1.5,color='red',marker='o',label='N = '+str(sum((Color1_sel>-100)*1)),zorder=10)


plt.show()