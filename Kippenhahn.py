#
#  Kippenhahn.py
#  
#
#  Created by Philip Rosenfield on 10/5/10.
#
# python Kippenhahn.py track.PMS
python_code_base = ('/Users/Phil/Dropbox/research/python')
import sys
sys.path.append(python_code_base)
from PadovaTracksUtils import *
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.axes_grid1 import host_subplot
import mpl_toolkits.axisartist as AA

tracks = ['/Users/Phil/Dropbox/research/Italy/tracks/lmcmix/PHIL_Z0.001_Y0.247/Z0.001Y0.247OUTA1.72ALFO0_F7_M5.00.PMS']
tracks.append('/Users/Phil/Dropbox/research/Italy/tracks/lmcmix/PHIL_Z0.001_Y0.247/Z0.001Y0.247OUTA1.72ENV_1.2ALFO0.5_F7_M5.00.PMS')
#track = '/Users/Phil/Dropbox/research/Italy/tracks/lmcmix/PHIL_Z0.001_Y0.247/Z0.001Y0.247OUTA1.72ENV_1.0ALFO0.5_F7_M5.00.PMS'
#track = '/Users/Phil/Dropbox/research/Italy/tracks/lmcmix/PHIL_Z0.001_Y0.247/Z0.001Y0.247OUTA1.72ENV_0.5ALFO0.5_F7_M5.00.PMS'

age1 = 7.9
age2 = 8.001
cols = ['black','blue']
#track = sys.argv[1]
for track,col in zip(tracks,cols):
    t = 0
    t = get_tracks(track)
    info = info_from_track_filename(track)
    
    age = Tracks.get_stable_col(t,'AGE')
    lage = np.log10(age)
    lx = Tracks.get_stable_col(t,'LX')
    qh1 = Tracks.get_stable_col(t,'QH1')
    qh2 = Tracks.get_stable_col(t,'QH2')
    ly = Tracks.get_stable_col(t,'LY')
    ycen = Tracks.get_stable_col(t,'YCEN')
    qhe1 = Tracks.get_stable_col(t,'QHE1')
    qhe2 = Tracks.get_stable_col(t,'QHE2')
    lc = Tracks.get_stable_col(t,'LC')
    qc1 = Tracks.get_stable_col(t,'QC1')
    qc2 = Tracks.get_stable_col(t,'QC2')
    conv = Tracks.get_stable_col(t,'CONV')
    ci1 = np.abs(Tracks.get_stable_col(t,'CI1'))
    cf1 = np.abs(Tracks.get_stable_col(t,'CF1'))
    ci2 = np.abs(Tracks.get_stable_col(t,'CI2'))
    cf2 = np.abs(Tracks.get_stable_col(t,'CF2'))
    ci3 = np.abs(Tracks.get_stable_col(t,'CI3'))
    cf3 = np.abs(Tracks.get_stable_col(t,'CF3'))
    ci4 = np.abs(Tracks.get_stable_col(t,'CI4'))
    cf4 = np.abs(Tracks.get_stable_col(t,'CF4'))
    ci5 = np.abs(Tracks.get_stable_col(t,'CI5'))
    cf5 = np.abs(Tracks.get_stable_col(t,'CF5'))
    qdisc = Tracks.get_stable_col(t,'QDISC')
    qinte = Tracks.get_stable_col(t,'QINTE')
    qhel  = Tracks.get_stable_col(t,'QHEL')
    logL = Tracks.get_stable_col(t,'LOG_L')
    logTe = Tracks.get_stable_col(t,'LOG_TE')
    he =  np.nonzero(ly > 0.0)[0]
    heb = np.nonzero((ly > 0.0) & (ycen > 0.0))[0]
    heb2 = np.nonzero((ly > 0.0) & (ycen == 0))[0]
    cfus = np.nonzero(lc > 0.0)[0]
    
    inds = [heb,heb2,cfus]
    #cols = ['red','purple','green']
    marks = ['-.','-','--']
    names = ['He core','He shell','C+O']
    pt1 = heb[0]
    #closest_match(num,somearray)
    pt2 = closest_match(logL[heb[pt1:]].min(),logL[heb])
    pt3 = logTe[heb].argmax()
    # pt4 = start of He out of core? 
    pt5 = heb2[0]
    ageBlueLoop = 10**lage[pt5] - 10**lage[pt1]
    print ageBlueLoop/age[pt1]
    
    strage = str('%.3e' % ageBlueLoop).replace('e+06',' Myr').replace('e+07','0 Myr').replace('e+08','00 Myr')
    #plt.figure(1)
    #for ind,mark,name in zip(inds,marks,names):
    #    plt.plot(logTe[ind],logL[ind],mark,linewidth=3,label=name,color=col)
    
    #plt.annotate('A',xy=(logTe[pt1],logL[pt1]),color=col)
    #plt.annotate('B',xy=(logTe[heb[pt2]],logL[heb[pt2]]),color=col)
    #plt.annotate('C',xy=(logTe[heb[pt3]],logL[heb[pt3]]),color=col)
    #plt.annotate('D',xy=(logTe[pt5],logL[pt5]),color=col)
    #try:
        #plt.title('Z = %.3f Y=%.3f ENV = %.1f M = %.2f D-A = %s' % (info['Z'],info['Y'],info['ENV'],info['M'], strage))
    #except KeyError:
        #plt.legend(loc=0)
        #plt.axis([logTe[he].max()+0.01,logTe[he].min(),logL[he].min(),logL[he].max()])
        #plt.title('Z = %.3f Y=%.3f ENV = 0.0 M = %.2f D-A = %s' % (info['Z'],info['Y'],info['M'], strage))
        #pass
    
    num = cols.index(col)+211
    host = host_subplot(num,axes_class=AA.Axes)
    plt.subplots_adjust(right=0.75)
    par1 = host.twinx()
    par2 = host.twinx()
    
    #plt.subplot(2,1,cols.index(col)+1)
    offset = 60
    new_fixed_axis = par2.get_grid_helper().new_fixed_axis
    par2.axis["right"] = new_fixed_axis(loc="right",
                                        axes=par2,
                                        offset=(offset, 0))
    par2.axis["right"].toggle(all=True)
    
    host.set_xlim(age1, age2) # age range
    host.set_ylim(0, 1) # m/M
    host.set_xlabel("Log Age (yrs)")
    host.set_ylabel("M/Mstar")
    par1.set_ylabel("Log L",color='green')
    par2.set_ylabel("Log Te",color='red')
    
    #plt.ylabel('M/Mstar')
    #plt.axis([lage.min(),lage.max(),0,1])
    #plt.axis([age1,age2,0,1])
    #try:
        #plt.title('Z = %.3f Y=%.3f ENV = %.1f M = %.2f D-A = %s' % (info['Z'],info['Y'],info['ENV'],info['M'], strage))
    #except KeyError:
        #plt.legend(loc=0)
        #plt.title('Z = %.3f Y=%.3f ENV = 0.0 M = %.2f D-A = %s' % (info['Z'],info['Y'],info['M'], strage))
    
    #plt.figure(2)
    p2, = host.plot(lage,logL,linewidth=3,color='green')
    p3, = host.plot(lage,logTe,linewidth=3,color='red')
    
    #plt.ylabel('Log L/Lsun')
    #plt.setp(plt.gca(), 'xticklabels', [])
    #ax.axis([age1,age2,logL.min(),logL.max()])
    #ax1 = ax.twinx()
    par1.set_ylim(logL.min(),logL.max())
    par2.set_ylim(logTe.min(),logTe.max())
    #host.legend()
    
    host.axis["left"].label.set_color('black')
    par1.axis["right"].label.set_color(p2.get_color())
    par2.axis["right"].label.set_color(p3.get_color())
    for t in par1.get_yticklabels():
        t.set_color('green')
    for t in par2.get_yticklabels():
        t.set_color('red')
        
    #plt.ylabel('Log Te', color='red')
    plt.draw()
    
    #for t in ax1.get_yticklabels():
    #    t.set_color('red')
    #plt.ylabel('Log Te', color='red')
    #plt.setp(plt.gca(), 'xticklabels', [])
    #ax1.axis([age1,age2,3.6,4.4])

plt.show()
    #print 'Wrote '+track+'_kip.png'
