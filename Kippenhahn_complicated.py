#
#  Kippenhahn.py
#  
#
#  Created by Philip Rosenfield on 10/5/10.
#
# python Kippenhahn.py track.PMS
import sys,os
sys.path.append(os.environ['PYTHON_CODES'])
from PadovaTracksUtils import *
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.axes_grid1 import host_subplot
import mpl_toolkits.axisartist as AA

#tracks = ['/Users/Phil/Dropbox/research/Italy/tracks/lmcmix/PHIL_Z0.001_Y0.247/Z0.001Y0.247OUTA1.72ALFO0_F7_M5.00.PMS']
#tracks.append('/Users/Phil/Dropbox/research/Italy/tracks/lmcmix/PHIL_Z0.001_Y0.247/Z0.001Y0.247OUTA1.72ENV_1.2ALFO0.5_F7_M5.00.PMS')
#tracks.append('/Users/Phil/Dropbox/research/Italy/tracks/lmcmix/PHIL_Z0.001_Y0.247/Z0.001Y0.247OUTA1.72ENV_1.0ALFO0.5_F7_M5.00.PMS')
#tracks.append('/Users/Phil/Dropbox/research/Italy/tracks/lmcmix/PHIL_Z0.001_Y0.247/Z0.001Y0.247OUTA1.72ENV_0.5ALFO0.5_F7_M5.00.PMS')
tracks = ['/Users/Phil/Dropbox/research/Italy/tracks/lmcmix/PHIL_Z0.002_Y0.252/Z0.002Y0.252OUTA1.72ALFO0_F7_M7.00.PMS']
#tracks.append('/Users/Phil/Dropbox/research/Italy/tracks/lmcmix/PHIL_Z0.002_Y0.252/Z0.002Y0.252OUTA1.72ENV_1.5ALFO0.5_F7_M7.00.PMS')
def kippenhahn_complicated(track_obj):
    t = track_obj
    #age1 = 7.9
    #age2 = 8.001
    age1 = 0.8
    age2 = 1
    axisNum =1
    info = info_from_track_filename(track_obj.name)
    try:
        env = info['ENV']
    except KeyError:
        env = 0.0
    age = Tracks.get_stable_col(t,'AGE')
    #lage = np.log10(age)
    lage = age/age.max()
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
    logrhc = Tracks.get_stable_col(t,'LOG_RHc')
    logTc =  Tracks.get_stable_col(t,'LOG_Tc')
    logPc =  Tracks.get_stable_col(t,'LOG_Pc')
    he =  np.nonzero(ly > 0.0)[0]
    heb = np.nonzero((ly > 0.0) & (ycen > 0.0))[0]
    heb2 = np.nonzero((ly > 0.0) & (ycen == 0))[0]
    cfus = np.nonzero(lc > 0.0)[0]
    
    pt1 = heb[0]
    pt2 = closest_match(logL[heb[pt1:]].min(),logL[heb])
    pt3 = logTe[heb].argmax()
    # pt4 = start of He out of core? 
    pt5 = heb2[0]
    ageBlueLoop = age[pt5] - age[pt1]
    print ageBlueLoop/age[pt1]
    
    strage = str('%.3e' % ageBlueLoop).replace('e+06',' Myr').replace('e+07','0 Myr').replace('e+08','00 Myr')
    
    num = len(tracks)*100+10+axisNum
    host = host_subplot(num,axes_class=AA.Axes)
    plt.subplots_adjust(right=0.6)
    par1 = host.twinx()
    par2 = host.twinx()
    
    offset = 60
    new_fixed_axis = par2.get_grid_helper().new_fixed_axis
    par2.axis["right"] = new_fixed_axis(loc="right",
                                        axes=par2,
                                        offset=(offset, 0))
    par2.axis["right"].toggle(all=True)
    host.set_title('Mass = %.2f Env = %.2f' % (info['M'],env))
    host.set_xlim(age1, age2) # age range
    host.set_ylim(0, 1) # m/M
#    host.set_xlabel("Log Age (yrs)")
    host.set_xlabel("per cent total lifetime")
    host.set_ylabel("M/Mstar")
    par1.set_ylabel("Log L",color='orange')
    par2.set_ylabel("Log T",color='purple')
    
    p1, = host.plot(lage,qh1,color='blue')
    p1, = host.plot(lage,qh2,color='blue',label='H fusion')
    p1, = host.plot(lage,qhe1,color='red')
    p1, = host.plot(lage,qhe2,color='red',label='He fusion')
    host.fill_between(lage,qhe1,qhe2,color='red',alpha=0.5)
    p1, = host.plot(lage,qc1,color='green')
    p1, = host.plot(lage,qc2,color='green',label='C+O fusion')
    host.fill_between(lage,qc1,qc2,color='green',alpha=0.5)
    host.fill_between(lage,qh1,qh2,color='blue',alpha=0.5)
    p1, = host.plot(lage,conv,color='grey',label='Conv zone')
    p1, = host.plot(lage,cf1,color='grey')
    p1, = host.plot(lage,qdisc,'-.',color='black',label='H discont.')
    p1, = host.plot(lage,qinte,'-',color='black',label='H shell base')
    p1, = host.plot(lage,qhel,'--',color='black',label='H depleted core')
    p1, = host.plot(lage,ycen,color='red',label='Y in core')
    host.fill_between(lage,cf1,ci1,color='grey',alpha=0.5)
    host.fill_between(lage,cf2,ci2,color='grey',alpha=0.5)
    host.fill_between(lage,cf3,ci3,color='grey',alpha=0.5)
    host.fill_between(lage,cf4,ci4,color='grey',alpha=0.5)
    host.fill_between(lage,cf5,ci5,color='grey',alpha=0.5)
    p2, = par1.plot(lage,logL,linewidth=3,color='orange',label='Log L')
    p2, = par1.plot(lage,np.log10(lx*10**logL),'--',linewidth=3,color='orange',label='Lx')
    p2, = par1.plot(lage,np.log10(ly*10**logL),'-.',linewidth=3,color='orange',label='Ly')
    p3, = par2.plot(lage,logTe,linewidth=3,color='purple',label='Log Te')
    p3, = par2.plot(lage,logTc,'--',linewidth=3,color='purple',label='Log Tc')    
    par1.set_ylim(logL.min(),logL.max())
    par2.set_ylim(logTe.min(),logTc.max())
    leg = host.legend(loc=2,bbox_to_anchor=(1.3,1))
    leg._drawFrame=False
    host.axis["left"].label.set_color('black')
    par1.axis["right"].label.set_color(p2.get_color())
    par2.axis["right"].label.set_color(p3.get_color())
    for t in par1.get_yticklabels():
        t.set_color('green')
    for t in par2.get_yticklabels():
        t.set_color('red')
        
    plt.draw()
    plt.show()
    #print 'Wrote '+track+'_kip.png'
