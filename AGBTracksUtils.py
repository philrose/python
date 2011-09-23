import os,sys
import numpy as np
sys.path.append(os.environ['PYTHON_CODES'])
from GenUtils import get_afile
import time
import matplotlib.pyplot as plt

class AGBTracks(object):
    '''
    This is just like Table class from GenUtils but it also allows
    for a call to AGB step: see get_row_by_num
    '''
    def __init__(self, data_array, col_keys,name):
        self.data_array = data_array
        self.key_dict = dict(zip(col_keys,range(len(col_keys))))
        self.name = name
        
    def get_row(self,i):
        return self.data_array[i,:]
    
    def get_row_bynum(self,i):
        row = np.nonzero(track.data_array[:,track.key_dict['step']] == i)[0]
        return self.data_array[row,:]
    
    def get_col(self,key):
        return self.data_array[:,self.key_dict[key]]

def get_numeric_data(filename):
    f=open(filename,'r')
    col_keys = f.readline().replace('#','').replace('lg','').strip().split()
    f.close()
    try:
        data = np.loadtxt(filename)
    except ValueError:
        print 'problem with',filename
        data = np.zeros(len(col_keys))
    return AGBTracks(data,col_keys,filename)
        
def make_iso_file(track,Qs,slopes,isofile):
    '''
    t_min:          age in yr
    logl_min,       logL
    logte_min,      logTe
    mass_min,       actual mass along track
    mcore_min,      core mass
    co_min,         C/O ratio
    per_min,        period in days
    ip_min,         1=first overtone, 0=fundamental mode
    mlr_min,        - mass loss rate in Msun/yr
    logtem_min,     keep equal to logTe
    x_min,          X
    y_min,          Y
    xcno_min        X_C+X_O+X_N
    slope           dTe/dL
    '''
    # cull agb track to quiessent, write out.
    rows = [q-1 for q in Qs[:-1]] # cutting the final point for the file.

    f=open(track.name,'r').readlines()
    col_keys = f[0].replace('#','').replace('lg','').strip().split()

    #check on P0 == per_min
    outcols = ['age/yr','L_*','T_*','M_*','M_c','C/O','P0','P1','Pmod','dM/dt','T_*','H','Y']
    outcols.extend([key for key in col_keys if (key.startswith('C1') or key.startswith('N1') or key.startswith('O1'))])
    col_nums = [int(col_keys.index(key)) for key in outcols]
    cno =  col_nums[-7:]
    for r in rows:
        try:
            data_row = map(float,f[r].split())
        except TypeError:
            print r
            print ''
            print rows
            print ''
            print f[r-1]
            print track.name
            sys.exit()
        CNO = np.sum([data_row[c] for c in cno])
        mdot = 10**(data_row[col_nums[9]])
        if data_row[col_nums[8]] == 0:
            period = data_row[col_nums[6]]
        else:
            period = data_row[col_nums[7]]
        
        try:
            isofile.write('%.3f %.4f %.4f %.5f %.5f %.4f %.4e %i %.4e %.4f %.6e %.6e %.6e %.4f\n'%(data_row[col_nums[0]],data_row[col_nums[1]],data_row[col_nums[2]],data_row[col_nums[3]],data_row[col_nums[4]],data_row[col_nums[5]],period,data_row[col_nums[8]],mdot,data_row[col_nums[2]],data_row[col_nums[11]],data_row[col_nums[12]],CNO,1./slopes[list(rows).index(r)]))
        except IndexError:
            print list(rows).index(r)
            print len(rows),len(slopes)
            print 1./slopes[list(rows).index(r)]
    return

def make_cmd_input(cmd_input_name,tracce):
    # take default... should have this in mk_sims.py
    cmdin = open(cmd_input_name,'w')
    cmdin.write('\n')
    cmdin.write('1 isotrack/griglie_tutto_corr.dat isotrack/bassazams_mista.dat #kind_tracks, file_isotrack, file_lowzams \n')
    cmdin.write('3 %s # kind_tpagb, file_tpagb \n'%(tracce.split('1.3/')[-1]))
    cmdin.write('1 isotrack/final/pne_wd_test.dat # kind_postagb, file_postagb DA VERIFICARE file_postagb \n')
    cmdin.write('0 tab_ifmr/weidemann.dat # ifmr_kind, file with ifmr\n\n')
    cmdin.write('################################explanation######################\n')
    cmdin.write('kind_tracks: 1= normal file \n')
    cmdin.write('file_isotrack: tracks for low+int mass \n')
    cmdin.write('file_lowzams: tracks for low-ZAMS \n')
    cmdin.write('kind_tpagb: 0= notrack/griglie_tutto_corr.dat isotrack/bassazams_mista.dat #kind_tracks, file_isotrack, file_lowzams \n')
    cmdin.write('3 %s # kind_tpagb, file_tpagb \n'%(tracce))
    cmdin.write('1 isotrack/final/pne_wd_test.dat # kind_postagb, file_postagb DA VERIFICARE file_postagb \n')
    cmdin.write('0 tab_ifmr/weidemann.dat # ifmr_kind, file with ifmr\n\n')
    cmdin.write('##########################')
    cmdin.write('file_postagb: PN+WD tracks\n\n')
    cmdin.write('kind_ifmr: 0= default\n')
    cmdin.write('           1= from file\n')
    cmdin.write('\n')
    cmdin.close()
    return cmd_input_name

def make_met_file(tracce,Zs,Ys,isofiles):
    t = open(tracce,'w')
    t.write('%i\n'%len(isofiles))
    [t.write('%.4f\t%.3f\t%s\n'%(Zs[i],Ys[i],isofiles[i].split('1.3/')[-1])) for i in range(len(isofiles))]
    t.close()

def get_TP_inds(ntp):
    un =  np.unique(ntp)
    iTPs = [list(ntp).index(u) for u in un] # this is the first step in each TP.             
    TPs = [np.arange(iTPs[i],iTPs[i+1]) for i in range(len(iTPs)-1)] # The indices of each TP.
    TPs.append(np.arange(iTPs[i+1],len(ntp))) # don't forget the last one.
    return TPs

def get_unique_inds(ntp):
    un =  np.unique(ntp)
    iTPs = [list(ntp).index(u) for u in un] # this is the first step in each TP.             
    TPs = [np.arange(iTPs[i],iTPs[i+1]) for i in range(len(iTPs)-1)] # The indices of each TP.
    TPs.append(np.arange(iTPs[i+1],len(ntp))) # don't forget the last one.
    return TPs,iTPs


def add_points_to_q_track(track,qs):
    # when to add an extra point for low masses
    # if logt[qs+1] is hotter than logt[qs]
    # and there is a point inbetween logt[qs] and logt[qs+1] that is cooler
    # than logt[qs]
    # add the coolest point.
    addpt=[]
    logt = AGBTracks.get_col(track,'T_*')
    step = AGBTracks.get_col(track,'step')
    status =AGBTracks.get_col(track,'status')
    Tqs = logt[qs]
    # need to use some unique array, not log t, since log t could repeat,
    # index would find the first one, not necessarily the correct one.
    Sqs = step[qs]-1. # makes it the same as qs.
    # takes the difference in logt(qs) to see if we get hotter or colder.
    # ht = np.nonzero(np.concatenate(([-1.],np.diff(Tqs)))>0)[0]
    # finds where the logt goes from getting colder to hotter...
    ht = np.where(np.diff(np.sign(np.diff(Tqs))))[0]+1
    ht = np.append(ht,ht+1) # between the first and second
    Sqs_ht = Sqs[ht]
    # the indices between each hot point.
    t_mids = [map(int,step[int(Sqs_ht[i]):int(Sqs_ht[i+1])]) for i in range(len(Sqs_ht)-1)]
    Sqs_ht = Sqs_ht[:-1]
    for i in range(len(Sqs_ht)):
        hot_inds = np.nonzero(logt[int(Sqs_ht[i])] > logt[t_mids[i]])[0]
        if len(hot_inds) > 0:
            # index of the min T of the hot index from above.
            addpt.append(list(logt).index(np.min(logt[[t_mids[i][hot_ind] for hot_ind in hot_inds]])))
            
    if len(addpt)>0: addpt = [a for a in addpt if status[a]==7.]
    Qs = np.sort(np.concatenate((addpt,qs)))
    return Qs,addpt

def find_dldt(track,TPs,addpt):
    status = AGBTracks.get_col(track,'status')
    #dl/dT seems somewhat linear for 0.2<phi<0.4 ...
    phi = AGBTracks.get_col(track,'PHI_TP')
    phi[0]= -1. # The first line in the agb track is 1. This isn't the quiessent...
    lin_rise = np.nonzero((status==7) & (phi<0.4) & (phi>0.2))[0]
    rising = [list(set(TP) & set(lin_rise)) for TP in TPs] 
    logl = AGBTracks.get_col(track,'L_*')
    logt = AGBTracks.get_col(track,'T_*')
    order = 1
    fits = [np.polyfit(logt[r],logl[r],order) for r in rising if len(r)>0]
    slopes = [fits[i][0] for i in range(len(fits))]
    if len(addpt)>0:
        addrise = [TPs.index(TP) for TP in TPs if len(set(addpt) & set(TP)) > 0 ]
        addslope = slopes[addrise[0]]
        Slopes = []
        for s in slopes:
            if s == addslope:
                Slopes.append(addslope)
            Slopes.append(s)
    else:
        Slopes = slopes[:]
    return rising,Slopes,fits


def diag_plots(track,logl,logt,slopes,Qs,addpt,rising,fits,plotpath='DIAGPLOTS/'):
    if not os.path.isdir(plotpath): os.mkdir(plotpath)
    x = np.linspace(min(logt),max(logt))
    dd=plt.figure()
    dd=plt.title('%s %s M=%.2f'%(track.name.split('/')[-3],track.name.split('/')[-2],mass))
    dd=[plt.plot(x,fits[i][0]*x+fits[i][1],'--',color='blue',lw=0.5) for i in range(len(fits))]
    dd=plt.plot(logt,logl,color='black')
    dd=[plt.plot(logt[r],logl[r],color='red',lw=2) for r in rising if len(r)>0]
    dd=plt.plot(logt[map(int,Qs)],logl[map(int,Qs)],color='green',lw=2)
    dd=[plt.plot(logt[q],logl[q],'o',color='green') for q in Qs]
    dd=[plt.plot(logt[add],logl[add],'o',color='purple') for add in addpt if len(addpt)>0]
    dd=plt.axis([max(logt)+0.01,min(logt)-0.01,min(logl)-0.01,max(logl)+0.01])
    plt.savefig(plotpath+'diag_'+track.name.split('/')[-1].replace('.dat','.png'))
    plt.close()
    #dd=plt.axis([3.74,3.35,2.4,4.5])
    plt.figure()
    plt.plot(age,logt)
    plt.plot(age[Qs],logt[Qs],'o')
    plt.xlabel('age')
    plt.ylabel('Log T')
    plt.savefig(plotpath+'diag_'+track.name.split('/')[-1].replace('.dat','_age_T.png'))
    plt.close()
    
    plt.figure()
    plt.plot(age,logl)
    plt.plot(age[Qs],logl[Qs],'o')
    plt.xlabel('age')
    plt.ylabel('Log L')
    plt.savefig(plotpath+'diag_'+track.name.split('/')[-1].replace('.dat','_age_L.png'))
    plt.close()
    
    return

def make_readme(agb_iso_track,agbmix):
    readmef = agb_iso_track+'readme.txt'
    if os.path.isfile(readmef): 
        readme = open(agb_iso_track+'readme.txt','a')
        readme.write('\n --- %s --- \n'%time.strftime("%Y-%m-%d %H:%M:%S"))
    else:
        readme = open(agb_iso_track+'readme.txt','w')
    readme.write('This directory has agb tracks parsed from AGBTracksUtils. Sources are \n %s \n and within\n'%agbmix)
    readme.close()

def plot_ifmr(imfrfile,name):
    mi,mf,z = np.loadtxt(imfrfile,unpack=True)
    zinds,unz = get_unique_inds(z)
    [plt.plot(mi[zinds[i]],mf[zinds[i]],label=str(z[unz[i]])) for i in range(len(unz))]
    [plt.plot(mi[zinds[i]],mf[zinds[i]],'o') for i in range(len(unz))]
    plt.legend(loc=2)
    plt.title(name)
    plt.xlabel(r'$M_i/M_{\odot}$')
    plt.ylabel(r'$M_f/M_{\odot}$')
    plt.savefig(imfrfile.replace('dat','png'))
    return

if __name__ == "__main__":
    
    #sources:
    agbmix = '/Users/Phil/research/Italy/WFC3SNAP/Sp2011/fromPaola/agb_full_wv/'
    name = 'agb_wv'
    ifmr = agbmix+name+'_ifmr.dat'
    
    # where to put stuff!
    trilegal_dir = '/Users/Phil/research/Italy/trilegal_1.3/'
    cmd_input_name = trilegal_dir+'cmd_input_'+name+'.dat'
    
    isotrack_dir = trilegal_dir+'isotrack/'
    tracce = isotrack_dir+'tracce_'+name+'.dat'
    # tracks, make the dir if needed.
    agb_iso_track = isotrack_dir+agbmix.split('/')[-2]+'/'
    if not os.path.isdir(agb_iso_track): 
        os.mkdir(agb_iso_track)
        print 'made dir:',agb_iso_track
    
    
    make_readme(agb_iso_track,agbmix)
    mets = [l for l in os.listdir(agbmix) if os.path.isdir(l) and l.startswith('Z')]
    mfile = open(ifmr,'w')
    mfile.write('# M_i M_f Z\n')
    
    isofiles,Zs,Ys = [],[],[]
    for met in mets:
        isofile = agb_iso_track+(met+'_'+name+'.dat').lower()
        out = open(isofile,'w')
        filenames = get_afile(agbmix+met+'/','*.dat')
        print 'found',len(filenames),'with metallicity',met
        for filename in filenames:
            addpt = []
            track =  get_numeric_data(filename)
            ntp = AGBTracks.get_col(track,'NTP')
            age = AGBTracks.get_col(track,'age/yr')
            phi = AGBTracks.get_col(track,'PHI_TP')
            logt = AGBTracks.get_col(track,'T_*')
            logl = AGBTracks.get_col(track,'L_*')
            
            TPs = get_TP_inds(ntp)
            phi[0]= -1. # The first line in the agb track is 1. This isn't the quiessent...
            
            # The Quessent phase is the the max phase in each TP, i.e., closest to 1.
            qs = [list(phi).index(np.max(phi[TP])) for TP in TPs]
            Qs = qs
            
            mass = float(filename.split('/')[-1].split('_')[1])
            if len(Qs) <= 9 and mass < 3.: Qs,addpt = add_points_to_q_track(track,qs)
            Qs = map(int,Qs)
            rising,slopes,fits = find_dldt(track,TPs,addpt)
            
            make_iso_file(track,Qs,slopes,out)
            diag_plots(track,logl,logt,slopes,Qs,addpt,rising,fits,plotpath=met+'/DIAGPLOTS/')
            M_s = AGBTracks.get_col(track,'M_*')
            mfile.write('%f %f %f\n'%(M_s[0],M_s[-1],float(met.replace('Z',''))))
            
        out.close()
        print 'wrote',isofile
        isofiles.append(isofile)
        Ys.append(AGBTracks.get_col(track,'Y')[0])
        Zs.append(float(met.replace('Z','')))
    
    metfile = make_met_file(tracce,Zs,Ys,isofiles)
    print 'wrote',tracce
    cmd_input = make_cmd_input(cmd_input_name,tracce)
    print 'wrote',cmd_input_name
    mfile.close()
    print 'wrote',ifmr
    plot_ifmr(ifmr,name)

'''
plt.figure(100)
dd=plt.plot(logt[map(int,Qs)],logl[map(int,Qs)],color='black',lw=2,label='%.2f'%(mass))

dd=plt.title('%s %s'%(filename.split('/')[-3],filename.split('/')[-2]))
dd=plt.axis([3.74,3.35,2.8,4.5])
plt.savefig('_'.join((filename.split('/')[-3],filename.split('/')[-2]))+'.png')
'''
