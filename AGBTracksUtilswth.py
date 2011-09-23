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
            print data_row[col_nums[0]]
            print data_row[col_nums[1]]
            print data_row[col_nums[2]]
            print data_row[col_nums[3]]
            print data_row[col_nums[4]]
           print  data_row[col_nums[5]]
           print period
           print data_row[col_nums[8]]
          print mdot
          print data_row[col_nums[2]]
         print data_row[col_nums[11]]
         print data_row[col_nums[12]]
           print CNO
         print 1./slopes[list(rows).index(r)]
    return

def make_cmd_input(cmd_input_name,tracce):
    # take default... should have this in mk_sims.py
    cmdin = open(cmd_input_name,'w')
    cmdin.write('\n')
    cmdin.write('1 is    otrack/griglie_tutto_corr.dat isotrack/bassazams_mista.dat #kind_tracks, file_isotrack, file_lowzams \n')
cmdin.write('3 %s # kind_tpagb, file_tpagb \n'%(tracce))
cmdin.write('1 isotrack/final/pne_wd_test.dat # kind_postagb, file_postagb DA VERIFICARE file_postagb \n')
cmdin.write('0 tab_ifmr/weidemann.dat # ifmr_kind, file with ifmr\n\n')
cmdin.write('################################explanation######################\n')
cmdin.write('kind_tracks: 1= normal file \n')
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

    def make_met_file(tracce,Zs,Ys,isofiles):￼

    def imfr():
    # make initial to final mass relationship
    prinn.write('3 %s # kind_tpagb, file_tpagb \n'%(tracce))
    cmdin.write('1 isotrack/final/pne_wd_test.dat # kind_postagb, file_postagb DA VERIFICARE file_postagb \n')
    cmdin.write('0 tab_ifmr/weidemann.dat # ifmr_kind, file with ifmr\n\n')
    cmdin.write('##########################')
    cmdin.write('file_postagb: PN+WD tracks\n\n')
    cmdin.write('kind_ifmr: 0= default\n')
    cmdin.write('           1= from file\n')
    cmdin.write('\n')
    cmdin.close()
return cmd_input_name

def make_met_file(tracce,Zs,Ys,isofiles):￼

def imfr():
# make initial to final mass relationship
print 'not done!'

def get_TP_inds(NTP):￼

def add_points_to_q_track(track,qs):￼

def find_dldt(track,TPs,addpt):￼

def diag_plots(track,logl,log    t,slopes,Qs,addpt,rising,fits,plotpath='DIAGPLOTS'):
    x = np.linspace(min(logt),max(logt))
    dd=plt.figure()
    dd=plt.title('%s %s M=%.2f'%(filename.split('/')[-3],filename.split('/')[-2],mass))
    dd=[plt.plot(x,fits[i][0]*x+fits[i][1],'--',color='blue',lw=0.5) for i in range(len(fits))]
    dd=plt.plot(logt,logl,color='black')
    dd=[plt.plot(logt[r],log    
    plt.figure()
    plt.plot(age,logt)
    plt.plot(age[Qs],logt[Qs],'o')
    plt.xlabel('age')
    plt.ylabel('Log T')
    plt.savefig(plotpath+'/diag_'+track.name.split('/')[-1].replace('.dat','age_T.png'))
    plt.close()
    
    plt.figure()
    plt.plot(age,logl)
    plt.plot(age[Qs],logl[Qs],'o')
    plt.xlabel('age')
    plt.ylabel('Log L')
    plt.savefig(plotpath+'/diag_'+track.name.split('/')[-1].replace('.dat','age_L.png'))
    plt.close()
    
    return


if __name__ == "__main__":
    
    #sources:
    agbmix = '/Users/Phil/research/Italy/WFC3SNAP/Sp2011/fromPaola/agb_full_wv/'
    name = 'agb_wv'
    
    # where to put stuff!
    trilegal_dir = '/Users/Phil/research/Italy/trilegal_1.3/'
    isotrack_dir = trilegal_dir+'isotrack/'
    agb_iso_track = isotrack_dir+agbmix.split('/')[-2]+'/'
    if not os.path.isdir(agb_iso_track): 
        os.mkdir(agb_iso_track)
        print 'made dir:',agb_iso_track
    
    readmef = agb_iso_track+'readme.txt'
    if os.path.isfile(readmef): 
        readme = open(agb_iso_track+'readme.txt','a')
        readme.write('\n --- %s --- \n'%time.strftime("%Y-%m-%d %H:%M:%S"))
    else:
        readme = open(agb_iso_track+'readme.txt','w')
    readme.write('This directory has agb tracks parsed from AGBTracksUtils. Sources are \n %s \n and within\n'%agbmix)
    readme.close()
    
    mets = [l for l in os.listdir(agbmix) if os.path.isdir(l) and l.startswith('Z')]
    
    isofiles,Zs,Ys = [],[],[]
    for met in mets:
        if met != 'Z0.019': continue
        isofile = agb_iso_track+(met+'_'+name+'.dat').lower()
        out = open(isofile,'w')
        filenames = get_afile(agbmix+met+'/','*.dat')
        print 'found',len(filenames),'with metallicity',met
        for filename in filenames:
            print filename
            skip=0
            addpt = []
            track =  get_numeric_data(filename)
            try: 
                ntp = AGBTracks.get_col(track,'NTP')
            except IndexError:
                print 'skipping',filename
                skip=1
                break
            
            TPs = get_TP_inds(ntp)
            
            age = AGBTracks.get_col(track,'age/yr')
            phi = AGBTracks.get_col(track,'PHI_TP')
            phi[0]= -1. # The first line in the agb track is 1. This isn't the quiessent...
            
            # The Quessent phase is the the max phase in each TP, i.e., closest to 1.
            qs = [list(phi).index(np.max(phi[TP])) for TP in TPs]
            Qs = qs
            logt = AGBTracks.get_col(track,'T_*')
            logl = AGBTracks.get_col(track,'L_*')
            
            mass = float(filename.split('/')[-1].split('_')[1])
            if len(Qs) <= 9 and mass < 3.: Qs,addpt = add_points_to_q_track(track,qs)
            Qs = map(int,Qs)
            rising,slopes,fits = find_dldt(track,TPs,addpt)
            
            make_iso_file(track,Qs,slopes,out)
            diag_plots(track,logl,logt,slopes,Qs,addpt,rising,fits,plotpath=met+'/DIAGPLOTS')
            plt.figure(100)
            dd=plt.plot(logt[map(int,Qs)],logl[map(int,Qs)],color='black',lw=2,label='%.2f'%(mass))
        
        out.close()    
        print 'wrote',isofile
        dd=plt.title('%s %s'%(filename.split('/')[-3],filename.split('/')[-2]))
        dd=plt.axis([3.74,3.35,2.8,4.5])
        plt.savefig('_'.join((filename.split('/')[-3],filename.split('/')[-2]))+'.png')
        if skip==0:
            isofiles.append(isofile)
            Ys.append(AGBTracks.get_col(track,'Y')[0])
            Zs.append(float(met.replace('Z','')))
        
    tracce = 'tracce_'+name+'.dat'
    metfile = make_met_file(tracce,Zs,Ys,isofiles)
    print 'wrote',tracce
    cmd_input_name = 'cmd_input_'+name+'.dat'
    cmd_input = make_cmd_input(cmd_input_name,tracce)
    print 'wrote',cmd_input_name





'''
    dd=plt.plot(age,logl,color='black')
    dd=[plt.plot(age[r],logl[r],color='red',lw=2) for r in rising if len(r)>0]
    dd=[plt.plot(age[q],logl[q],'o',color='green') for q in Qs]
    dd=plt.plot(age[addpt],logl[addpt],'o',color='purple')
    plt.axis([0,1,np.min(logl),np.max(logt)])

 
    
def get_agb_tracks(filename):
     As opposed to AGBTracks, this will just read into dictionary
    f=open(filename,'r')
    col_keys = f.readline().replace('#','').replace('lg','').strip().split()
    f.close()
    
    data = np.loadtxt(filename,unpack=True)
    data_array ={}
    for i,key in enumerate(col_keys):
        data_array[key] = data[i]
    
    return data_array
   

'age/yr','L_*','T_*','M_*','M_c','C/O','P0','P1','Pmod','dM/dt','T_*','H','Y','CNO'

10**'dM/dt'


def write_quies_track()
step, status, ntp, age, M, dM/dt, logl, L_H, logt, M_c, Y, Z, phi, C/O, Tbot, Pmod, P, DUST = np.loadtxt(agbtrack,unpack=True)

out = open(outfile,'w')
out.write(
''%())
 0.000 2.9643 3.6538 0.55000 0.49450 0.4784 1.7109E+01 1 2.0893e-08 3.6538 7.629569E-01 2.339915E-01 0.000708379
age,logl,logt,M,M_c,co,P,Pmod,Mdot,logt,1-Y-Z,Y,CNO??

#step status NTP       age/yr        M_* lg dM/dt     lg L_*   lg L_H   lg T_*        M_c    Y       Z     PHI_TP    C/O     Tbot    Pmod     P    DUST

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
xcno_min ;      X_C+X_O+X_N    



CNO = [key for key in track.key_dict if (key.startswith('C1') or key.startswith('N1') or key.startswith('O1'))]






1. Seperate by TP
2. Calculate Quies.
3. Calculate dL/dT
4. Write out necessary line.

#diag plots
x = np.linspace(min(logt),max(logt))
plt.plot(logt,logl,color='black')
[plt.plot(x,fits[i][0]*x+fits[i][1],color='blue',lw=2) for i in range(len(fits))]
[plt.plot(logt[r],logl[r],color='red',lw=2) for r in rising]
plt.axis([3.65,3.5,2,4])


    
'''