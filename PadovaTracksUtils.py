#
#  PadovaTracksUtils.py
#  
import numpy as np
import sys,os,re

class TrilegalTab:
	def __init__(self,
					col_heads,  #col_heads
					data_array, #slice of data array returned by create_data_array
					col_keys,   #dictionary returned by get_col_keys
				):
                self.data_array = data_array
                self.key_dict = dict(zip(col_keys,range(len(col_keys))))
		
		def get_row(self,i):
			return self.data_array[i,:]
        
		def get_col(self,key):
			return self.data_array[:,self.key_dict[key]]

def get_data(filename):
	f = open(filename,'r')
	lines = f.readlines()
	f.close()
	d = {}
	col_heads = lines[0].replace('#','')
	col_keys = col_heads.split()
	for key in col_keys:
		d[key]=[]
	for i in range(len(lines)):
		if lines[i].startswith('#'): continue
		data = lines[i].split()
		for key,i in zip(col_keys,range(len(col_keys))):
			d[key].append(is_numeric(data[i]))
	return d
    
def read_tagged_data(filename):
    # Hopefully you pasted the line with the fiter names on the header.
    f = open(filename,'r')
    lines = f.readlines()
    f.close()
    i = 0
    d = {}
    for line in lines:
        if line.startswith('#'):
            line = line.replace('#','')
            col_keys = line.split()
            i += 1
    for key in col_keys:
        d[key] = []
    if i == 0: 
        print 'I don''t know the names of the damn filters.' 
    for line in lines:
        if line.startswith('#'): continue
        row = line.split()
        for key,i in zip(col_keys,range(len(col_keys))):
            d[key].append( is_numeric(row[i]))
    return d

def track_string(track):
    better = track.split('tracks/')[-1].replace('/',' ').replace('PHIL','').replace('_',' ').replace('Z','Z=').replace('Y','Y=').strip()
    if re.search('C12O16',better):
        better = better.replace('C12O16','')
        rate = 'cf88'
    else:
        rate = 'bu96'
    better += ' C12O16: '+rate
    return better

def info_from_track_filename(filename):
    # just gives info from tracks, 
    # if you want to have env=0 in the dict, use 
    # info_from_track_filename_env 
    # which calls this.
    filename=filename.split('/')[-1]
    a = filename.split('.PMS')[0]
    a = a.replace('_',' ').replace('R1','')
    s = ''.join(c for c in a if not c.isdigit())
    s = s.replace('.',' ').split()
    d = {}
    x = a[:]
    s.append(' ')
    for i in range(len(s)-1):
        if re.search(s[i],a):
            x = x.replace(s[i],'')
            y = x.split(s[i+1])[0]
            x = x.split(s[i+1])[-1]
            d[s[i]]=float(y)
    return d

def info_from_track_filename_env(file):
    info = info_from_track_filename(file)
    try:
        env = info['ENV']
    except KeyError:
        info['ENV'] = 0.0
    return info

def read_tracks(filename,loud=False):
    print 'Use get_tracks'
    sys.exit()
    return d

class Tracks(object):
    def __init__(self, data_array, col_keys,name):
        self.data_array = data_array
        self.key_dict = dict(zip(col_keys,range(len(col_keys))))
        self.name = name
        
    def get_row(self,i):
        return self.data_array[i,:]
        
    def get_model_row(self,i):
        row = np.nonzero(track.data_array[:,track.key_dict['MODE']] == i)[0]
        return self.data_array[row,:]
    
    def get_col(self,key):
        return self.data_array[:,self.key_dict[key]]
    
    def get_stable_col(self,key):
        lnuc = self.data_array[:,self.key_dict['LNUC']]
        ind = np.nonzero((lnuc < 1.2) & (lnuc > 0.8))[0]
        try:
            col = self.data_array[:,self.key_dict[key]]
            scol = col[ind]
        except KeyError:
            print key,'not found'
            scol = 0
        return scol

def get_tracks(filename,loud=False):
    if loud: 
        from time import time
        t1 = time()
    f=open(filename,'r')
    lines = f.readlines()
    f.close()
    
    for i in range(len(lines)):
        if lines[i].strip().startswith('BEGIN TRACK'):
            start_index = i
            break
    try:
        col_keys = lines[start_index+1].split()
    except UnboundLocalError:
        print 'No BEGIN TRACK ', filename
        return 0
    Ncols = len(col_keys)
    Nrows = len(lines[(start_index+2):])
    data = np.ndarray(shape=(Nrows,Ncols), dtype=float)
    row = 0
    for line in lines[(start_index+2):]:
        line = line.strip()
        try: # sometimes there's a sting...
            float(line[0])
        except ValueError:
            continue
        except IndexError:
            continue
        try:
            data[row] = map(float,line.split())
        except ValueError:
            try:
                data[row] = map(float,fortransucks(line))
            except ValueError:
                print 'This is beyond the fortan error... skipping: '
                print line
                print 'from ',filename
        row += 1
    TrackDict = Tracks(data,col_keys,filename)
    if loud == True: 
        print 'get_tracks:',Nrows,'tracks from',\
                os.path.split(filename)[1],'in',time()-t1,'s'
    return TrackDict
    
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
        tmp2.append(item)
    return tmp2

def OLDget_Hfus_HeB(t1,Filter1,Filter2):
    if t1 == 0: return 0
    try:
        mag1 = Tracks.get_col(t1,Filter1)
        mag2 = Tracks.get_col(t1,Filter2)
    except KeyError:
        mag1 = Tracks.get_col(t1,Filter1+'1')
        mag2 = Tracks.get_col(t1,Filter2+'1')
        
    QHE1 = Tracks.get_col(t1,'QHE1') # inner m/M of He fusion
    age = Tracks.get_col(t1,'AGE') 
    LY = Tracks.get_col(t1,'LY')
    LY = LY[1:]
    QHE1 = Tracks.get_stable_col(t1,'QHE1') # inner m/M of He fusion
    LX = Tracks.get_stable_col(t1,'LX') # inner m/M of He fusion
    QH1 = Tracks.get_stable_col(t1,'QH1') # inner m/M of H fusion
    Hfus = np.nonzero((QH1 == 0.) & (LX > 0.))[0] # L from H is > 0, in core, and not before star starts
    HeB = np.nonzero(LY > 0.)[0]

    return mag1,mag2,age,HeB,Hfus
    
def plot_tracks(tracks):
    # e.g.
    from GenUtils import get_afile
    tracks = [get_tracks(t) for t in get_afile('/Users/Phil/research/Italy//tracks/ell00/C_Z0.07_Y_0.389/','*phat')]
    T = [Tracks.get_stable_col(t,'LOG_TE') for t in tracks]
    L = [Tracks.get_stable_col(t,'LOG_L') for t in tracks]
    M = [info_from_track_filename(i.name)['M'] for i in tracks]
    x = [plt.plot(t,l,color='black') for t,l in zip(T,L)]
    x = [plt.annotate('%.2f'%m,xy=(max(t),l[argmax(t)]),color='blue',ha='right') for m,t,l in zip(M,T,L)]
    plt.axis([5,3.2,-2,4.5])
    
    '''
    ts = ['/Users/Phil/research/Italy//tracks/ell00/C_Z0.07_Y_0.389/Z0.07Y0.389OUTA1.72ENV_0.05_F7_M1.10.PMS.HB.phat','/Users/Phil/research/Italy//tracks/ell00/C_Z0.07_Y_0.389/Z0.07Y0.389OUTA1.72ENV_0.05_F7_M1.10.PMS.phat']
 
    Tt = [Tracks.get_stable_col(t,'LOG_TE') for t in tts]
    Ll = [Tracks.get_stable_col(t,'LOG_L') for t in tts]
    x = [plt.plot(t,l,color='red') for t,l in zip(Tt,Ll)]
    '''

def get_msto(track_obj):
    QH1 = Tracks.get_stable_col(track_obj,'QH1')
    model = Tracks.get_stable_col(track_obj,'MODE')
    msto = np.nonzero(QH1 > 0.0)[0][0]
    return model[msto]

def get_track_by_mass(tracks,mass):
    masses = [info_from_track_filename(track.name)['M'] for track in tracks]
    try:
        t = tracks[masses.index(mass)]
    except ValueError:
        print 'get_track_by_mass: mass = ',mass,'not found, returning all tracks'
        print masses
        t = tracks
    return t

def get_Hfus_HeB(t1,Filter1,Filter2):
    if t1 == 0: return 0
    try:
        mag1 = Tracks.get_col(t1,Filter1)
        mag2 = Tracks.get_col(t1,Filter2)
    except KeyError:
        mag1 = Tracks.get_col(t1,Filter1+'1')
        mag2 = Tracks.get_col(t1,Filter2+'1')
    lnuc = Tracks.get_col(t1,'LNUC')
    QHE1 = Tracks.get_col(t1,'QHE1') # inner m/M of He fusion
    age = Tracks.get_col(t1,'AGE') 
    LY = Tracks.get_col(t1,'LY')
    QHE1 = Tracks.get_col(t1,'QHE1') # inner m/M of He fusion
    LX = Tracks.get_col(t1,'LX') # inner m/M of He fusion
    QH1 = Tracks.get_col(t1,'QH1') # inner m/M of H fusion
    Hfus = np.nonzero((QH1 == 0.) & (LX > 0.8) & (lnuc < 1.1) & (lnuc > 0.9))[0] # L from H is > 0, in core, and not before star starts
    #HeB = np.nonzero((LY > 0.) & (lnuc < 1.2) & (lnuc > 0.8))[0]    
    #Hfus = np.nonzero((QH1 == 0.) & (LX > 0.))[0] # L from H is > 0, in core, and not before star starts
    HeB = np.nonzero(LY > 0.)[0]    
    return mag1,mag2,age,HeB,Hfus

def track2polygon(tmag1,tmag2):
    ''' 
    Take a cmd track and basically connects the ends to make it a polygon
    '''
    color = tmag1-tmag2
    verts = np.column_stack((color,tmag2))
    return verts

def spread_angst(trilegal_output,file_ast,outfile,filt1,filt2):
    '''
    correct trilegal output with artifical star tests. Usage:
    trilegal_output ast_file new_output_file, filt1, filt2
    adapted from Leo's spread_angst.c
    '''
    from GenUtils import tablefile2dict
    astfilt1,astfilt2 = file_ast.split('_')[-3],file_ast.split('_')[-2]
    if filt1 != astfilt1: 
        print 'AST filter1 (%s) doesn\'t match spread angst input (%s)'%(astfilt1,filt1)
    if filt2 != astfilt2: 
        print 'AST filter2 (%s) doesn\'t match spread angst input: (%s)'%(astfilt2,filt2)
        
    DELTAMAG = 0.2
    # Read AST
    t1 = time.time()
    ast_mag1,ast_mag2,ast_dmag1,ast_dmag2 = np.loadtxt(file_ast,unpack=True)
    ast_maxmag1 = np.max(ast_mag1)
    ast_maxmag2 = np.max(ast_mag2)
    nast = len(ast_mag1)
    # read trilegal out
    triout = read_table(trilegal_output)
    '''
    triout = np.loadtxt(trilegal_output,unpack=True)
    colheads = open(trilegal_output,'r').readline().replace('#','').split()
    
    try:
        mag1 = triout[colheads.index(filt1)]#+triout[colheads.index('m-M0')]
    except ValueError:
        print filt1,'not found in',trilegal_output
    try:
        mag2 = triout[colheads.index(filt2)]#+triout[colheads.index('m-M0')]
    except ValueError:    
        print filt2,'not found in',trilegal_output
    #obs = np.nonzero((mag1<ast_maxmag1) & (mag2<ast_maxmag2))[0]
    #random.shuffle(obs)    
    
    gc  = triout[colheads.index('Gc')]
    logage = triout[colheads.index('logAge')]
    mh = triout[colheads.index('[M/H]')]
    mini = triout[colheads.index('m_ini')]
    logl = triout[colheads.index('logL')]
    logte = triout[colheads.index('logTe')]
    logg = triout[colheads.index('logg')]
    dm0 = triout[colheads.index('m-M0')]
    av = triout[colheads.index('Av')]
    mr = triout[colheads.index('m2/m1')]
    try:
        mcore  = triout[colheads.index('Mcore')]
        co  = triout[colheads.index('C/O')]
        per = triout[colheads.index('Per')]
        pmode = triout[colheads.index('mode')]
        logmdot = triout[colheads.index('logML')]
        mact = triout[colheads.index('Mact')]
    except ValueError:
        'Error in trilegal_output file, was trilegal run with -a flag?'
        sys.exit()
    '''
    
    try:
        mag1 = triout.get_col(filt1)#+triout.get_col('m-M0')
    except ValueError:
        print filt1,'not found in',trilegal_output
    try:
        mag2 = triout.get_col(filt2)#+triout.get_col('m-M0')
    except ValueError:    
        print filt2,'not found in',trilegal_output
        
    #obs = np.nonzero((mag1<ast_maxmag1) & (mag2<ast_maxmag2))[0]
    #random.shuffle(obs)    
    
    try:
        mcore  = triout.get_col('Mcore')
        co  = triout.get_col('C/O')
        per = triout.get_col('Per')
        pmode = triout.get_col('mode')
        logmdot = triout.get_col('logML')
        mact = triout.get_col('Mact')
    except ValueError:
        'Error in trilegal_output file, was trilegal run with -a flag?'
        sys.exit()

    gc  = triout.get_col('Gc')
    logage = triout.get_col('logAge')
    mh = triout.get_col('[M/H]')
    mini = triout.get_col('m_ini')
    logl = triout.get_col('logL')
    logte = triout.get_col('logTe')
    logg = triout.get_col('logg')
    dm0 = triout.get_col('m-M0')
    av = triout.get_col('Av')
    mr = triout.get_col('m2/m1')

    
    fp_out = open(outfile,'w')
    print 'writing to',outfile
    # header
    fp_out.write('# gc logage mh mini logl logte logg dm0 av mr %s %s diff_%s diff_%s mcore co per pmode logmdot mact\n'%(filt1,filt2,filt1,filt2))
    
    # at least knock out the ones that don't get recovered at all...
    points = np.column_stack((mag1,mag2))
    verts = get_verts(ast_mag1,ast_mag2,dx=0.1,dy=0.1)
    mask = nxutils.points_inside_poly(points, verts)
    ind = np.nonzero(mask)[0]
    nind = np.nonzero(abs(mask-1))[0]
    mag1 = mag1[ind]
    mag2 = mag2[ind]
    
    # make nice little boxes of all the asts
    vertss = [np.column_stack(([ast_mag1[m]-DELTAMAG,ast_mag1[m]+DELTAMAG,ast_mag1[m]+DELTAMAG,ast_mag1[m]-DELTAMAG],[ast_mag2[m]+DELTAMAG,ast_mag2[m]+DELTAMAG,ast_mag2[m]-DELTAMAG,ast_mag2[m]-DELTAMAG])) for m in range(nast)]
    # shuffle the inds.
    ind = range(len(mag1))
    random.shuffle(ind)
    # combine each mag1 mag2 randomly
    pointss = [np.column_stack((mag1[i],mag2[i])) for i in ind]
    
    lost,found = [],[]
    for points,verts in zip(pointss,vertss):
        mask = nxutils.points_inside_poly(points, verts)
        found.extend(np.nonzero(mask)[0])
        lost.extend(np.nonzero(abs(mask-1))[0])
    
    t1 = time.time()
    for i in range(len(mag1)):
        # look for a suitable AST
        j = 0
        stop = 0
        ip = int(np.random.random()*float(nast)) # random position
        while stop == 0:
            j+=1
            if (j>=nast): 
                dmag1 = 99.
                dmag2 = 99.
                print 'dude, code better.'
                stop = 1
            if (ip>=nast): ip = 0
            if ( (mag1[i]<ast_mag1[ip]+DELTAMAG) & (mag1[i]>ast_mag1[ip]-DELTAMAG) & (mag2[i]<ast_mag2[ip]+DELTAMAG) & (mag2[i]>ast_mag2[ip]-DELTAMAG)):
                dmag1 = ast_dmag1[ip]
                dmag2 = ast_dmag2[ip]
                stop = 1
            ip+=1
        fp_out.write('%i %.2f %.2f %.5f %.3f %.3f %.3f %.2f %.3f %.2f %.3f %.3f %.3f %.3f %.2f %.2f %.1f %i %.2f %.3f\n'%(gc[i], logage[i], mh[i], mini[i], logl[i], logte[i], logg[i], dm0[i],av[i], mr[i], mag1[i], mag2[i], dmag1, dmag2,mcore[i], co[i], per[i], pmode[i], logmdot[i], mact[i]))
        if (ip>=nast): ip = 0
        if float(i)/len(mag1)*100 % 10 == 0: print float(i)/len(mag1)*100 
    t2 = time.time()
    fp_out.close()
    return

