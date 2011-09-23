"""
Functions useful in using ANGST pipeline output.
to use:
import sys
sys.path.append('/astro/users/philrose/python/')
sys.path.append('/Users/Phil/Dropbox/research/angstdir/python/')
from BinFitsFunctions import *

"""
import pyfits
import numpy as np
import re,os,subprocess
import sys
sys.path.append(os.environ['PYTHON_CODES'])
match_root = os.environ['MATCH_ROOT']
noisyCMD_path=os.path.join(match_root,'src/noisyCMD')
zcombine_path=os.path.join(match_root,'src/zcombine')
# calcsfh is run on wolverine...
calcsfh_path = '/astro/net/angst/projects/src/match2.3.1/match2.3/calcsfh'


# Examples
# Loadtxt:
# M1, M2, M, lZ, A = np.loadtxt(fcmd,dtype=float,unpack=True)

filename = '/Users/Phil/research/Italy/WFC3SNAP/Sp2011/fromJason/9755_IC2574-SGS_F435W_F555W_F814W_gst.param'
def parse_noisyCMD(filename):
    outfile = filename.split('.')[0]+'.dat'
    out = open(outfile,'w')
    npar = open(filename,'r').readlines()
    tmin,tmax,logz,sfr = np.transpose([map(float,line.split()) for line in npar[9:]])
    logzmin = logz-0.05
    logzmax = logz+0.05
    for i in range(len(tmax)):
        if tmax[i] == 10.15:
            tmax[i] = 10.13
    x=[out.write('%.3f %.3f %.7f %.7f %.12f 0.00000\n'%(logzmin[o],logzmax[o],tmin[o],tmax[o],sfr[o])) for o in range(len(tmin))]
    out.close()
    print 'parsed noisyCMD to',outfile
    return


def plot_processed_zctmp(processed_zctmps):
    sfhs = processed_zctmps
    for sfh in sfhs:
        zlow,zhigh,tmin,tmax,sfr,sfrerr = np.loadtxt(sfh,unpack=True)
        to = 10**(tmin-9)
        tf = 10**(tmax-9)
        agebins,sfrbins = lh_bins(to,tf,sfr)
        dump = plt.plot((to+tf)/2.,sfr,',',lw=3)
    column = 'SFR'
    data = read_table('/Users/Phil/research/Italy/WFC3SNAP/noAGB/zcmerged/SCL-DE1.final.mc')    
    col = Table.get_col(data,column)
    colep = Table.get_col(data,column+'err+')
    colem = Table.get_col(data,column+'err-')
    To = 10**(Table.get_col(data,'To')-9)
    Tf = 10**(Table.get_col(data,'Tf')-9)
    Agebins,colbins = lh_bins(To,Tf,col)

    plt.plot(Agebins,colbins,drawstyle='steps',color='black',lw=3)
    plt.errorbar((To+Tf)/2.,col,yerr=[colem,colep],elinewidth=1,linewidth=0,color='black')

    
def process_zctmp(extrastr='',zctmp=None,zc=None,processed=None):
    from GenUtils import get_afile
    if zctmp == None:
        zctmp = get_afile(os.getcwd()+'/','*'+extrastr+'*.zctmp')[0]
    if zc == None:
        zc =  get_afile(os.getcwd()+'/','*'+extrastr+'*.zc')[0]
    if processed == None:
        process = zc+'.dat'
    fed = 0
    p = open(process,'w')
    zct = open(zctmp,'r').readlines()
    metalbins = map(float,zct[1].strip().split())    
    halfbins = np.diff(metalbins)/2
    halfbin = halfbins[0] # this is just the equally spaced step in metallicity
    zc_data = read_zc(zc)
    sfr = zc_data['SFR']
    sfr_err = zc_data['SFRerr1']
    for j in range(2,len(zct)):
        #zct_line = time1 time2 sfr_z1 sfr_z2 ...
        zct_line = map(float,zct[j].strip().split())
        # sfr array is shifted two indices...
        if (sfr[j-2] != 0.): 
            fracerr = sfr_err[j-2]/sfr[j-2]
        else:
            fracerr = 0.
        l = 0
        for k in range(len(metalbins)):
            l+=1
            if zct_line[1] == 10.15: 
                zct_line[1]= 10.13
                fed = 1
            p.write('%.3f %.3f %.7f %.7f %.12f %.12f \n' % \
                            (metalbins[k]-halfbin, metalbins[k]+halfbin, \
                             zct_line[0], zct_line[1], \
                             zct_line[l+1], fracerr*zct_line[l+1]))
    
    p.close()
    print 'process_zctmp wrote',process
    if fed == 1: 'warning, found age of 10.15, changed to 10.13 to work with TRILEGAL 2.0'
    return process
    
def read_exclude_gates(filename):
    f = map(float,open(filename,'r').readline().strip().split()[1:-1])
    return np.column_stack((np.array(f[0::2]),np.array(f[1::2]-np.array(f[0::2]))))

def read_fitslist(filename=None):
    if filename==None:
        fits = np.genfromtxt('fits.list',autostrip=True,dtype="S")
    else:
        fits = np.genfromtxt(filename,autostrip=True,dtype="S")
    if len(fits) == 0:
        print 'Nothing read.'
    return fits                         

def this_fits(filename,fits=None):
    """
    Matches a filename PropID_Target_Filter1_Filter2* with a fits file
    (from ./fits.list) unless a fits (= np.ndarray of string values) is given
    """
    i=0
    thisfits = ''
    PropID,Target,Filter1,Filter2,filename = extract_title(filename)
    if fits == None:
        fits = read_fitslist()
    for fit in fits:
        if re.search(Target,fit):
            if re.search(Filter1,fit):
                thisfits = fit
                i+=1
    if i > 1:
        print i,' matches on ',filename
    return thisfits

def read_zc(filename):
    """
    reads a zcombine output file into a dictionary. This skips header lines in
    the file.
    """
    header,footer,bestfit,zcdata = [],[],[],[]
    file = open(filename,'r')
    Galaxy = extract_title(filename)[1]
    lines = file.readlines()    
    
    nlines = len(lines)
    for i in range(nlines):
        if lines[i].startswith('\n'): continue
        if lines[i].startswith('Found'):
            header.append(lines[i])
        elif lines[i].startswith('Best'):
            bestfit.append(lines[i+1])
            bestfit.append(lines[i+2])
        elif lines[i].startswith('background'):
            footer.append(lines[i])
        else:
            try:
                zcdata.append(map(float,lines[i].split()))
            except ValueError:
                continue
    zcdata = np.transpose(zcdata)
    data = {'To':zcdata[0],
            'Tf':zcdata[1],
            'Mag':zcdata[2],
            'SFR':zcdata[3],
            'SFRerr1':zcdata[4],
            'SFRerr2':zcdata[5],
            'Zave':zcdata[6],
            'Zaveerr1':zcdata[7],
            'Zaveerr2':zcdata[8],
            'Zspread':zcdata[9],
            'Zsprederr1':zcdata[10],
            'Zsprederr2':zcdata[11],
            'CSFH':zcdata[12],
            'CSFHerr1':zcdata[13],
            'CSFHerr2':zcdata[14],
            'Galaxy':Galaxy,
            'Bestfits':bestfit,
            'Header':header,
            'Footer':footer
          }
    return data

def run_zc(sfhfile,zcpars=None):
    """
    runs zcombine using the output star formation history file from calcsfh.
    zcpars is the zcombine parameter file containing time binning information.
    If none is specified, zcombine will be run with the same time binning as
    given in the match input parameter file. Other options include zc1 and
    zc2, these files will be built on the fly.
    zc1 will run zcombine in one bin, giving average sfh, z etc.
    zc2 will run zcombine in two bins, but stupidly. It just creates the
    zcpars file as if the time binning is from Phil's typical angst runs, it
    COULD read in the matchpars file, find where 1Gyr is and make it,
    I'll add that when I change my initial time bins...
    """
    pathzc = '/astro/net/angst/projects/src/match2.3.1/match2.3/src/zcombine'
    sfh = sfhfile
    msg = sfh.split('.sfh')[0] +'.msg'
    #cut = get_zcombcut(msg)
    zctmp = sfh +'.zctmp'
    zc = sfh +'.zc'
    if zcpars != None:
        print 'Running zcombine with ',zcpars,' paramater file'
    if zcpars == None:
        zcpars = 'foo'
        print 'Running zcombine in full time resolution'
        zc = sfh +'.fullres.zc'
        zctmp = sfh +'.fullres.zctmp'
    if zcpars == 'zc1':
        print 'Running zcombine with one time bin'
        zcpars = sfh + '.1zcpars'
        write_zcpars(zcpars)
        zc = sfh + '.1bin'
        zctmp = sfh +'.1bin.zctmp'
    if zcpars == 'zc2':
        print 'Running zcombine with two time bins'
        zcpars = sfh + '.2zcpars'
        write_zcpars(zcpars)
        zc = sfh + '.2bin'
        zctmp = sfh +'.2bin.zctmp'
    
    data = subprocess.Popen([pathzc,zcpars,' 9 1 ',zctmp,sfh,'-1'],
                            stdout=subprocess.PIPE).communicate()[0]
    f = open(zc,'w')
    f.writelines(data)
    f.close()
    print 'Wrote ',zc
    return zc
    
def write_zcpars(filename):
    """
    writes a zcombine parameter file for one bin or two ... for two bins
    it does it stupidly.
    It just creates the zcpars file as if the time binning is from Phil's
    typical angst runs,
    it COULD read in the matchpars file, find where 1Gyr is and make it,
    I'll add that when I
    change my initial time bins...

    Finds the one or two bins from the filename containing .1 or .2
    """
    match = get_file('.matchpars')[0]
    matchpars = read_matchpars(match)
    Ntbins = matchpars['Ntbins']
    lines = np.zeros(Ntbins+2,dtype=int)
    if filename.endswith('.1zcpars'):
        lines[0] = 1
    if filename.endswith('.2zcpars'):
        lines[0] = 2
        lines[25:] = 1
    np.savetxt(filename,np.transpose(lines),fmt='%i')

def get_zcombcut(msg):
    """
    gets a fit cut needed to run zcombine from the screen output of calcsfh
    basically a standard deviation returns the number as a string
    """
    data = read_matchmsg(msg)
    bestfit = data['BestFit']
    fits = data['fit']
    nfits = float(len(fits))
    cut = sum(fits-bestfit)/nfits + bestfit
    return str(cut)

def read_matchmsg(filename):
    """
    reads the screen output of calcsfh into a dictionary. If calcsfh was run
    with a -zinc flag there will be three more columns, currently this isn't
    supported. This is because Phil is lazy.
    """
    file = open(filename,'r')
    lines = file.readlines()
    Nstars = int(lines[5].split()[0])
    Nfakes = int(lines[6].split()[0])
    dataline=[]
    for line in lines:
        if re.match('Av',line): # like starts with
            tmp = line.replace('=',',')
            tmp = tmp.replace(':',',')
            tmp = tmp.replace('\n','')
            tmp = tmp.split(',')
            dataline.append(map(float,tmp[1:8:2]))
        if re.match('Best',line):
            bestline = line.replace('=',',')
            bestline = bestline.replace('\n','')
            bestline = bestline.split(',')
    msgdata = np.transpose(dataline)
    data = { 'Av': msgdata[0],
             'imf': msgdata[1],
             'dmod': msgdata[2],
             'fit': msgdata[3],
             'BestAv': float(bestline[1]),
             'BestDmod': float(bestline[3]),
             'BestFit': float(bestline[5])
        }
    return data

def read_matchpars(filename):
    """
    reads calcsfh parameter file into a dictionary. Doesn't work for -zinc
    flag because Phil is lazy.
    """
    file = open(filename,'r')
    IMF,dmodmin,dmodmax,dmodstep,Avmin,Avmax,Avstep = file.readline().split()
    logZmin,logZmax,dlogZ= file.readline().split()
    BF,Bad0,Bad1= file.readline().split()
    Ncmds= file.readline().split()
    Mag1step,Colorstep,fake_sm,Colormin,Colormax,Colors=file.readline().split()
    Mag1min, Mag1max, Mag1name = file.readline().split()
    Mag2min, Mag2max, Mag2name  = file.readline().split()
    Ntbins= file.readline()
    lines = file.readlines()
    bgline2 = lines.pop()
    bgline1 = lines.pop()
    times = []
    for line in lines:
        times.append(map(float,line.split()))
    
    tbins = np.transpose(times)
    data = {'IMF':float(IMF),
            'dmodmin': float(dmodmin), 
            'dmodmax': float(dmodmax), 
            'dmodstep': float(dmodstep), 
            'Avmin': float(Avmin), 
            'Avmax': float(Avmax), 
            'Avstep': float(Avstep), 
            'logZmin': float(logZmin), 
            'logZmax': float(logZmax), 
            'dlogZ': float(dlogZ), 
            'BF': float(BF), 
            'Bad0': float(Bad0), 
            'Bad1': float(Bad1), 
            'Mag1step': float(Mag1step), 
            'Colortep': float(Colorstep), 
            'fake_sm': float(fake_sm), 
            'Colormin': float(Colormin), 
            'Colormax': float(Colormax), 
            'Colors': Colors, 
            'Mag1min':float( Mag1min), 
            'Mag1max': float(Mag1max), 
            'Mag1name': Mag1name, 
            'Mag2min': float(Mag2min), 
            'Mag2max': float(Mag2max), 
            'Mag2name': Mag2name, 
            'Ntbins': int(Ntbins), 
            'To': tbins[0], 
            'Tf': tbins[1], 
            'bgline2': bgline2, 
            'bgline1': bgline1 
        }
    return data

def readbintab_plus(file):
    fits = pyfits.open(file)
    data = fits[1].data
    camera = fits[0].header['CAMERA']
    Mag1name = 'MAG1_'+camera
    Mag2name = 'MAG2_'+camera
    Mag1 = data.field(Mag1name)
    Mag2 = data.field(Mag2name)
    ra = data.field('RA')
    dec = data.field('DEC')
    return Mag1,Mag2,ra,dec,data

def readbintab(file):
    names=extract_title(file)
    fits = pyfits.open(file)
    data = fits[1].data
    camera = fits[0].header['CAMERA']
    if re.search('WFC3',camera): camera = camera.replace('WFC3-','')
    Mag1name = 'MAG1_'+camera
    Mag2name = 'MAG2_'+camera
    Mag1 = data.field(Mag1name)
    Mag2 = data.field(Mag2name)
    ra = data.field('RA')
    dec = data.field('DEC')
    bintab = {'Mag1': Mag1,
              'Mag2': Mag2,
              'ra': ra,
              'dec': dec,
              'names': names}
    return bintab

def readbintab_fake(file):
    fits = pyfits.open(file)
    names=extract_title(file)
    data = fits[1].data
    Mag1in = data.field('MAG1IN')
    Mag2in = data.field('MAG2IN')
    Mag1out = data.field('MAG1OUT')
    Mag2out = data.field('MAG2OUT')
    ra = data.field('RA')
    dec = data.field('DEC')
    Mag1diff = Mag1out - Mag1in
    Mag2diff = Mag2out - Mag2in
    bintab = {'Mag1in': Mag1in,
              'Mag2in': Mag2in,
              'Mag1diff': Mag1diff,
              'Mag2diff': Mag2diff,
              'ra': ra,
              'dec':dec,
              'names': names}
    return bintab

def extract_title(file):
    tmp = file.split('/')
    filename = tmp[-1]
    split = filename.split('_')
    PropID = split[0]
    Target = split[1]
    Filter1 = split[2]
    Filter2tmp = split[3]
    tmp = Filter2tmp.split('.')
    Filter2 = tmp[0]
    return PropID,Target,Filter1,Filter2,filename

# read_fake
# just do this:
# mag1, mag2, mass, age, logz = np.loadtxt(filename,unpack=True)

def read_zc_fake(filename):
    file = open(filename,'r')
    lines = file.readlines()
    
    nfitstmp = lines[1].split()
    nfits = int(nfitstmp[1])
    
    bests = lines[4].split(',')
    Av,AvErr = zc_strings(bests[0])
    IMF,IMFErr = zc_strings(bests[1])
    dmod, dmodErr = zc_strings(bests[2])
    nlines = len(lines)
    rows = lines[6:nlines-2]
    
    to,tf,mag,SFR,SFRerr1,SFRerr2= [],[],[],[],[],[]
    Zave,Zaveerr1,Zaveerr2,zspread,zsprederr1=[],[],[],[],[]
    zsprederr2,cSFH,cSFHerr1,cSFHerr2= [],[],[],[]
    for row in rows:
        data = map(str,row.split())
        to.append(data[0])
        tf.append(data[1])
        mag.append(data[2])
        SFR.append(data[3])
        SFRerr1.append(data[4])
        SFRerr2.append(data[5])
        Zave.append(data[6])
        Zaveerr1.append(data[7])
        Zaveerr2.append(data[8])
        zspread.append(data[9])
        zsprederr1.append(data[10])
        zsprederr2.append(data[11])
        cSFH.append(data[12])
        cSFHerr1.append(data[13])
        cSFHerr2.append(data[14])
    
    return to,tf,dmod,Av,SFR,IMF,zspread,Zave

def zc_strings(somestring):
    # ex: 'Av=0.150+/-0.042,'
    # or this:  'Av=0.141+0.059-0.041'
    # sometimes you get this:
    #' IMF=1.350+-0.000-0.000'
    # this will return 1.350,['', '0.000', '0.000']
    # which is awesome.
    q = somestring.split('+/-')
    if q[0] == somestring:
        q = somestring.split('+')
        q1 = q[1:]
        err = q1[0].split('-')
    else:
        err = q[1]
    r = q[0].split('=')
    datum = r[1]
    return datum,err

def read_matchpars_fake(filename):
    file = open(filename,'r')
    lines = file.readlines()
    
    filttmp = lines[4].split()
    filters = filttmp[-1]
    Vmaxtmp = lines[5].split()
    Vmax = Vmaxtmp[1]
    Imaxtmp = lines[6].split()
    Imax = Imaxtmp[1]
    return filters,Vmax,Imax

def get_trgb(filename):
    from LatexUtils import ReadLatexTable
    PropID,Target,Filter1,Filter2,filename = extract_title(filename)
    #tab5 = ReadLatexTable('/astro/users/philrose/python/tab5.tex')
    tab5 = ReadLatexTable('/Users/Phil/Dropbox/research/python/tables/tab5.tex')
    for i in range(len(tab5['Target Name'])):
        if Target == tab5['Target Name'][i]:
            trgb = tab5['m_TRGB'][i]
    return trgb

def write_match(filename,Mag1,Mag2):
    match = Mag1,Mag2
    np.savetxt(filename,np.transpose(match),fmt=' %7.4f %7.4f',delimiter=' ')
    print 'Wrote '+filename

def write_matchfake(filename,Mag1in,Mag2in,Mag1diff,Mag2diff):
    matchfake = Mag1in,Mag2in,Mag1diff,Mag2diff
    np.savetxt(filename,
               np.transpose(matchfake),
               fmt=' %7.4f %7.4f %7.4f %7.4f',
               delimiter=' ')
    print 'Wrote '+filename

def DirsToGalaxy(dirs):
	"""
	Takes either a list of PropID_Target or a string PropID_Target
	and tries to make it match the ANGST paper names (eg UGC => U)

        This is really just the calling function. See GalStrings for crazy
        shit I have to do to have a consistent galaxy name.
	"""
	if type(dirs) is list:
            gals = []
            for dir in dirs:
                gals.append(GalStrings(dir))
        if type(dirs) is str:
            gals = GalStrings(dirs)
        return gals

def GalStrings(gal):
	gal = gal.split('_')[1]
	gal = gal.replace('GC','')
	gal = gal.replace('SO','')
	ngal = gal
	# Now some stuff I had to add for the different catalogue names...
        # only works on a case by case basis,
        # I had to add them when I realized the pipeline dir
	# names don't always match the primary names from catalogues.
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
	return ngal

def write_qsub(param,phot,fake,qsubfile,zinc=True,mc=False,cwd=None):
    flags = ''
    if zinc == True: flags = '-zinc'
    qsub = open(qsubfile,'w')
    if cwd == None: cwd = os.getcwd()
    fits = phot.split('/')[-1].split('.match')[0]
    sfh = fits+'.sfh'
    log = fits+'.log'
    msg = fits+'.msg'
    if mc==True:
        if not re.search('mc',cwd):
            log = 'mc/'+log
            msg = 'mc/'+msg
            sfh = 'mc/'+sfh
    lines = []
    lines.append('#PBS -l nodes=1:ppn=1 \n')
    lines.append('#PBS -j oe \n')
    lines.append('#PBS -o '+ cwd+'/'+log + '\n')
    lines.append('#PBS -l walltime=12:00:00 \n')
    if mc == False:
        lines.append('#PBS -M philrose@astro.washington.edu \n')
    
    lines.append('#PBS -m abe \n')
    lines.append('#PBS -V \n')
    lines.append('cd '+cwd+'\n')
    if mc == False:
        lines.append('%s %s %s %s %s %s > %s \n' %(calcsfh_path,param,phot,fake,sfh,flags,msg))
    else:
        lines.append('%s %s %s %s %s -allstars -logterrsig=0.03 -mbolerrsig=0.41 %s > %s \n' %(calcsfh_path,param,phot,fake,sfh,flags,msg))     
    qsub.writelines(lines)

def read_zctmp(filename):
    f = open(filename,'r')
    ncol,nrow = map(int,f.readline().split())
    data = {'To':[],'Tf':[],'sfr':[],'logz':[]}
    data['logz'] = map(float,f.readline().split())
    for line in f:
        to = float(line.split()[0])
        tf = float(line.split()[1])
        sfr = map(float,line.split()[2:])
        data['To'].append(to)
        data['Tf'].append(tf)
        data['sfr'].append(sfr)
    f.close()
    return data
