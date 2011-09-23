import sys,os,re,shutil,subprocess
import numpy as np


def fix_mass(file):
    infile = open(file,'r')
    lines = infile.readlines()
    infile.close()
    mass1 = get_mass(lines[0])
    mass2 = get_mass(lines[1])
    mass3 = float(lines[3].strip().split()[1])
    if mass1 != mass2: 
        print 'Something is fucked dude.'
        return
    if mass3 != mass2:
        lines[3] = lines[3].replace(str(mass3),str(mass2))
        shutil.copy(file,file+'.bak')
        outfile = open(file,'w')
        outfile.writelines(lines)
        outfile.close()
        print 'Fixed mass of ',file
    return

def get_crit_pt(dataline):
        critpt =dataline.strip().replace('\'','').split()[1]
        if critpt.startswith('C='): 
            critpt = 'C='
        return critpt

def fix_crit_pts(coms,num):
    # take a the critical points out of isofile that aren't 
    # found in lcom write transition iso to outiso
    # but if there is a transition from C= to F need to copy
    # the C= file and change the C= to an F
    if num == -1: extra = '.b'
    if num == 1: extra = '.f'
    for j in range(len(coms)-1):
        if( num == -1) and (j == 0): continue
        othercom = open(coms[j+num],'r').readlines()
        com = open(coms[j],'r').readlines()
        header = com[:4]
        cpdata = header
        data = com[4:]
        otherdata = othercom[4:]

    # get critical points
        critpts = []
        for i in range(len(data)-1):
            critpts.append(get_crit_pt(data[i]))
            
        other_critpts = []
        for i in range(len(otherdata)-1):
            other_critpts.append(get_crit_pt(otherdata[i]))
    
    # find missing critical points forward
        missing_crits = [a for a in critpts if not a in other_critpts]
        f_check = [a for a in other_critpts if not a in critpts]
        #print missing_crits,'f',f_check
        if len(missing_crits) != 0:
            if missing_crits[0] == 'F': continue
            print coms[j],' has ',missing_crits,' and ',coms[j+num],' does not.'
            isofile = coms[j].replace('com','iso')
            # write the lines in isofile that aren't tagged with 
            # missing critical points to outiso
            iso = open(isofile,'r').readlines()
            outiso = isofile.replace('.iso',extra+'tran.iso')
            f = open(outiso,'w')
            for line in iso:
                cut = 0
                for critpt in missing_crits:
                    if re.search(critpt,line.strip()): 
                        cut = 1
                        crit_full = line.strip().split(' ')[5]
                        if critpt == 'C=':
                            print 'exchanging a C= for an',f_check[0]
                            f.write(line.replace(crit_full,'F       '))
                        else:
                            print 'taking ',critpt,' out.'
                            f.write(line.replace(crit_full,'     '))
                if cut == 0:
                    f.write(line)                          
            f.close()
            print 'wrote ',outiso
    return

def sort_by_mass(unsort):
    # unsort is a filename *.PMS...
    masses,sbm = [],[]
    for file in unsort:
        masses.append(get_mass(file))
    for i in range(len(unsort)):
        sbm.append(unsort[np.argsort(masses)[i]])    
    return sbm

def get_mass(filestring):
    return float(filestring.split('_M')[-1].split('.PMS')[0])

def get_param_from_PMS(filename,paramstr):
    filename=filename.split('/')[-1]
    a = filename.split('.PMS')[0]
    s = ''.join(c for c in a if not c.isdigit())
    s = s.replace('.',' ').split()
    a = a.replace('_','')
    Z = a.split(s[1])[0].replace(s[0],'')
    Y = a.split(Z)[-1].split(s[2][0])[0].replace(s[1],'')
    if paramstr == 'Z': 
        return float(Z)
    elif paramstr == 'Y': 
        return float(Y)
    else:
        print paramstr,' not found in ',prefix
        return 'NaN'

def write_isotrack(filename,isochrone,Z,Y):
    cmd_input = open(filename,'w')
    cmd_input.write('1 \n')
    cmd_input.write('%.3f %.3f \n' % (Z,Y))
    cmd_input.write('isotrack/tutto/low_z0001_corr.dat  \n')
    cmd_input.write('isotrack/tutto/hb_z0001.dat \n')
    cmd_input.write('isochrones/%s \n' % 
                    (isochrone))
    cmd_input.close()
    print 'isotrack input file written: ',filename

def write_sfr_const(filename,Z):
    sfr_const = open(filename,'w')
    ages = [1e+06,1.5e+06,2e+06,3e+06,4e+06,6e+06,1e+07,1.5e+07,
            2e+07,3e+07,4e+07,6e+07,1e+08,1.259e+08,1.585e+08,
            1.995e+08,2.512e+08,3.162e+08,3.981e+08,5.012e+08,
            6.31e+08,7.943e+08,1e+09,1.259e+09,1.585e+09,
            1.995e+09,2.512e+09,3.162e+09,3.981e+09,5.012e+09,
            6.31e+09,7.943e+09,1e+10,1.259e+10,1.26e+10]
    sfr = 1.0000
    for age in ages:
        sfr_const.write('%.2e %.4f %.3f \n' % (age,sfr,Z))
    sfr_const.close()
    print 'sfr const written, ',filename,' put in tab_sfr/.'

def write_cmd_input(filename,isotrack):
    shutil.copy('/home/philrose/cmd_input_pr.dat',filename)
    lines =open('/home/philrose/cmd_input_pr.dat','r').readlines()
    out = open(filename,'w')
    tmp = lines[1].split(' ')
    tmp[1] = 'isotrack/'+isotrack
    lines[1] = ' '.join(tmp)
    out.writelines(lines)
    print 'cmd input written: ',filename

def write_trilegal_input(filename,tab_sfr):
    shutil.copy('/home/philrose/input_const_phil',filename)
    lines =open('/home/philrose/input_const_phil','r').readlines()
    out = open(filename,'w')
    tmp = lines[-3].split(' ')
    tmp[0] = 'tab_sfr/'+tab_sfr
    lines[-3] = ' '.join(tmp)
    out.writelines(lines)
    print 'trilegal input written, ',filename

def write_trilegal_output(filename):
    out = open(filename,'w')
    out.write(' \n')
    out.close()
    print 'wrote empty file for trilegal output, ',filename


class trilegal(object):
    def __init__(self,data_array,col_keys):
        self.data_array = data_array
        self.key_dict = dict(zip(col_keys,range(len(col_keys))))

    def get_row(self,i):
        return self.data_array[i,:]
    
    def get_col(self,key):
        return self.data_array[:,self.key_dict[key]]
    
def get_trilegal_out(filename):
    lines = open(filename,'r').readlines()
    col_heads = lines[0].strip().replace('#','')
    col_keys = col_heads.split()
    nrows = len(lines)
    ncols = len(col_keys)
    print 'trilegal output read. nrows:',nrows,'ncols:',ncols
    data = np.ndarray(shape=(nrows,ncols),dtype=float)    
    row = 0
    for line in lines:
        if line.startswith('#'): continue
        line = line.strip()
        data[row] = map(float,line.split())
        row +=1

    return trilegal(data,col_keys)
