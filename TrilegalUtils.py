#
#  TrilegalUtils.py
#  

import matplotlib.pyplot as plt
import numpy as np
from GenUtils import is_numeric,check_mkdir,Table,read_table
import time
def color_color_diagnostic(trilegal_output_file,filter1,filter2,filter3,filter4,ax=None,**pltkwargs):
    print 'reading',trilegal_output_file
    tstart = time.time()
    t = read_table(trilegal_output_file)
    tend = time.time()
    print 'reading took %.2f seconds'%(tend-tstart)
    #mag1 = np.array(t[filter1])
    #mag2 = np.array(t[filter2])
    #mag3 = np.array(t[filter3])
    #mag4 = np.array(t[filter4])
    
    mag1 = Table.get_col(t,filter1)
    mag2 = Table.get_col(t,filter2)
    mag3 = Table.get_col(t,filter3)
    mag4 = Table.get_col(t,filter4)
    
    if ax == None: ax = plt.axes()
    
    ax.plot(np.array(mag1)-np.array(mag2),np.array(mag3)-np.array(mag4),',',label=trilegal_output_file,**pltkwargs)
    ax.set_ylabel('$%s-%s$'%(filter1,filter2))
    ax.set_xlabel('$%s-%s$'%(filter3,filter4))
    return ax
      
  
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
    print 'dont use get_data, use GenUtils.read_table'
    sys.exit()
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

#me
#mc_dir = '/Users/Phil/research/Italy/WFC3SNAP/noAGB/mcTests/10503_SCL-DE1/mc/'

#Jason
mc_dir = '/Users/Phil/research/Italy/WFC3SNAP/JASON/sfh_overlap/10503_SCL-DE1'

def mc_tests(ID,dir_name,outdir="MC_TESTS",bestfit_loc = "MODELS/"):
    '''
    ID is the Galaxy name
    dir_name is PropID_GalaxyName
    Will put mc tests in MC_TESTS/*/dirname
    '''
    
    # Make new folders
    folders = ['PARS/','INPUT/','OUTPUT/']
    outdir = check_mkdir(outdir)
    dummy = [check_mkdir(os.path.join(outdir,folder)) for folder in folders]
    this_outdir = check_mkdir(os.path.join(outdir,folder,dir_name))
    
    # A place to write the commands if this needs to be run again
    cmds_out = check_mkdir(os.path.join(outdir,'MC_COMMANDS/'))
    cmds_file = open('%s/%s_MC_commands.dat'%(cmds_out,ID),'r')
    
    # Place for the output table
    table_out = check_mkdir(os.path.join(outdir,'TABLES/'))
    out = open('%s/%s_MC.dat'%(table_out,ID),'w')
    out.write('# mc_ID p_value NRGB_data NAGB_data NRGB_model NAGB_model mass_model N_wind Flux1_wind Flux2_wind\n')

    #best fits: bfs[0] = pars for run_trilegal, bfs[1] = input for trilegal bfs[2] = output for trilegal
    bfs = [get_afile(bestfit_loc+folder,'*'+ID+'*')[0] for folder in  folders]
    
    # find model from file!!! This is not general! ! 
    model = bfs[-1].split('model_')[-1].replace('.dat.dat','.dat')
    
    # Switch to the MC test directories
    new_place = [bf.replace(bestfit_loc,'MC_TESTS/').replace(folder,folder+dir_name+'/') for bf in bfs]
    
    # Best fit pars is the template, we'll swap out sfh file
    in_pars = open(bfs[0]).readlines()
    
    # Load SFHs
    sfhs = get_afile(mc_dir+'/','SFR*mc*dat')
    
    for i,sfh in enumerate(sfhs):
        mcid = sfh.split('/')[-1].split('.')[1] # need to change jason/me
        new_names = [np.replace(ext,'.'+mcid+ext) for np,ext in zip(new_place,['.pars','.dat','.dat.dat'])]
        # New pars file for run_trilgal.py
        pfile=open(new_names[0],'w')
        [pfile.write(inpar) for inpar in in_pars[:-1]] #sfh line is at the bottom.
        pfile.write("%-18s %s\n"%('object_sfr  ',sfh)) 
        pfile.close()
        
        cmd="/Users/Phil/research/Italy/WXTRILEGAL/run_trilegal.py "
        #cmd="/home/philrose/WXTRILEGAL/run_trilegal.py "
        cmd+="-e code_2.0/main "
        cmd+="-a "
        cmd+="-i %s "%new_names[1]
        cmd+="-o %s "%new_names[2]
        cmd+="-f %s "%model
        cmd+=new_names[0]
        
        cmds_file.write('%s \n'%cmd)
        print 'running TRILEGAL:',model,ID
        p = subprocess.Popen(cmd, shell=True,stdout=subprocess.PIPE)
        p = subprocess.Popen(cmd, shell=True)
        sts = os.waitpid(p.pid, 0)[1]
        
        fak_file=get_fakFILE(ID,jason=jason)
        ast_file = new_names[2].split('/')[0]+'/ast_'+new_names[2].split('/')[-1]
        cmd="AST/code/spread_angst <<EOF \n"
        cmd+=fak_file+"\n"
        cmd+=new_names[2]+"\n"
        cmd+=ast_file+"\n"
        cmd+="EOF \n"
        print "  ... completeness using %s"%fak_file
        print "  %s -> %s"%(new_names[2],ast_file)
        print 'Running spread_angst...'
        p = subprocess.Popen(cmd, shell=True,stdout=subprocess.PIPE)
        sts = os.waitpid(p.pid, 0)[1]
        os.system("wc -l %s %s|head -2"%(new_names[2],ast_file))
        
        cmds_file.write('%s \n'%cmd)
        
        synthcmd = read_table(ast_file)
        s_mag2 = synthcmd.get_col('mag2') + synthcmd.get_col('diff_mag2'.strip())
        s_mag2 = s_mag2[np.nonzero(abs(synthcmd.get_col('diff_mag2'.strip())) < 90.)[0]]
        s_mag1 = synthcmd.get_col('mag1') + synthcmd.get_col('diff_mag1'.strip())
        s_mag1 = s_mag1[np.nonzero(abs(synthcmd.get_col('diff_mag1'.strip())) < 90.)[0]]
        s_color = s_mag1-s_mag2
        Norm = trgb + 1.5
        ind,nB_AGB,nNorm,ps_nNorm,ps_nB_AGBm,hist,bins,s_hist_normed,p_value,normalization = calc_LF(mag2,s_mag2,Norm,trgb)
        Nstars, flux_rates = flux_from_mass_loss(synthcmd,rates,[filt1,filt2],ast_inds=ind,rel_flux=True)
        out.write('# mc_ID p_value NRGB_data NAGB_data NRGB_model NAGB_model mass_model N_wind Flux1_wind Flux2_wind\n')
        out.write('%s %.3f %i %i %i %i %e %i %e %e \n' %(mcid,p_value,nNorm,nB_AGB,ps_nNorm,ps_nB_AGBm,object_mass,Nstars[0],flux_rates[0][1],flux_rates[1][0]))
    
        os.remove(ast_file)
        print 'deleted',ast_file
        os.remove(new_names[2])
        print 'deleted',new_names[2]

    out.close()
    cmds_file.close()
