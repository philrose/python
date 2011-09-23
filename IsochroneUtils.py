import numpy as np
import sys,os
sys.path.append(os.environ['PYTHON_CODES'])
from GenUtils import is_numeric

class Isochrone(object):
	def __init__(self,
					metallicity, #metallicity
					age,        #age
					data_array, #slice of data array returned by create_data_array
					col_keys,   #dictionary returned by get_col_keys
				):
		self.metallicity = metallicity
		self.age = age
		self.data_array = data_array
		self.key_dict = dict(zip(col_keys,range(len(col_keys))))

	def get_row(self,i):
		return self.data_array[i,:]

	def get_col(self,key):
		return self.data_array[:,self.key_dict[key]]
	
def create_key(metallicity,age):
	#create a unique dictionary key for each isochrone
	return 'ISO_%.2g_%.2g' % (metallicity,age) #(or whatever makes sense)

def get_all_isochrones(filename):
    #first go through the file and find how many total rows and columns
    #of data you have, how many isochrones you have, and the starting
    #index & metadata of each isochrone
    #  pretend these are the values you get
    f=open(filename,'r')
    lines = f.readlines()
    f.close()

    start_indices = []
    metallicities = []
    ages = []
    Nrows = 0
    N_isochrones = 0
    for i in range(len(lines)):
        # The start of the column head # log age/yr
        if lines[i].startswith('# l'):
            start_indices.append(i+1)
        if not lines[i].startswith('#'):
            Nrows += 1
        if lines[i].startswith('#\tI'):
            N_isochrones += 1
            line = lines[i]
            line = line.split()
            metallicities.append(float(line[4]))
            ages.append(float(line[7]))
            
    colhead = lines[start_indices[0]-1].strip()
    colhead = colhead.replace('#','')
    colhead = colhead.replace('/','')
    col_keys = colhead.split()
    col_keys[0] = 'LogAge'

    Ncols = len(col_keys)

    #now go back through the file and read in all the rows and columns of data
    data = np.ndarray(shape=(Nrows,Ncols), dtype=float)
    row = 0
    for i in range(len(lines)):
        #print i
        if lines[i].startswith('#'): continue
        items = [is_numeric(m) for m in lines[i].split()]
        if len(lines[i].split()) != Ncols: items.append('')
        data[row] = items
        row += 1

    IsoDict = {}

    #this will help with slicing: see below
    start_indices = np.concatenate((start_indices,[-1]))

    #create your dictionary of isochrone objects
    for i in range(N_isochrones):
        key = create_key(metallicities[i],ages[i])
        IsoDict[key] = Isochrone(metallicities[i],
                                    ages[i],
                                    data[start_indices[i]:start_indices[i+1]],
                                    col_keys)
    return IsoDict, col_keys