#
#  scrap.py
#  
#
#  Created by Philip Rosenfield on 9/28/10.
#  Copyright (c) 2010 __MyCompanyName__. All rights reserved.
#
class DataTable(object):
	def __init__(self,
					data_array, #slice of data array returned by create_data_array
					col_keys,   #dictionary returned by get_col_keys
				):
		self.data_array = data_array
		self.key_dict = dict(zip(col_keys,range(len(col_keys))))

	def get_row(self,i):
		return self.data_array[i,:]

	def get_col(self,key):
		return self.data_array[:,self.key_dict[key]]

def load_data_table(filename,comment_char):
    f=open(filename,'r')
    lines = f.readlines()
    f.close()
    
    for i in range(len(lines)):
        if lines[i].startswith(comment_char):
            line = lines[i].replace(comment_char,'')
            col_keys = line.split()
        break
    
    Ncols = len(col_keys)
    Nrows = len(lines[i+1:])
    data = np.ndarray(shape=(Nrows,Ncols))
    for line in lines[i+1:]:
        row = line.split()
        for j in range(Ncols):
            data[j] = is_numeric(row[j])
    
    Data = DataTable(data,col_keys)
    return Data