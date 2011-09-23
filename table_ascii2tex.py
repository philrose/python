#
#  table_ascii2tex.py
#  
#
#  Created by Philip Rosenfield on 11/24/10.
#
import sys
import numpy as np
import os

def ascii2tex(file):
    table = open(file,'r').readlines()
    tex = file.replace(os.path.splitext(file)[1],'.tex')
    out = open(tex,'w')
    print 'Writing to ',tex
    ncols=0
    for line in table:
        if line.startswith('#'): continue
        n = len(line.split())
        if n > ncols: ncols = n

    for line in table:
        if line.startswith('#'): 
            out.write( '%s \\\ \n' %line.strip().replace('#','%'))
            continue
        data = line.strip().split()
        data_line = ' & '.join(data)
        ndata = len(data)
        missing = ncols-ndata+1
        missing_vals ='&'.join(missing*' ')
        out.write('%s %s \\\ \n' % (data_line,missing_vals))

    out.close()
    return tex

if __name__ == '__main__':
    file =ascii2tex(sys.argv[1])
