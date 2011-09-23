import pyfits
import numpy as np
import re,os,subprocess,sys
sys.path.append('/astro/users/philrose/python/')
from GenUtils import *

def write_latex_table(dict,*keyorder,**kw):
    # if python 3.0 could use filename='' etc... oh well.
    # kw are below, none are optional, none have default values.
    # ex:
    # write_latex_table(dict,'Catalog Name','Filters', filename=None,
    #                  colfmt=None,tablecaption=None,tablelabel=None,
    #                  commentsFile=None)
    if kw['filename'] == None:
        filename = 'Table.tex'
    else:
        filename = kw['filename']
    if kw['colfmt'] == None:
        colfmt='{}'
    if kw['tablecaption'] == None:    
        tablecaption='{}'
    if kw['tablelabel'] == None:
        tablelabel='{}'
    if kw['commentsFile'] == None:
        commentsFile=None
    f = open(filename,'w')
    f.write('\\begin{deluxetable}'+colfmt+'\n')
    f.write('\\setlength{\\tabcolsep}{0.02in}\n')
    f.write('\\tabletypesize{\\scriptsize}\n')
    f.write('\\tablecaption{'+tablecaption+'}\n')
    f.write('\\tablehead{ ' )
    dic = DictToLatex(dict)
    colhead,scolhead = LatexTableHeader(keyorder)
    for col in colhead:
        f.write(col)
    for scol in scolhead:
        f.write(scol)
    f.write('}\n')
    f.write('\startdata \n')
    nrows = len(dic[dic.keys()[0]])
    nkeys = len(keyorder)
    j = 0
    for i in range(nrows):
        for key in keyorder:
            j += 1
            if j == nkeys:
                f.write(str(dic[key][i])+' ')
            else:
                f.write(str(dic[key][i])+' & ')
        f.write('\\\\ \n')
        j = 0
    f.write('\enddata \n')
    if commentsFile != None:
        f.write(LatexTableComments(commentsFile))
    f.write('\\label'+tablelabel+' \n')
    f.write('\end{deluxetable}')


def DictToLatex(dic):
    # take err and split ... and make $##\pm###$ 
    for key in dic.keys():
        if re.search('err',key):
            errkey = key
            measkey = errkey.replace(' err',"")
            for i in range(len(dic[measkey])):
                dic[measkey][i] = '$'+str(dic[measkey][i])+'\\pm'+str(dic[errkey][i])+'$'
            del dic[errkey]
    return dic

def LatexTableHeader(keys):
    nkeys = len(keys)
    colhead = []
    scolhead = []
    for key in keys:
        # fold at first space to \\ \colhead{}
        tkey = key.split(' ')
        if len(tkey) > 1:
            col = '\\colhead{'+LatexSuperSubScript(tkey[0])+'} & \n'
            scol = '\\colhead{'+ LatexSuperSubScript(tkey[1])+'} & \n'
            fold = 1
        else:
            col = '\\colhead{'+LatexSuperSubScript(tkey[0])+'} & \n'
            scol = '\\colhead{} & \n'
        colhead.append(col)
        scolhead.append(scol)
    if fold == 1:
        colhead[-1] = colhead[-1].replace('&','\\\\')
    return colhead,scolhead

def LatexSuperSubScript(key):
    # go from M_TRGB to M_{TRBG}
    # go from 10^-4 to 10^{-4}
    # put on $ $ if math
    math = 0
    if re.search('_',key):
        math += 1
        tmp = key.split('_')
        tmp1 = '{'+tmp[-1]+'}'
        key = tmp[0]+'_'+tmp1
    if re.search('\^',key):
        math += 1
        tmp = key.split('^')
        tmp1 = '{'+tmp[-1]+'}'
        key = tmp[0]+'^'+tmp1
    if math > 0:
        key = '$'+key+'$'
    return key

def LatexTableComments(filename):
    # print from some text file I've written and supplied the name of
    file = open(filename,'r')
    comment = file.read()
    comment = comment.replace('\n',' ')
    comments = '\\tabelcomments{'+comment + '}'
    return comments

def ReadLatexTable(filename):
    """
    Reads in a latex table that has \startdata and \enddata
    Will figure out what type of datum is being loaded, int, float, str,
    (even binary, hex).
    If string contains \\pm it will split that data and make a new key with
    the correct column heading + ' err'.
    
    try:
    t = ReadLatexTable(filename)
    t.keys()
    """
    file = open(filename,'r')
    data = {}
    colheads = []
    subhead = 0
    i = 0
    lines = file.readlines()
    # run through all the lines, grab the column headings
    # find where the data start and end
    for line in lines:
        if line.startswith('%'): 
            i+=1
            continue
        if re.search('colhead',line):
            if re.search('\\\\',line):
                subhead = 1
            colheads.append(TexTableString(line))
        if re.search('startdata',line):
            dstart = i+1
        if re.search('enddata',line):
            dend = i
        i +=1
    #print dstart,dend
    # if there is a line break in the column heads, fold it.
    # won't work for three lines of column headings...
    if subhead == 1:
        tmp,tmp2 = [],[]
        ncols = len(colheads)
        for j in range(ncols/2):
            tmp = colheads[j]+' '+colheads[j+ncols/2]
            tmp2.append(tmp.rstrip())
        colheads = tmp2
    for colhead in colheads:
        data[colhead] = []
    Ncols = len(colheads)
    Nrows = dend-dstart
    #data = np.ndarray(shape=(Nrows,Ncols), dtype=float)
    for k in range(dstart,dend):
        if lines[k].startswith('%'): continue # skip comments
        if lines[k].startswith('\\'): continue # skip latex syntax lines
        ldata = TextTableLine(lines[k])
        if ldata == [""]: continue # skip empty lines
        for l in range(len(ldata)):
            dat = ldata[l].lstrip()
            dat = dat.rstrip()
            if re.search('tab3',filename): # angst table 3 has lots 
                if dat == '':              # of missing entries
                    tldata = TextTableLine(lines[k-1])
                    dat = tldata[l].lstrip()
                    dat = dat.rstrip()
                    if dat == '':
                        tldata = TextTableLine(lines[k-2])
                        dat = tldata[l].lstrip()
                        dat = dat.rstrip()
            datum = is_numeric(dat)
            data[colheads[l]].append(datum) # data dictionary!
    data = LatexPlusMinus(data) # if \pm, add a err key to dictionary.
    return data

def TextTableLine(line):
    """
    Made to be used with ReadLatexTable
    Cuts the typical crap out of a line from table data in a
    Latex table.
    Splits the data on & and sends back a list.
    """
    return line.replace('\n','').replace('\\\\','').replace('\t','').replace('$','').replace('\\,','').split('&')

def LatexPlusMinus(data):
    """
    Made to be used with ReadLatexTable
    Splits up strings within in dictionary that contain \pm
    before the \pm is updated as the key values
    after the \pm is added as a new key title key + ' err'
    """
    for key in data.keys():
        if re.search('\\pm',str(data[key][0])):
            for m in range(len(data[key])):
                if re.search('\\pm',str(data[key][m])):
                    keyerr = key + ' err'
                    try:
                        dat,er = data[key][m].split('\\pm')
                    except ValueError:
                        data[key][m] = None
                        data[keyerr].append(None)
                    else:
                        datum = is_numeric(dat)
                        err = is_numeric(er)
                        data[key][m] = datum
                        if m == 0:
                            data[keyerr] = [err]
                        else:
                            data[keyerr].append(err)
    return data

def TexTableString(line):
    """
    Made to be used with ReadLatexTable
    Takes latex like:
    \colhead{Exptime} &
    and returns
    Exptime
    If     \colhead{$m_{TRGB}$} &
    will give back m_TRGB
    Not sure what it will do if there are more than two sets of {}
    """
    a = line.split('{')
    if len(a) > 2:
        b = ''.join(a[1:])
    else:
        b = line.split('{')[-1]
    info = b.split('}')[0]
    info = info.replace('\\','')
    info = info.replace('$','')
    info = info.replace(')','')
    info = info.replace('(','')
    info = info.strip()
    return info
