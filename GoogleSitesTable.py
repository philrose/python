# GoogleSitesTable.py
"""
python GoogleSitesTable.py imglist outfile
"""
import sys,subprocess
sys.path.append('/astro/users/philrose/python/')

from MatchUtils import *
if len(sys.argv) < 3:
    outfile = 'GoogleSitesTable.html'
else:
    outfile = sys.argv[2]

file = open(outfile,'w')
file.write('<table style=\'border-collapse: collapse; border-color: rgb(136,136,136);border-width=1px;\' border=\'1\' bordercolor=\'#888888\' cellspaceing=\'0\'>\n<tbody>\n')

if len(sys.argv) < 2:
    print 'finding all pngs in this dir'
    import os,glob
    imgs = glob.glob1(os.getcwd(),'*.png')
else:
    infile = open(sys.argv[1],'r')
    imgs = infile.readlines()
    infile.close()

for img in imgs:
    img = img.replace('\n',"")
    print img
    info = subprocess.Popen([r"file",img],stdout=subprocess.PIPE).communicate()[0]
    try:
        wh = info.split(',')[1]
    except IndexError:
        print 'problem with ',img
    w,h = wh.split('x')
    w = int(w.strip())
    h = int(h.strip())
    try:
        PropID,Target,Filter1,Filter2,filename = extract_title(img)
        title = PropID+' '+Target+' '+Filter1+'-'+Filter2
    except IndexError:
        title = img.replace('_',' ').replace('.png','')
    file.write('<tr><td style=\'text-align: center;\'><h3>'+title+'</a></h3></td></tr>\n')
    file.write('<tr><td><div style=\'display: block; text-align: center; margin-right: auto;\'>\n')
    file.write('<a href=\'http://www.astro.washington.edu/users/philrose/Research/'+img+'\'>\n')
    file.write('<img src=\'http://www.astro.washington.edu/users/philrose/Research/'+img+'\' \'border=\'0\' height=\''+str(h/2)+'\' width=\''+str(w/2)+'\'></a></div></td></tr>\n')
    
file.write('</tbody></table>\n')
file.close()

print 'wrote file: ', outfile 
print 'tar -cvf pngs.tar ',' '.join(img for img in imgs)
print 'scp pngs.tar philrose@portal.astro.washington.edu:/www/astro/users/philrose/html/Research/.'
print 'ssh philrose@portal.astro.washington.edu'
print 'cd /www/astro/users/philrose/html/Research/'
print 'tar -xvf pngs.tar'
    
