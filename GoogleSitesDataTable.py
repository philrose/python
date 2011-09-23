#
#  GoogleSitesDataTable.py
#  
#
#  Created by Philip Rosenfield on 9/23/10.
#
import sys

filename = sys.argv[1]

f = open(filename,'r')
lines = f.readlines()
f.close()

o = open(filename+'.html','w')
o.write('<br><table style="border-collapse: collapse; border-color: rgb(136, 136, 136); border-width: 1px;" border="1" bordercolor="#888888" cellspacing="0">\n<tbody>\n')
for line in lines:
    if line.startswith('#'):
        line = line.replace('#','')
        heads = line.split()
        o.write('<tr>')
        continue
        for head in heads:
            o.write('<th> %s </th>' % head)
        o.write('</tr> \n')
    row = line.split()
    o.write('<tr>')
    for item in row:
        o.write('<td> %s </td>' % item.strip())
    o.write('</tr>\n')

o.write('</tbody>\n</table>\n')
o.close()
print 'wrote '+filename+'.html'