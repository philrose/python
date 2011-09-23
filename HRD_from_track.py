#
#
import sys,os
sys.path.append('/Users/Phil/Dropbox/research/python')
from PadovaTracksUtils import *

if len(sys.argv) < 1:
    print 'python HRD_from_track.py track.PMS'
    sys.exit()

track = sys.argv[1]

plt.figure()
plt.xlabel('Log Te')
plt.ylabel('Log L')
t1 = read_tracks(track)
plt.plot(t1['LOG_TE'],t1['LOG_L'],color = 'black')
plt.axis([max(t1['LOG_TE']),min(t1['LOG_TE']),min(t1['LOG_L']),max(t1['LOG_L'])])
plt.title(track)
plt.show()