'''
wrapper to run Leo's bccode.
BCDIR must be set in .cshrc
This works better than call_bccode.py.
'''

import sys,os
sys.path.append(os.environ['PYTHON_CODES'])
from GenUtils import get_afile
BCDIR = os.environ['BCDIR']

if len(sys.argv) < 3: 
    print 'usage python color_tracks.py filter res track_location'
    print 'if res = 0 full resoltion res>0 set resolution'
    print 'ex phat 0.0 /Users/Phil/research/PHAT/tracks/S11_Z0.017_Y0.279/'
    sys.exit()

filter = sys.argv[1]
res = float(sys.argv[2])
track_loc(sys.argv[3])

tracks = get_afile(track_loc,'*PMS*')

Tracks = [track+'_'+str(res) for track in tracks if res !=0]
if len(Tracks) == 0: Tracks = tracks


[os.system(' '.join((os.path.join(BCDIR,'main'),filter,track,str(res)))) for track in Tracks]