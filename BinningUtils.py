#
#  BinningUtils.py
#  
#
#  Created by Philip Rosenfield on 9/23/10.
#
import numpy as np
def bin_up(x,y,nbinsx=None,nbinsy=None):
    """
    Adapted from contour_plus writen by Oliver Fraser. 
    Commented out a means to bin up by setting dx,dy. I found it more useful
    to choose the number of bins in the x and y direction (nbinsx,nbinsy)
    
    # xrange = np.arange(xmin, xmax, dx)
    # yrange = np.arange(ymin, ymax, dy)
    # nbinsx = xrange.shape[0]
    # nbinsy = yrange.shape[0]
    
    """
    npts = len(x)
    xmin  = float( x.min() ) # cast these for matplotlib's w/o np support
    xmax  = float( x.max() )
    ymin  = float( y.min() )
    ymax  = float( y.max() )
    if nbinsx == None:
        nbinsx = 60
    if nbinsy == None:
        nbinsy = 60
    dx = (xmax - xmin) / nbinsx
    dy = (ymax - ymin) / nbinsy
    xrange = np.linspace(xmin, xmax, nbinsx)
    yrange = np.linspace(ymin, ymax, nbinsy)
    Z    = np.zeros((nbinsy, nbinsx))
    xbin = np.zeros(npts)
    ybin = np.zeros(npts)
    # assign each point to a bin
    for i in range(npts):
        xbin[i] = min(int((x[i] - xmin) / dx), nbinsx-1)
        ybin[i] = min(int((y[i] - ymin) / dy), nbinsy-1)
        # and count how many are in each bin
        Z[ ybin[i] ][ xbin[i] ] += 1
        
    return Z,dx,dy

def PtsWithinContour(contourset,levels,Z):
    inds = []
    for level in levels:
        Pind = []
        pinds = []
        paths_inside= []
        paths = contourset.collections[level].get_paths()
        npaths = len(paths)
        for i in range(npaths):
            path_inside = []
            path = paths[i]
            #print 'Path ',i
            #for verts in path.vertices:
            # x23 = plt.plot(verts[0],verts[1],'o',color=colors[i])
            pind = np.nonzero(nxutils.points_inside_poly(Z, path.vertices))[0]
            for j in range(npaths):
                #print 'Going over Path ',j,' to see if it is inside of Path ',i
                if j != i: # skip over identical
                    # are there contours within contours?
                    mask = nxutils.points_inside_poly(paths[j].vertices, path.vertices)*1
                    if sum(mask) > 0:
                        path_inside.append(j)
                        paths_inside.append(j)
            for k in path_inside:
                sind = np.nonzero(nxutils.points_inside_poly(Z, paths[k].vertices))[0]
                for s in sind:
                    try:
                        dels = np.nonzero(pind == s)[0][0]
                        pind = np.delete(pind,dels)
                    except IndexError:
                        continue
            pinds.append(pind)
        for l in np.delete(range(npaths),paths_inside):
            Pind.append(pinds[l])
        ind = concatenate(Pind)
        inds.append(ind)
    return inds



def fitellipseCaller(verts):
    theta = np.arange(0,2.*np.pi,0.05)
    c,A,B,alpha = fitellipse(verts,'linear')
    xvals = c[0] + A*cos(theta)*cos(alpha)-B*sin(theta)*sin(alpha)
    yvals = c[1] + A*cos(theta)*sin(alpha)+B*sin(theta)*cos(alpha)
    return xvals, yvals,c


def inside(x,y,u,v,verts=False):
    """
    returns the indices of u,v that are within the boundries of x,y.
    """
    if verts != True:
        verts = get_verts(x,y,nbinsx=60,nbinsy=60,smooth=1)
    else:
        verts =  np.column_stack((x,y))
    points = np.column_stack((u,v))
    mask = nxutils.points_inside_poly(points, verts)
    ind = np.nonzero(mask)[0]
    return ind
    
def get_verts(x,y,dx=None,dy=None,nbinsx=None,nbinsy=None,smooth=None):
    if smooth == None:
        smooth = 0
    else:
        smooth = 1
    ymin = y.min()
    ymax = y.max()
    xmin = x.min()
    xmax = x.max()
    
    if dx == None and dy == None:
        dx = (xmax-xmin)/nbinsx
        dy = (ymax-ymin)/nbinsy
        
    if nbinsx == None and nbinsy == None:
        nbinsy = (ymax-ymin)/dy
        nbinsx = (xmax-xmin)/dx
    
    ymid = []
    min_x = [] 
    max_x = []
    for j in range(nbinsy):
        yinner = ymin+j*dy
        youter = ymin+(j+1)*dy
        # counter intuitive because I'm dealing with mags...
        ind = np.nonzero((y > yinner) & (y < youter))[0]  
        if len(ind) > 0:
            if smooth == 1:
                min_x.append(np.average(x[ind])-3.*np.std(x[ind]))
                max_x.append(np.average(x[ind])+3.*np.std(x[ind]))
                ymid.append((yinner+youter)/2.)
            else:
                min_x.append(min(x[ind]))
                max_x.append(max(x[ind]))
                ymid.append((yinner+youter)/2.)
    
    max_x.reverse() 
    ymidr = ymid[:]
    ymidr.reverse()
    
    # close polygon
    max_x.append(min_x[0])
    
    # close polygon
    ymidr.append(ymid[0])
    
    # get verticies of polygon
    xs = np.concatenate((min_x,max_x))
    ys = np.concatenate((ymid,ymidr))
    verts = np.column_stack((xs,ys))
    
    return verts
