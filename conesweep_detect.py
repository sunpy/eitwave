
#from visualize import visualize
#from sim import wave2d
#from skimage.transform import hough
import numpy as np
import sunpy
import os

def linesampleindex(a,b,np = 1000):
    """ Get the indices in an array along a line"""
    x, y = np.linspace(a[0],b[0],np), np.linspace(a[1],b[1],np)
    xi = x.astype(np.int)
    yi = y.astype(np.int)
    return xi, yi

def cone(a,b,c,d,nperp=100,npar=1000):
    """ define a simple trapezoidal cone (a,b) opens up to (c,d) """
    abxi, abyi = linesampleindex(a,b,np = nperp)
    cdxi, cdyi = linesampleindex(c,d,np = nperp)
    sxi = np.zeros((nperp,npar))
    syi = np.zeros((nperp,npar))
    for i in range(0,nperp):
        sxi[i,:], syi[i,:] = linesampleindex( (abxi[i],abyi[i]), (cdxi[i],cdyi[i]), np = npar)
    return sxi, syi

def sumalongcone(img,a,b,c,d,nperp=100,npar=1000):
    """sum along the perpendicular direction a -> b, in parallel with c -> d"""
    sxi, syi = cone(a,b,c,d,nperp=nperp, npar=npar)
    zi = img(sxi,syi)
    return zi.sum(axis=0)

def conevals(center,orientation,r1,openangle1,r2,openangle2):
    """ define the four points of the cone given the center of
    the circle it emanates from, the orientation of the cone, and the
    inner and outer opening angles.  Note that r2 must be greater than r1 
    and openangle2 >= openangle1 """
    
    ax = center[0] + r1*np.cos(orientation - openangle1)
    ay = center[1] + r1*np.sin(orientation - openangle1)
    
    bx = center[0] + r1*np.cos(orientation + openangle1)
    by = center[1] + r1*np.sin(orientation + openangle1)
    
    cx = center[0] + r2*np.cos(orientation - openangle1)
    cy = center[1] + r2*np.sin(orientation - openangle1)
    
    dx = center[0] + r2*np.cos(orientation + openangle2)
    dy = center[1] + r2*np.sin(orientation + openangle2)

    return (ax,ay), (bx,by), (cx,cy), (dx,dy)


def htLine(distance,angle,img):
    shape = img.shape
    ny = shape[0]
    nx = shape[1]
    eps = 1.0/float(ny)

    if abs(np.sin(angle)) > eps:
        gradient = - np.cos(angle) / np.sin(angle)
        constant = distance / np.sin(angle)
        for x in range(0,nx):
            y = gradient*x + constant
            if y <= ny-1 and y >= 0:
                img[y,x] = 255
            else:
                img[:,distance] = 255

    return img

def main():

    
    # Lots of big images.  Need to be smart about how to handle the data
    
    # load in the data with a single EIT wave
    directory = "/home/ireland/eitwave_data/jp2/20110601_02_04/"
    extension = ".jp2"
    # get the file list
    lst = []
    for f in os.listdir(directory):
        if f.endswith(extension):
            lst.append(f)
    lll = sorted(lst)
    # accumulate every "naccum" files
    naccum = 2
    
    # accumulate using superpixels
    nsuper = 4
    
    # number of files to consider
    #nfiles = len(lll)
    
    nfiles = 10
    
    # storage for all the maps
    maps = []
    
    # accumulate naccum files, super-pixel them, then add them together.
    j = 0
    while j+naccum <= nfiles:
        i = 0
        print('\n Starting new accumulation:')
        while i < naccum:
            filename = os.path.join(directory,lll[j+i])
            print('  File %(#)i out of %(nfiles)i' % {'#':i+j, 'nfiles':nfiles})
            print('  Reading in file '+filename)
            map1 = (sunpy.make_map(filename)).superpixel((nsuper,nsuper))
            if i == 0:
                m = map1
            else:
                m = m + map1
            i = i + 1
        j = j + naccum
        maps.append(m)

    #z = wave2d.transform(params,maps)
    
    # number of running differences
    ndiff = len(maps)-1
    
    # Each JP2 file has a maximum of 255
    #maxval = 255 * nsuper * nsuper * naccum
    
    # difference threshold as a function of the maximum value
    #diffthresh = 0.01 * maxval #300
    
    diffs = []

    # angle for each sweep
    sweepstep = 10.0
    norient = 360.0/sweepstep
    
    # half cone width
    theta = sweepstep/2.0
    
    # inner radius in pixels
    innerRadius = 10
    
    # maximum outer radius in pixels
    outerRadius = 300
    
    # maximum number of points to sample radially
    npar = 300
    
    # storage for the summing
    
    flareloc = (500,500)
    summed = np.zeros((ndiff,norient,npar))
    
    # perform the analysis
    for i in range(0,ndiff):
        
        # take the difference
        diffmap = maps[i+1] -maps[i]  #> diffthresh
    
        # keep
        diffs.append(diffmap)
    
        # extract the image from the storage array
        img = diffmap

        # perform the cone sweep
        for j in range(0,norient):

            orientation = j*sweepstep

            # calculate the convalues
            a, b, c, d = conevals(flareloc, orientation, innerRadius, theta, outerRadius, theta)
            
            summed[i,j,:] = sumalongcone(diffmap,a,b,c,d, npar = npar)
            
    
        # Perform the hough transform on each of the difference maps
        #transform,theta,d = hough(img)
    
        # Filter the hough transform results and find the best lines
        # in the data
        #indices =  (transform >votethresh).nonzero()
        #distances = d[indices[0]]
        #theta = theta[indices[1]]
        #n =len(indices[1])
        #print n
    
        # Perform the inverse transform to get a series of rectangular
        # images that show where the wavefront is.
        #invTransform = sunpy.make_map(maps[i+1])
        #invTransform.data = np.zeros(imgShape)
        #for i in range(0,n):
        #    nextLine = htLine( distances[i],theta[i], np.zeros(shape=imgShape) )
        #    invTransform = invTransform + nextLine
    
        # Dump the inverse transform back into a series of maps
        #detection.append(invTransform)
    
    
    #visualize(detection)

if __name__ == '__main__':
    main()

