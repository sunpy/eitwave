from scipy import stats
import numpy as np
import matplotlib.pyplot as plt
import copy

def velocity_histogram(velocity):
    """Plot a histogram of the wavefront velocity values for each eit wave image frame"""
    for i in range(0,len(velocity)):
        if velocity[i] != []:
            print('Plotting velocity histogram for frame ' + str(i))
            vel=filter(None,velocity[i])
            plt.hist(vel,range=[-40,10],bins=50)
            plt.title('Velocity Histogram for frame ' + str(i))
            plt.xlabel('Velocity (deg/frame)')
            plt.ylabel('N')
            plt.show()
            plt.pause(0.5)
        else:
            print('No velocity data for frame ' + str(i) + '. Skipping.')

def width_histogram(width):
    """Plot a histogram of the wavefront width values for each eit wave image frame"""
    for i in range(0,len(width)):
        if width[i] != []:
            print('Plotting width histogram for frame ' + str(i))
            wid=filter(None,width[i])
            plt.hist(wid,range=[0,20],bins=60)
            plt.title('Width Histogram for frame ' + str(i))
            plt.xlabel('Width (deg)')
            plt.ylabel('N')
            plt.show()
            plt.pause(0.5)
        else:
            print('No width data for frame ' + str(i) + '. Skipping')

def mean_width_vs_time(width,maps):

    w=[]
    times=[]
    m=maps[0]
    base_time=m.date
    for i in range(0,len(width)):
        m=maps[i]
        sec=m.date-base_time
        times.append(sec.seconds)
        if width[i] == []:
            w.append(np.nan)
        else:
            noz=np.nonzero(width[i])
            noz=np.array(noz)
            w_tmp=np.array(width[i])
            w.append(w_tmp[noz].mean())
    x=times
    #ignore NaN values
    u=np.isnan(w)
    u=np.invert(u)
    #convert from list to numpy array to avoid type error
    w=np.array(w)
    x=np.array(x)
    #keep everything that's not NaN
    w=w[u]
    x=x[u]
    x2=np.linspace(0,x[-1],x[-1])

    #make a plot of the width vs time.
    plt.title('Mean wavefront width')
    plt.xlabel('Elapsed time (s)')
    plt.ylabel('Mean width (deg)')
    plt.plot(x,w,'b-',ms=6)
    plt.plot(x,w,'bs',ms=6)
    #plt.plot(x2,p(x2),'b-',linewidth=2.0)
    #plt.plot(x,line,'r-')
    #plt.annotate('v = '+ str(vel) + ' km/s' + ' + ' + str(acc) + ' km/s^2',[200,86])
    plt.show()

def amplitude_vs_time(wavefront,maps,column):
    a=[]
    times=[]
    m=maps[0]
    base_time=m.date
    for i in range(0,len(wavefront)):
        m=maps[i]
        sec=m.date-base_time
        times.append(sec.seconds)
        if wavefront[i] == []:
            a.append(np.nan)
        else:
            w_tmp=wavefront[i][0:,column]
            if w_tmp.max() == 0:
                a.append(np.nan)
            else:
                a.append(w_tmp.max())
            
    x=times
    #ignore NaN values
    u=np.isnan(a)
    u=np.invert(u)
    #convert from list to numpy array to avoid type error
    a=np.array(a)
    x=np.array(x)
    #keep everything that's not NaN
    a=a[u]
    x=x[u]
    x2=np.linspace(0,x[-1],x[-1])

    #make a plot of the width vs time.
    plt.title('Amplitude vs time')
    plt.xlabel('Elapsed time (s)')
    plt.ylabel('Amplitude (arb.)')
    plt.plot(x,a,'b-',ms=6,linewidth=2.0)
    plt.plot(x,a,'bs',ms=4)
    #plt.plot(x2,p(x2),'b-',linewidth=2.0)
    #plt.plot(x,line,'r-')
    #plt.annotate('v = '+ str(vel) + ' km/s' + ' + ' + str(acc) + ' km/s^2',[200,86])
    plt.show()

def width_vs_time(width,maps,column):
    w=[]
    times=[]
    m=maps[0]
    base_time=m.date
    for i in range(0,len(width)):
        m=maps[i]
        sec=m.date-base_time
        times.append(sec.seconds)
        if width[i] == []:
            w.append(np.nan)
        else:
            if width[i][column] == []:
                w.append(np.nan)
            else:
                w.append(width[i][column])
    x=times
    #ignore NaN values
    u=np.isnan(w)
    u=np.invert(u)
    #convert from list to numpy array to avoid type error
    w=np.array(w)
    x=np.array(x)
    #keep everything that's not NaN
    w=w[u]
    x=x[u]
    x2=np.linspace(0,x[-1],x[-1])

    #make a plot of the width vs time.
    plt.title('Wavefront width')
    plt.xlabel('Elapsed time (s)')
    plt.ylabel('Width (deg)')
    plt.plot(x,w,'bs',ms=10)
    #plt.plot(x2,p(x2),'b-',linewidth=2.0)
    #plt.plot(x,line,'r-')
    #plt.annotate('v = '+ str(vel) + ' km/s' + ' + ' + str(acc) + ' km/s^2',[200,86])
    plt.show()
    
    

def width_vs_longitude(width,ref_map):
    """Plot the width of the wavefront as a function of longitude for a given image frame"""
    #convert list to numpy array
    width2=np.array(width)
    #generate the x-axis
    xlimits=ref_map.xrange
    xsize=ref_map.shape[1]
    xcoord=np.linspace(xlimits[0],xlimits[1],xsize)

    #need to remove null values in numpy array to avoid error. Get indices of nonzero elements
    indices=width2.nonzero()

    #make plot
    plt.title('Wavefront width vs longitude')
    plt.xlabel('Heliographic longitude (deg.)')
    plt.ylabel('wavefront width (deg.)')
    plt.ylim( (0,3))
    plt.plot(xcoord[indices],width2[indices],'bs')
    #plt.plot(xcoord[indices],width2[indices],'-',linewidth=2.0)
    plt.show()

def velocity_polyfit(position, maps, column):
    """Perform a polynomial fit to the wavefront position as a function of time, for a given image column"""
    pos=[]
    times =[]
    m=maps[0]
    base_time=m.date
    for i in range(0,len(position)):
        m=maps[i]
        sec=m.date-base_time
        times.append(sec.seconds)
        if position[i] == []:
            pos.append(np.nan)
        else:
            if position[i][column] == []:
                pos.append(np.nan)
            else:
                #if i == 1:
                    #position[i][column] = []:
                    #   pos.append(np.nan)
                    # else:
                pos.append(position[i][column])

    
                #xlen=len(pos)+1
                #x=np.linspace(1,xlen,xlen)
    x=times
    #ignore NaN values
    u=np.isnan(pos)
    u=np.invert(u)
    #convert from list to numpy array to avoid type error
    pos=np.array(pos)
    x=np.array(x)
    #keep everything that's not NaN
    pos=pos[u]
    x=x[u]
    x2=np.linspace(0,x[-1],x[-1])

    #perform a polynomial fit to the data (pos). x is the independent variable,
    #and the third parameter is the order of the polynomial.
    w=np.polyfit(x,pos,2,full=True)
    #print w
    #designate p as the polynomial function
    p=np.poly1d(w[0])

    #try linear regression routine from scipy.stats
    #slope, intercept, r_value, p_value, std_err = stats.linregress(x,pos)
    #line = slope*x + intercept

    #print slope, intercept, std_err
    #work out the velocity in sensible units
    vel=w[0][1]
    acc=w[0][0]
    #vel=copy.copy(slope)
    #in km/s
    vel=(vel*1.21e4)
    acc=(acc*1.21e4)
    vel=round(vel,1)
    acc=round(acc,1)
    #plot the polynomial fit over the original position data
    plt.title('Velocity fit',fontsize=20)
    plt.tick_params(labelsize=20)
    plt.xlabel('Elapsed time (s)',fontsize=20)
    plt.ylabel('Position (deg)',fontsize=20)
    plt.plot(x,pos,'g.',ms=10)
    plt.plot(x2,p(x2),'b-',linewidth=2.0)
    #plt.plot(x,line,'r-')
    plt.annotate('v = '+ str(vel) + ' km/s',[500,86],fontsize=20)
    plt.annotate('a = '+ str(acc) + 'km/s^2',[500,83],fontsize=20)
    plt.show()

    
