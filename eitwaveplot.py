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
    plt.title('Velocity fit')
    plt.xlabel('Elapsed time (s)')
    plt.ylabel('Position (deg)')
    plt.plot(x,pos,'g.',ms=10)
    plt.plot(x2,p(x2),'b-',linewidth=2.0)
    #plt.plot(x,line,'r-')
    plt.annotate('v = '+ str(vel) + ' km/s' + ' + ' + str(acc) + ' km/s^2',[200,86])
    plt.show()

    
