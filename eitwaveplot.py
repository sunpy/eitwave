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
    

def velocity_polyfit(position, column):
    """Perform a polynomial fit to the wavefront position as a function of time, for a given image column"""
    pos=[]
    for i in range(0,len(position)):
        if position[i] == []:
            pos.append(np.nan)
        else:
            if position[i][column] == []:
                pos.append(np.nan)
            else:
                if i == 1:
                    #position[i][column] = []:
                    pos.append(np.nan)
                else:
                    pos.append(position[i][column])

    xlen=len(pos)+1
    x=np.linspace(1,xlen,xlen)
    #ignore NaN values
    u=np.isnan(pos)
    u=np.invert(u)
    #convert from list to numpy array to avoid type error
    pos=np.array(pos)
    #keep everything that's not NaN
    pos=pos[u]
    x=x[u]

    #perform a polynomial fit to the data (pos). x is the independent variable,
    #and the third parameter is the order of the polynomial.
    w=np.polyfit(x,pos,1,full=True)
    #print w
    #designate p as the polynomial function
    p=np.poly1d(w[0])

    #try linear regression routine from scipy.stats
    #slope, intercept, r_value, p_value, std_err = stats.linregress(x,pos)
    #line = slope*x + intercept

    #print slope, intercept, std_err
    #work out the velocity in sensible units
    vel=w[0][0]
    #vel=copy.copy(slope)
    #in km/s
    vel=(vel*1.21e4 / 12.0)
    vel=round(vel,2)
    #plot the polynomial fit over the original position data
    plt.title('Velocity linear fit')
    plt.xlabel('Frame')
    plt.ylabel('Position (deg)')
    plt.plot(x,pos,'g.')
    plt.plot(x,p(x),'r-')
    #plt.plot(x,line,'r-')
    plt.annotate('Velocity = '+ str(vel) + ' km/s',[10,86])
    plt.show()

    