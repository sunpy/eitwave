import scipy
import numpy as np
import matplotlib.pyplot as plt

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
                pos.append(position[i][column])

    x=np.linspace(1,19,19)

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
    w=np.polyfit(x,pos,1)
    #designate p as the polynomial function
    p=np.poly1d(w)

    #plot the polynomial fit over the original position data
    plt.title('Velocity linear fit')
    plt.xlabel('Frame')
    plt.ylabel('Position (deg)')
    plt.plot(x,pos,'g.')
    plt.plot(x,p(x),'r-')
    plt.show()

    
