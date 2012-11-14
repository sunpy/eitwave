import scipy
import numpy
import matplotlib.pyplot as plt

def velocity_histogram(velocity):
    for i in range(0,len(velocity)):
        if velocity[i] != []:
            print('Plotting velocity histogram for frame ' + str(i))
            vel=filter(None,velocity[i])
            plt.hist(vel,range=[-40,10],bins=50)
            plt.title('Velocity Histogram for frame ' + str(i))
            plt.xlabel('Velocity deg/frame')
            plt.ylabel('N')
            plt.show()
            plt.pause(0.5)
        else:
            print('No velocity data for frame ' + str(i) + '. Skipping.')

def velocity_polyfit(velocity, column):
    #wrong! Want to input positions, not velocities.
    vel=[]
    for i in range(0,len(velocity)):
        if velocity[i] == []:
            vel.append(0)
        else:
            if velocity[i][column] == []:
                vel.append(0)
            else:
                vel.append(velocity[i][column])

    x=np.linspace(1,19,19)
    w=np.polyfit(x,vel,1)
    p=np.poly1d(w)

    plt.title('Velocity linear fit')
    plt.xlabel('Frame')
    plt.ylabel('Position (deg)')
    plt.plot(x,vel,'g.')
    plt.plot(x,p(x),'r-')
    plt.show()

    
