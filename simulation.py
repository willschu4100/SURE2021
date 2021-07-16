from numpy import empty,linspace,tri,copy,diff,zeros,interp,exp,sign,argmax,max
import math
import matplotlib.pyplot as plt
#import matplotlib.animation as animation
from matplotlib.animation import FuncAnimation
from pandas import DataFrame
def RunSimulation(asc_rate, z0, dz, payloads, delay, startTime, dt, I, LightningTimes, LightningMag, LightningPos, LightningWidth):
    t = 0
    jj = 0
    balloon_z = empty(payloads)
    n = 360
    tmax = 4600
    e0 = 8.854187817e-12
    z = linspace(0,16,n+1)
    #threshold = 2.84e5*math.exp(-z[1:]/8)   this probably isn't necessary
    I = I(z)
    a = tri(n, k=n) + tri(n) + tri(n,k=-1)
    A = a*(1/(2*e0))
    release_time = [[] for i in range(payloads)]
    for j in range(payloads):
        release_time[j].append(delay*j + startTime)
    for j in range(payloads):
        balloon_z[j] = z0
    stored_data = [[] for i in range(payloads)]
    stored_height = [[] for i in range(payloads)]
    stored_time = [[] for i in range(payloads)]
    def calcField(sigma):
        return A.dot(sigma)
    E = zeros(n)
    sigma = zeros(n)
    E_field = []
    sigma_array = []
    
    def TimeDischarge(sigma, E, pos, mag, width):
        for i in range(n):
            if(z[i] >= pos):
                #lightning_dis = -mag * 1e-8 * exp(-width*(z - z[i])**2)
                lightning_dis = -mag * 1e-8 * exp(-(z - z[i])**2/(2*width**2))
                sigma[:] = sigma - sign(E[i]) * diff(lightning_dis)
                E = calcField(sigma)
                break
    sigma = - diff(I) * delay #this sometimes behaves strangely in the animation
    while t < tmax:
        sigma = sigma - diff(I)*dt
        E = calcField(sigma)
        E_field.append(copy(E))
        i = 0    
        for i in range (len(LightningTimes)):
            if t >= LightningTimes[i]:
                TimeDischarge(sigma,E, LightningPos[i], LightningMag[i], LightningWidth[i])
                LightningTimes.pop(i)
                LightningMag.pop(i)
                LightningPos.pop(i)
                LightningWidth.pop(i)
                break
        sigma_array.append(copy(sigma))
        balloon_data = empty((tmax-delay)//dt)
        for j in range(payloads):
            if t < release_time[j][0]:
                balloon_data[j] = 0
                balloon_z[j] = 0
                stored_data[j].append(balloon_data[j])
                stored_height[j].append(balloon_z[j])
                stored_time[j].append(t)
            if t > release_time[j][0]:
                balloon_data[j] = interp(balloon_z[j],z[1:],E_field[jj])
                balloon_z[j] = balloon_z[j] + asc_rate*dt
                stored_data[j].append(balloon_data[j])
                stored_height[j].append(balloon_z[j])
                stored_time[j].append(t)      
        jj = jj + 1  
        t = t + dt
        
    return(stored_data, stored_height, stored_time, E_field, z)

def RunAnimation(stored_data, stored_height, E_field, z):
    fig = plt.figure()
    ax = plt.axes(xlim=(-2e5,2e5), ylim=(0,16))
    line, = ax.plot([], [], lw=2)
    line2, = ax.plot([], [], lw=2)
    line3, = ax.plot([], [], lw=2)
    #plot(threshold, z[1:])
    #plot(-threshold, z[1:])
    plt.xlabel("E (V/m)")
    plt.ylabel("Altitude (km)")
    plt.close()
    def init():
        line.set_data([], [])
        line2.set_data([], [])
        line3.set_data([], [])
        return line,
    def animate(i):
        line.set_data(E_field[i],z[1:])
        line2.set_data(stored_data[0][:i],stored_height[0][:i])
        line3.set_data(stored_data[1][:i],stored_height[1][:i])
        #to get this to work, i might be able to make it so payloads launched
        #after the first payload have an altitude of zero until their launch time. not sure how to get this to work however
        return line,
    return FuncAnimation(fig, animate, init_func=init, frames=360, interval=75, blit=True)
#to animate multiple payloads, set payload data as a list of payload data. loop over list, make plots for each, store results
#this probably isn't worth it however
def MakePlots(stored_data, stored_height, stored_time, payloads):
    for j in range(payloads):
        plt.plot(stored_data[j],stored_height[j], label="Payload " + str(j+1))
        #plt2=plt.twinx(plt)
        #plt2.plot(stored_data[j],stored_time[j])
    plt.xlabel("E (V/m)")
    plt.ylabel("Altitude (km)")
    plt.legend(loc="upper left")
    plt.ylim(stored_height[0][0],11)
    #plt.axes().set_aspect(30000)
def PlotMultipleSims(stored_data1, stored_height1, stored_time1, stored_data2, stored_height2, stored_time2, payloads):
    for j in range(payloads):
        plt.plot(stored_data1[j],stored_height1[j], label="Simulation 1")
        plt.plot(stored_data2[j],stored_height2[j], label="Simulation 2")
    plt.xlabel("E (V/m)")
    plt.ylabel("Altitude (km)")
    plt.legend(loc="upper left")
    plt.ylim(stored_height1[0][0],11)        
def MakePlotsReversed(stored_data, stored_height, stored_time, payloads):
    for j in range(payloads):
        plt.plot(stored_height[j], stored_data[j], label="Payload " + str(j+1))
    plt.xlabel("Altitude (km)")
    plt.ylabel("E (V/m)")
    plt.legend(loc="upper left")
    plt.xlim(stored_height[0][0],11)
    plt.axes().set_aspect(1/30000)
def SaveAnimationGif(animation): #Caution: this is experimental
    animationToSave = animation
    animationToSave.save('TestAnimation.gif')
def SaveAnimationHTML(animation): #Caution: this is experimental. It also does not work
    animationToSave = animation
    animationToSave.save('TestAnimation.html')

def InferCurrent(stored_data, stored_height, stored_time): #This is still a work in progress. I'm not even sure what to do with this. i is some generic layer. i could probably find some relationship between layers and height and use that instead of i to make it easier to understand
    e0 = 8.854187817e-12
    #use stored data from multiple payloads to find some change in electric field over some change in time, deltaE/deltaT is proportional to current density, J
    df = DataFrame(list(zip(stored_data[0], stored_time[0], stored_height[0])), columns = ['e', 't', 'z'])
    df2 = DataFrame(list(zip(stored_data[1], stored_time[1], stored_height[1])), columns = ['e', 't', 'z'])
    #dataE = interp(stored_data, stored_height)#Interpolate payload data at altitude, interpolate time at altitude, use those to figure out e field change
    #dataTime = interp(stored_time, stored_height)
    gaussianPeak = argmax(abs(df['e']))
    gaussianPeak2 = argmax(abs(df2['e']))
    gaussianPeakE = - max(abs(df['e']))
    gaussianPeakE2 = - max(abs(df2['e']))
    centerGaussian = df.iloc[gaussianPeak]['z']
    dE = gaussianPeakE2 - gaussianPeakE
    dT = df2.iloc[gaussianPeak2]['t'] - df.iloc[gaussianPeak]['t'] 
    dEdT = dE/dT
    #width = #i'm not sure how to find the width. 
    magnitude = -e0 * dEdT
    #print(magnitude * exp(('z' - centerGaussian)**2 / (2 * width**2))
    print(centerGaussian, dEdT, magnitude)