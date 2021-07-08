from numpy import empty,linspace,tri,copy,diff,zeros,interp,exp,sign
import math
import matplotlib.pyplot as plt
#import matplotlib.animation as animation
from matplotlib.animation import FuncAnimation
def RunSimulation(asc_rate, z0, dz, payloads, delay, dt, I, LightningTimes, LightningMag, LightningPos, LightningWidth):
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
        release_time[j].append(60*j + delay)
    for j in range(payloads):
        balloon_z[j] = 3.22376
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
                lightning_dis = -mag * 1e-8 * exp(-width*(z - z[i])**2)
                sigma[:] = sigma - sign(E[i]) * diff(lightning_dis)
                E = calcField(sigma)
                break
    sigma = - diff(I) * delay #this is behaving strangely
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
    #plot(threshold, z[1:])
    #plot(-threshold, z[1:])
    plt.xlabel("E (V/m)")
    plt.ylabel("Altitude (km)")
    plt.close()
    def init():
        line.set_data([], [])
        line2.set_data([], [])
        return line,
    def animate(i):
        line.set_data(E_field[i],z[1:])
        line2.set_data(stored_data[0][:i],stored_height[0][:i])
        return line,
    return FuncAnimation(fig, animate, init_func=init, frames=360, interval=75, blit=True)
def MakePlots(stored_data, stored_height, stored_time):
        payloads = 1
        for j in range(payloads):
            plt.plot(stored_data[j],stored_height[j], label="Payload " + str([j+1]))
        plt.xlabel("E (V/m)")
        plt.ylabel("Altitude (km)")
        plt.legend(loc="upper left")
        plt.ylim(3.22376,11)
    
def SaveAnimationGif(animation): #Caution: this is experimental
    animationToSave = animation
    animationToSave.save('TestAnimation.gif')
def SaveAnimationHTML(animation): #Caution: this is experimental
    animationToSave = animation
    animationToSave.save('TestAnimation.html')