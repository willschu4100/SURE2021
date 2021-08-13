from numpy import empty,linspace,tri,copy,diff,zeros,interp,exp,sign,argmax,max, array
import math
import matplotlib.pyplot as plt
#import matplotlib.animation as animation
from matplotlib.animation import FuncAnimation
from pandas import DataFrame
import scipy.optimize as op
def RunSimulation(asc_rate, z0, dz, payloads, delay, startTime, dt, I, LightningTimes, LightningMag, LightningPos, LightningWidth):
    t = 0
    jj = 0
    balloon_z = empty(payloads)
    n = 360
    tmax = 4600
    e0 = 8.854187817e-12
    z = linspace(0,16,n+1)
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
                lightning_dis = -mag * 1e-8 * exp(-(z - z[i])**2/(2*width**2))
                if(E[i] == 0):
                    sigma[:] = sigma - diff(lightning_dis)
                elif(E[i] != 0):
                    sigma[:] = sigma - sign(E[i]) * diff(lightning_dis)
                E = calcField(sigma)
                break
    sigma = - diff(I) * delay
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
            #the following code is needed for using the animation with multiple payloads, not necessary otherwise
            #if t < release_time[j][0]: 
                #balloon_data[j] = 0
                #balloon_z[j] = 0
                #stored_data[j].append(balloon_data[j])
                #stored_height[j].append(balloon_z[j])
                #stored_time[j].append(t)
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
    ax = plt.axes(xlim=(-140000,45000), ylim=(stored_height[0][0],11))
    line, = ax.plot([], [], lw=2, label = "Electric Field")
    line2, = ax.plot([], [], lw=2, label = "Payload 1 Data")
    #line3, = ax.plot([], [], lw=2, label = "Payload 2 Data")
    plt.xlabel("E (V/m)")
    plt.ylabel("Altitude (km)")
    plt.close()
    plt.legend(loc="upper left")
    def init():
        line.set_data([], [])
        line.set_color("orange")
        line2.set_data([], [])
        #line2.set_color("blue")
        #line3.set_data([], [])
        return line,
    def animate(i):
        line.set_data(E_field[i],z[1:])
        line2.set_data(stored_data[0][:i],stored_height[0][:i])
        #line3.set_data(stored_data[1][:i],stored_height[1][:i]) #this animates multiple payloads. there's probably a better way to do it
        return line,
    return FuncAnimation(fig, animate, init_func=init, frames=260, interval=75, blit=True)
def MakePlots(stored_data, stored_height, stored_time, payloads):
    for j in range(payloads):
        plt.plot(stored_data[j],stored_height[j], label="Payload " + str(j+1))
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
def SaveAnimationGif(animation):
    animationToSave = animation
    animationToSave.save('TestAnimation.gif')


def CalculateCurrent(stored_data, stored_height, stored_time, delay, dt, asc_rate, LightningTimes): 
    e0 = 8.854187817e-12
    n = 10
    I = []
    payload1 = DataFrame(list(zip(stored_data[0], stored_time[0], stored_height[0])), columns = ['e', 't', 'z'])#measurements taken by payload 1
    for a in range(int(delay/dt)):
        payload1 = payload1.iloc[:-1, :] #since payload 1 launches before payload 2, it has extra entries at the end that don't matter
    payload2 = DataFrame(list(zip(stored_data[1], stored_time[1], stored_height[1])), columns = ['e', 't', 'z'])#measurements taken by payload 2
    posUpper=0 #position of payload 1 when lightning strike occurs. set to 0 by default
    posLower=0 #position of payload 2 when lightning strike occurs. set to 0 by default
    layerUpper = 0 #layer of payload 1's dataset when lightning strike occurs. set to 0 by default
    layerLower = 0 #layer of payload 2's dataset when lightning strike occurs
    correctionFactor = []
    for k in range(len(payload1)):
        if(len(LightningTimes) >= 1): #this is mainly just to prevent errors when the lightning time list is empty
            if(payload2.iloc[k]['t'] > LightningTimes[0]):
                posLower = round(payload2.iloc[k]['z'], 4)   #the next few lines calculate where the lightning strike occurs for each payload
                posUpper = round(payload1.iloc[k + int(delay/dt)]['z'], 4)
                layerUpper = k + int(delay/dt)
                layerLower = k 
                LightningTimes.pop(0)
                #print("pop")
                dECurrentLower = payload2.iloc[layerLower-1,0] - payload2.iloc[layerLower-2,0]#change in electric field right before lightning strikes. this change in electric field is solely due to current, and is used to estimate the change in electric field solely due to lightning afterwards
                dECurrentUpper = payload1.iloc[layerUpper+1,0] - payload1.iloc[layerUpper,0]
                dELower = (payload2.iloc[layerLower,0] - payload2.iloc[layerLower-1,0])#change in electric field right as lightning strikes, as recorded by payload 2. this change in electric field is due to both lightning discharging the storm and current charging the storm
                dELightningLower = dELower - dECurrentLower#change in electric field right as lightning strikes, solely due to lightning. this subtracts the change in electric field solely due to current to give an estimate of the change in electric field due to lightning
                dEUpper = (payload1.iloc[layerUpper,0] - payload1.iloc[layerUpper-1,0])
                dELightningUpper = dEUpper -  dECurrentUpper
                #print("dECurrentLower is", dECurrentLower, "dECurrentUpper is", dECurrentUpper)
                #print("dELower is" ,dELower, "dEUpper is", dEUpper)
                #print("dELightningLower is", dELightningLower, "dELightningUpper is", dELightningUpper)
        if((payload2.iloc[k]['z'] >= posLower) & (payload2.iloc[k]['z'] < posUpper)):
            correctionFactor = interp(payload2.iloc[k]['z'], [posLower, posUpper], [dELightningLower, dELightningUpper]) #correction factor is an interpolation, using the upper and lower positions of the payloads when lightning strikes, and using the change in electric fields to directly change the electric field itself
            payload2.iloc[k,0] = payload2.iloc[k,0] -  correctionFactor
    for i in range(len(payload2)):
        I.append(-e0 * (payload2.iloc[i]['e'] - payload1.iloc[i]['e']) / delay)
    z = linspace(0, len(payload1)*asc_rate*dt, len(payload2))
    plt.plot(z, I, label = "Inferred Current")
    
def CurveFit(stored_data, stored_height, delay, asc_rate, z0, dz, payloads, startTime, dt, LightningTimes, LightningMag, LightningPos, LightningWidth, guessMagLower, guessWidthLower, guessLocLower, guessMagMid, guessWidthMid, guessLocMid, guessMagHigh, guessWidthHigh, guessLocHigh):#this is probably unnecessary. CurveFitLightning doesn't handle guesses of 0 for lightning parameters very well, but this can be fixed by giving it very small guesses
    def FitSim(x, mag, width, loc, mag2, width2, loc2, mag3, width3, loc3):
        def myI(zz):
            return mag*exp(-(zz-loc)**2/width) + mag2*exp(-(zz-loc2)**2/width2) + mag3*exp(-(zz-loc3)**2/width3)
        data, a, b, c, d = RunSimulation(asc_rate, z0, dz, payloads, delay, startTime, dt, myI, LightningTimes, LightningMag, LightningPos, LightningWidth)
        return data[0]
    #return FitSim(0, 5e-10, .5, 5.7)
    ans= op.curve_fit(FitSim, stored_height[0], stored_data[0], p0 = [guessMagLower, guessWidthLower, guessLocLower, guessMagMid, guessWidthMid, guessLocMid, guessMagHigh, guessWidthHigh, guessLocHigh])
    fit = ans[0][0]
    fit2 = ans[0][1]
    fit3 = ans[0][2]
    fit4 = ans[0][3]
    fit5 = ans[0][4]
    fit6 = ans[0][5]
    fit7 = ans[0][6]
    fit8 = ans[0][7]
    fit9 = ans[0][8]
    return plt.plot(stored_height[0], FitSim(stored_height[0], fit, fit2, fit3, fit4, fit5, fit6, fit7, fit8, fit9))

def CurveFitLightning(stored_data, stored_height, delay, asc_rate, z0, dz, payloads, startTime, dt, LightningTimes, LightningMag, LightningPos, LightningWidth, guessMagLower, guessWidthLower, guessLocLower, guessMagMid, guessWidthMid, guessLocMid, guessMagHigh, guessWidthHigh, guessLocHigh):
    Ltime = LightningTimes
    Lmag = LightningMag
    Lpos = LightningPos
    Lwidth = LightningWidth
    def FitSim(x, mag, width, loc, mag2, width2, loc2, mag3, width3, loc3, LightningMag, LightningPos, LightningWidth):
        def myI(zz):
            return mag*exp(-(zz-loc)**2/width) + mag2*exp(-(zz-loc2)**2/width2) + mag3*exp(-(zz-loc3)**2/width3)
        #print("sim")
        LightningTimes = Ltime
        LightningMag = Lmag
        LightningPos = Lpos
        LightningWidth = Lwidth
        data, a, b, c, d = RunSimulation(asc_rate, z0, dz, payloads, delay, startTime, dt, myI, [LightningTimes], [LightningMag], [LightningPos], [LightningWidth])
        #print("sim done")
        return data[0]
    ans= op.curve_fit(FitSim, stored_height[0], stored_data[0], p0 = [guessMagLower, guessWidthLower, guessLocLower, guessMagMid, guessWidthMid, guessLocMid, guessMagHigh, guessWidthHigh, guessLocHigh, LightningMag, LightningPos, LightningWidth])
    fit = ans[0][0]
    fit2 = ans[0][1]
    fit3 = ans[0][2]
    fit4 = ans[0][3]
    fit5 = ans[0][4]
    fit6 = ans[0][5]
    fit7 = ans[0][6]
    fit8 = ans[0][7]
    fit9 = ans[0][8]
    fit10 = ans[0][9]
    fit11 = ans[0][10]
    fit12 = ans[0][11]
    #print(fit, fit2, fit3, fit4, fit5, fit6, fit7, fit8, fit9, fit10, fit11, fit12)
    return plt.plot(stored_height[0], FitSim(stored_height[0], fit, fit2, fit3, fit4, fit5, fit6, fit7, fit8, fit9, fit10, fit11, fit12))
    