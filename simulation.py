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
                #print("lightning")
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
    #print(stored_data)
    return(stored_data, stored_height, stored_time, E_field, z)

def RunAnimation(stored_data, stored_height, E_field, z):
    fig = plt.figure()
    ax = plt.axes(xlim=(-140000,45000), ylim=(stored_height[0][0],11))
    line, = ax.plot([], [], lw=2, label = "Electric Field")
    line2, = ax.plot([], [], lw=2, label = "Payload Data")
    #line3, = ax.plot([], [], lw=2)
    #plot(threshold, z[1:])
    #plot(-threshold, z[1:])
    plt.xlabel("E (V/m)")
    plt.ylabel("Altitude (km)")
    plt.close()
    plt.legend(loc="upper left")
    def init():
        line.set_data([], [])
        line.set_color("orange")
        line2.set_data([], [])
        line2.set_color("blue")
        #line3.set_data([], [])
        return line,
    def animate(i):
        line.set_data(E_field[i],z[1:])
        line2.set_data(stored_data[0][:i],stored_height[0][:i])
        #line3.set_data(stored_data[1][:i],stored_height[1][:i])
        #to get this to work, i might be able to make it so payloads launched
        #after the first payload have an altitude of zero until their launch time. not sure how to get this to work however
        return line,
    return FuncAnimation(fig, animate, init_func=init, frames=260, interval=75, blit=True)
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

def InferCurrent(stored_data, stored_height, stored_time, delay, dt): 
    width = 0
    e0 = 8.854187817e-12
    df = DataFrame(list(zip(stored_data[0], stored_time[0], stored_height[0])), columns = ['e', 't', 'z'])
    df2 = DataFrame(list(zip(stored_data[1], stored_time[1], stored_height[1])), columns = ['e', 't', 'z'])
    gaussianPeak = argmax(abs(df['e']))
    gaussianPeakE = - max(abs(df['e']))
    gaussianPeakE2 = - max(abs(df2['e']))
    centerGaussian = df.iloc[gaussianPeak]['z']
    dE = gaussianPeakE2 - gaussianPeakE
    dT = df2.iloc[gaussianPeak]['t']+delay - df.iloc[gaussianPeak]['t'] 
    dEdT = dE/dT
    magnitude = -e0 * dEdT
    edge = gaussianPeakE - gaussianPeakE * .606 #.9545#.6827 #99.7% of normally distributed variables fall within 3 standard deviations of the mean
    for i in range(len(df)):
        #print(df2.iloc[i]['e'], df2.iloc[i]['z'])
        if (abs(df.iloc[i]['e']) < abs(edge)) & (df.iloc[i]['z'] < centerGaussian):
            if (abs(df.iloc[i+1]['e']) > abs(edge)):
                distance = centerGaussian-df.iloc[i]['z']
                width = distance
    #print(magnitude * exp(('z' - centerGaussian)**2 / (2 * width**2))
    print(centerGaussian, dEdT, magnitude, width, gaussianPeak)
    return(centerGaussian, magnitude, width)

def Infer2Currents(stored_data, stored_height, stored_time, delay, dt): 
    width = 0
    e0 = 8.854187817e-12
    df = DataFrame(list(zip(stored_data[0], stored_time[0], stored_height[0])), columns = ['e', 't', 'z'])
    df2 = DataFrame(list(zip(stored_data[1], stored_time[1], stored_height[1])), columns = ['e', 't', 'z'])
    gaussianPeakNegative = argmax(abs(df['e']))
    gaussianPeakPositive = argmax(df['e'])
    gaussianPeakNegativeE = - max(abs(df['e']))
    gaussianPeakNegativeE2 = - max(abs(df2['e']))
    gaussianPeakPositiveE = max(df['e'])
    gaussianPeakPositiveE2 = max(df2['e'])
    centerGaussianNegative = df.iloc[gaussianPeakNegative]['z']
    centerGaussianPositive = df.iloc[gaussianPeakPositive]['z']
    dE = gaussianPeakNegativeE2 - gaussianPeakNegativeE
    dEPositive = gaussianPeakPositiveE2 - gaussianPeakPositiveE
    dT = df2.iloc[gaussianPeakNegative]['t']+delay - df.iloc[gaussianPeakNegative]['t'] 
    dTPositive = df2.iloc[gaussianPeakPositive]['t']+delay - df.iloc[gaussianPeakPositive]['t']
    dEdT = dE/dT
    dEdTPositive = dEPositive / dTPositive
    magnitude = -e0 * dEdT
    magnitudePositive = -e0 * dEdTPositive
    edge = gaussianPeakNegativeE - gaussianPeakNegativeE * .6827 #.9545#.6827 #99.7% of normally distributed variables fall within 3 standard deviations of the mean
    edgePositive = gaussianPeakPositiveE - gaussianPeakPositiveE * .6827
    for i in range(len(df)):
        if (df.iloc[i]['e'] > edge) & (df.iloc[i]['z'] < centerGaussianNegative):
            if (df.iloc[i+1]['e'] < edge):
                distance = centerGaussianNegative-df.iloc[i]['z']
                width = distance
    for i in range(len(df)):
        if (df.iloc[i]['e'] < edgePositive) & (df.iloc[i]['z'] > centerGaussianPositive):
            if (df.iloc[i-1]['e'] > edgePositive):
                distancePositive = df.iloc[i]['z'] - centerGaussianPositive
                widthPositive = distancePositive
    print(centerGaussianNegative, magnitude, width, centerGaussianPositive, magnitudePositive, widthPositive)
    return(centerGaussianNegative, magnitude, width, centerGaussianPositive, magnitudePositive, widthPositive)



def CurrentAltitude(stored_data, stored_height, stored_time, delay, dt, asc_rate, LightningTimes):
    e0 = 8.854187817e-12
    I = []
    #df = DataFrame(list(zip(stored_data[0], stored_data[1], stored_time[0], stored_time[1], stored_height[0], stored_height[1])), columns = ['e', 'e2', 't', 't2', 'z', 'z2'])
    df1 = DataFrame(list(zip(stored_data[0], stored_time[0], stored_height[0])), columns = ['e', 't', 'z'])
    df2 = DataFrame(list(zip(stored_data[1], stored_time[1], stored_height[1])), columns = ['e', 't', 'z'])
    for i in range(len(df2)):
        pos1 = 0
        pos2 = 0
        current1 = 0
        current2 = 0
        for j in range(len(LightningTimes)):
            for k in range(len(df1)):
                if(df1.iloc[k]['t'] > LightningTimes[j]):
                    pos1 = df1.iloc[k]['z']
                    current1 = k-1
                    break
            for l in range(len(df2)):
                if(df2.iloc[l]['t'] > LightningTimes[j]):
                    pos2 = df2.iloc[l]['z']
                    current2 = l-1
                    break
        if((df2.iloc[i]['z'] < pos2) | (df2.iloc[i]['z'] >= pos1)):
            I.append(-e0 * (df2.iloc[i]['e'] - df1.iloc[i]['e']) / delay)
        else:
            I.append(0)
        #if((df2.iloc[i]['z'] >= pos2) & (df2.iloc[i]['z'] < pos1)):
            #I.append(0)
        #else:
            #I.append(-e0 * (df2.iloc[i]['e'] - df1.iloc[i]['e']) / delay)
    z = linspace(0, len(df1)*asc_rate*dt, len(df2))
    plt.plot(z, I, label = "Inferred Current")
    print(len(df1))
    #return df1, df2


def CurveFit(stored_data, stored_height, delay, asc_rate, z0, dz, payloads, startTime, dt, LightningTimes, LightningMag, LightningPos, LightningWidth, guessMagLower, guessWidthLower, guessLocLower, guessMagMid, guessWidthMid, guessLocMid, guessMagHigh, guessWidthHigh, guessLocHigh):
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
    print(fit, fit2, fit3, fit4, fit5, fit6, fit7, fit8, fit9, fit10, fit11, fit12)
    return plt.plot(stored_height[0], FitSim(stored_height[0], fit, fit2, fit3, fit4, fit5, fit6, fit7, fit8, fit9, fit10, fit11, fit12))
    