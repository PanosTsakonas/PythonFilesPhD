import matplotlib.pyplot as plt
import numpy as np
from numpy import pi,cos,sin,exp
from scipy.signal import butter, filtfilt
from matplotlib.pyplot import plot,show,legend,title,xlabel,ylabel
from scipy.optimize import curve_fit as cf


data=open("C:/Users/u1857308/OneDrive - University of Warwick/PhD/MATLAB Code/Unfilt.csv","r")#importing the  unfiltered data from MATLAB
f=np.loadtxt(data,delimiter=",")#Loading the csv file
t=np.zeros(len(f))#time matrix of length f
th=np.zeros(len(f))#angle matrix of length f
data.close()#closing the MATLAB file
#A python array starts from 0 and let f[i,j] be a numpy array then the i corresponds to the row and the j to the collumn
for i in range (0,len(f)):
    t[i]=f[i,0]
    th[i]=f[i,1]*pi/180#Angles in radians

#Import the moment of inertia
m=open("//doozer.ads.warwick.ac.uk/User64/u/u1857308/Documents/flexion_moments.csv","r")
I1=np.loadtxt(m,delimiter=",")
m.close()
Imcp=I1[0,0]
Ipip=I1[0,1]
Idip=I1[0,2]

#Butterworth 4th order Low pass filter with 20 Hz cut-off
b1,a1=butter(4,20/200,btype="low")
th=filtfilt(b1,a1,th)

#The fitting process#

#Define the fitting function
def fun(x,a,b,c,d,w):
    return a+c*exp(-b*x)*(d*sin(w*x)+cos(w*x))

p0=(0.4,25,0.2,1,56)
#Curve fitting
fit,GoF=cf(fun,t,th,absolute_sigma=True,method="lm")

a,b,c,d,w=fit
print(fit)
plt.figure(1)
plot(t,th)
plot(t,fun(t,a,b,c,d,w))
title("DIP joint data")
xlabel("Time (s)")
ylabel("Angle (rad)")


plt.figure(2)
plot(t,np.abs(th-fun(t,a,b,c,d,w))*180/pi)
xlabel("Time (s)")
ylabel("Angle (rad)")
title("Residuals")
show()

B=2*Idip*b
print(B)
K=((w*2*Idip)**2+B**2)/(4*Idip)
print(K)
