import numpy as np
from numpy import zeros, sin, cos,pi
import os
from scipy import optimize
from matplotlib.pyplot import plot,show
import pandas as pd
#import the data from the excel file
#a=input("Participant number:")
a="1"
s='C:/Users/u1857308/OneDrive - University of Warwick/PhD/Gait Lab Data//Digit Weights/P_'+a+'.xlsx'
f=pd.read_excel(s)
n=(5,3)
L=np.zeros(n)
R=np.zeros(n)
M=zeros(n)
I=zeros(n)
Digit=["Ldip (mm)","Rdip (mm)","Lpip (mm)","Rpip (mm)", "Lmcp (mm)","Rmcp (mm)", "Mdip (kg)", "Mpip (kg)","Mmcp (kg)"]

#Leave this you generate the required mass, length and radius matrices for DIP, PIP, MCP joints
D=pd.DataFrame([f[Digit[0]], f[Digit[1]],f[Digit[2]],f[Digit[3]],f[Digit[4]],f[Digit[5]],f[Digit[6]],f[Digit[7]],f[Digit[8]]])
m=np.concatenate(D)
n1=0
print(len(f[Digit[0]]))
P=zeros((len(f[Digit[0]]),9))
for i in range (0,5):
    for j in range (0,9):
        P[i,j]=m[j+n1]
    n1=j+1+n1

for i in range (0,5):
    for j in range(0,3):
        L[i,2-j]=P[i,2*j]*10**-3
        R[i,2-j]=P[i,2*j+1]*10**-3
        M[i,2-j]=P[i,6+j]

#Create the moment of inertia variables from the data collected
Icmc=M[0,0]*(R[0,0]**2 /4+L[0,0]**2 /3)+(M[0,1]+M[0,2])*L[0,0]**2
Imp=M[0,1]*(R[0,1]**2 /4+L[0,1]**2 /3)+M[0,2]*L[0,1]**2
Iip=M[0,2]*(R[0,2]**2/4+L[0,2]**2/3)
Imcp2=M[1,0]*(R[1,0]**2/4+L[1,0]**2/3)+(M[1,1]+M[1,2])*L[1,0]**2
Ipip2=M[1,1]*(R[1,1]**2/4+L[1,1]**2/3)+M[1,2]*L[1,1]**2
Idip2=M[1,2]*(R[1,2]**2/4+L[1,2]**2/3)
Imcp3=M[2,0]*(R[2,0]**2/4+L[2,0]**2/3)+(M[2,1]+M[2,2])*L[2,0]**2
Ipip3=M[2,1]*(R[2,1]**2/4+L[2,1]**2/3)+M[2,2]*L[2,1]**2
Idip3=M[2,2]*(R[2,2]**2/4+L[2,2]**2/3)
Imcp4=M[3,0]*(R[3,0]**2/4+L[3,0]**2/3)+(M[3,1]+M[3,2])*L[3,0]**2
Ipip4=M[3,1]*(R[3,1]**2/4+L[3,1]**2/3)+M[3,2]*L[3,1]**2
Idip4=M[3,2]*(R[3,2]**2/4+L[3,2]**2/3)
Imcp5=M[4,0]*(R[4,0]**2/4+L[4,0]**2/3)+(M[4,1]+M[4,2])*L[4,0]**2
Ipip5=M[4,1]*(R[4,1]**2/4+L[4,1]**2/3)+M[4,2]*L[4,1]**2
Idip5=M[4,2]*(R[4,2]**2/4+L[4,2]**2/3)

#Now define the damper and spring constants for each segment
Kcmca=5
Bcmca=0.1
Kcmc=10
Bcmc=0.1
Kmp=10
Bmp=0.1
Kip=10
Bip=0.1
Kmcp2a=0.5
Bmcp2a=0.01
Kmcp2=10
Bmcp2=0.1
Kpip2=10
Bpip2=0.1
Kdip2=10
Bdip2=0.1
Kmcp3a=0.
Bmcp3a=0.01
Kmcp3=10
Bmcp3=0.1
Kpip3=10
Bpip3=0.1
Kdip3=10
Bdip3=0.1
Kmcp4a=0.56
Bmcp4a=0.06
Kmcp4=10
Bmcp4=0.1
Kpip4=10
Bpip4=0.1
Kdip4=10
Bdip4=0.1
Kmcp5a=0.96
Bmcp5a=0.03
Kmcp5=10
Bmcp5=0.1
Kpip5=10
Bpip5=0.1
Kdip5=10
Bdip5=0.1

#Now define the equilibrium agles for all segments. Give values in degrees
theqcmca=0
theqcmc=0
theqmp=0
theqip=0
theqmcp2a=0
theqmcp2=0
theqpip2=10*pi/180
theqdip2=0
theqmcp3a=0
theqmcp3=0
theqpip3=0
theqdip3=0
theqmcp4a=0
theqmcp4=0
theqpip4=0
theqdip4=0
theqmcp5a=0
theqmcp5=0
theqpip5=0
theqdip5=0


theq=[theqcmca,theqcmc,theqmp,theqip,theqmcp2a,theqmcp2,theqpip2,theqdip2,theqmcp3a,theqmcp3,theqpip3,theqdip3,theqmcp4a,theqmcp4,theqpip4,theqdip4,theqmcp5a,theqmcp5,theqpip5,theqdip5]
K=[Kcmca,Kcmc,Kmp,Kip,Kmcp2a,Kmcp2,Kpip2,Kdip2,Kmcp3a,Kmcp3,Kpip3,Kdip3,Kmcp4a,Kmcp4,Kpip4,Kdip4,Kmcp5a,Kmcp5,Kpip5,Kdip5]
B=[Bcmca,Bcmc,Bmp,Bip,Bmcp2a,Bmcp2,Bpip2,Bdip2,Bmcp3a,Bmcp3,Bpip3,Bdip3,Bmcp4a,Bmcp4,Bpip4,Bdip4,Bmcp5a,Bmcp5,Bpip5,Bdip5]

t=[1,2,3,4,5,6,7,8]
#Once everything is defined, the the next step will be to pre-define the size of the
#muscle level activation matrices

aFDP=zeros((1,len(t)))
aFDS=zeros((1,len(t)))
aEDC=zeros((1,len(t)))
aEIP=zeros((1,len(t)))
aEDM=zeros((1,len(t)))
aLUMI=zeros((1,len(t)))
aLUMM=zeros((1,len(t)))
aLUMR=zeros((1,len(t)))
aLUML=zeros((1,len(t)))
aPI3=zeros((1,len(t)))
aPI2=zeros((1,len(t)))
aPI1=zeros((1,len(t)))
aDI2=zeros((1,len(t)))
aDI3=zeros((1,len(t)))
aDI4=zeros((1,len(t)))
aDIMC1=zeros((1,len(t)))
aDIMC2=zeros((1,len(t)))
aADP=zeros((1,len(t)))
aAPB=zeros((1,len(t)))
aAPL=zeros((1,len(t)))
aEPB=zeros((1,len(t)))
aEPL=zeros((1,len(t)))
aFPB=zeros((1,len(t)))
aFPL=zeros((1,len(t)))
aOPP=zeros((1,len(t)))

#Total 26 muscles are used in this model

lb=zeros((1,25));
ub=np.ones((1,25))
init=np.random.rand(1,25)

