#This code calculates the sensitivity analysis of the IBK model using the state space representation
import numpy as np
pi=np.pi
import scipy.signal as signal
from scipy.signal import StateSpace as SS
import matplotlib.pyplot as plt
from SALib.sample import saltelli as sl
from SALib.analyze import sobol as sb


theq= 20*pi/180
th0= 50*pi/180
K= 0.03963
B1= 0.0011936
I1= 3.45125E-05


A= [[0, 1],
   [-K/I1, -B1/I1]]

B= [[0],
   [K/I1]]
B2= [[0], [1/I1]]
C= [[1, 0]]

D=0

#Creates the state space representation for the free response of IBK model
sys = SS(A, B, C, D)
sys1 = SS(A, B2, C, D)
init = 0
last = 0.5
t = np.linspace(init,last, 2000)
N = len(t)

#initial conditions
x0 = [th0, 0]
u = theq*np.ones(N)
u1 = np.zeros(N)+K*theq*np.ones(N)
t, y, x = signal.lsim(sys, u, t, x0)
t, y1, x1 =signal.lsim(sys1, u1, t, x0)

plt.figure(1)
plt.plot(t, y*180/pi)
plt.title("Free response of IBK state space model")
plt.xlabel("Time (s)")
plt.ylabel("angle (degrees)")
plt.show()

plt.figure(2)
plt.plot(t, (y-y1)*180/pi)
plt.title("Difference between the two methods of state space")
plt.xlabel("Time (s)")
plt.ylabel("angle (degrees)")
plt.show()

del A, B, C, D, B2, x, x1,t

problem={
   'num_vars':3,
   'names':['k','b','i1'],
   'bounds':[[K-4*0.0036522, K+4*0.0036522],
             [B1-4*0.00010063, B1+4*0.00010063],
             [I1-I1*20/100, I1+I1*20/100]]
}

param_values = sl.sample(problem, 1024)

def ibk(x,k,b,i1):
   w=np.sqrt(4*k*i1-b**2)/(2*i1)
   return theq+(x0[0]-theq)*np.exp(-b*x/(2*i1))*((b/(np.sqrt(4*k*i1-b**2)))*np.sin(w*x)+np.cos(w*x))

x = np.linspace(init,last, 2000)
y=np.array([ibk(x,*params) for params in param_values])

# analyse
sobol_indices = [sb.analyze(problem, Y) for Y in y.T]


S1s = np.array([s['S1'] for s in sobol_indices])

fig = plt.figure(3,figsize=(10, 6), constrained_layout=True)
gs = fig.add_gridspec(3, 3)

ax0 = fig.add_subplot(gs[:, 0])
ax1 = fig.add_subplot(gs[0, 1])
ax2 = fig.add_subplot(gs[1, 1])
ax3 = fig.add_subplot(gs[2, 1])
for i, ax in enumerate([ax1, ax2, ax3]):
    ax.plot(x, S1s[:, i],
            label=r'S1$_\mathregular{{{}}}$'.format(problem["names"][i]),
            color='black')
    ax.set_xlabel("Time (s)")
    ax.set_ylabel("First-order Sobol index")

    ax.set_ylim(0, 1.04)

    ax.yaxis.set_label_position("right")
    ax.yaxis.tick_right()

    ax.legend(loc='upper right')

ax0.plot(x, np.mean(y, axis=0), label="Mean", color='black')

# in percent
prediction_interval = 95

ax0.fill_between(x,
                 np.percentile(y, 50 - prediction_interval/2., axis=0),
                 np.percentile(y, 50 + prediction_interval/2., axis=0),
                 alpha=0.5, color='black',
                 label=f"{prediction_interval} % prediction interval")

ax0.set_xlabel("Time (s)")
ax0.set_ylabel("angle (rad)")
ax0.legend(title=r"Underdamped IBK model", loc='upper center')._legend_box.align = "left"

plt.show()