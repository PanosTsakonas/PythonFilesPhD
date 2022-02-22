import numpy as np
import sympy
from numpy import pi
from sympy import cos,sin
from matplotlib.pyplot import plot,show

from scipy.integrate import RK45
import sympy.physics.mechanics as mech
import matplotlib.pyplot as plt
mech.init_vprinting()

#Create symbolic arguments for the system
K,B,m1,m2,m3,l1,l2,l3=sympy.symbols('K,B,m1,m2,m3,l1,l2,l3',positive=True)
g=9.81
t=sympy.symbols('t',positive=True)

th1=sympy.Function('th1')(t)
th2=sympy.Function("th2")(t)
th3=sympy.Function('th3')(t)

#Creating the polar coordinate system
x1=sympy.Rational(1,2)*m1*l1*cos(th1)
y1=sympy.Rational(1,2)*m1*l1*sin(th1)
x2=sympy.Rational(1,2)*m2*l2*cos(th2)+l1*m1*cos(th1)
y2=sympy.Rational(1,2)*m2*l2*sin(th2)+l1*m1*sin(th1)
x3=sympy.Rational(1,2)*m3*l3*cos(th3)+l1*m1*cos(th1)+l2*m2*cos(th2)
y3=sympy.Rational(1,2)*m3*l3*sin(th3)+l1*m1*sin(th1)+l2*m2*sin(th2)

#Calculating the moments of inertia at the center of mass for each segment
I1=sympy.Rational(1,12)*m1*l1**2
I2=sympy.Rational(1,12)*m2*l2**2
I3=sympy.Rational(1,12)*m3*l3**2

#Calculating Kinetic energy of the system
Ktrans=sympy.Rational(1,2)*(m1*(x1.diff(t)**2+y1.diff(t)**2)+m2*(x2.diff(t)**2+y2.diff(t)**2)+m3*(x3.diff(t)**2+y3.diff(t)**2))
Krot=sympy.Rational(1,2)*(I1*th1.diff(t)**2+I2*th2.diff(t)**2+I3*th3.diff(t)**2)
Kin=Krot+Ktrans



#Calculating the Potential energy
Vgrav=-(m1*g*y1+m2*g*y2+m3*g*y3)
Vspring=sympy.Rational(1,2)*K*(th1**2+th2**2+th3**2)
V=Vgrav+Vspring

#Lagrangian Equation

L=Kin-V

#Rayleigh Dissipation Function
R=sympy.Rational(1,2)*B*(th1.diff(t)**2+th2.diff(t)**2+th3.diff(t)**2)

#Euler-Lagrange Equations

dt_th1d=L.diff(th1.diff(t)).diff(t)
dt_th2d=L.diff(th2.diff(t)).diff(t)
dt_th3d=L.diff(th3.diff(t)).diff(t)

dL_th1=L.diff(th1)
dL_th2=L.diff(th2)
dL_th3=L.diff(th3)

dR_th1d=R.diff(th1.diff(t))
dR_th2d=R.diff(th2.diff(t))
dR_th3d=R.diff(th3.diff(t))

#Grouping terms together
TH1=dt_th1d-dL_th1+dR_th1d
TH2=dt_th2d-dL_th2+dR_th2d
TH3=dt_th3d-dL_th3+dR_th3d

#Get the length of individual segments

L1=0.031
L2=0.024
L3=0.026
Ra=0.012
b=0.1
k=10

dens=1.1*10**(-3)*10**6
M1=pi*Ra**2*L1*dens
M2=pi*Ra**2*L2*dens
M3=pi*Ra**2*L3*dens
#System inputs
Q_1=0
Q_2=0
Q_3=0

#Solve the system
#accel_th1=sympy.solve(sympy.Eq(TH1,Q_1),th1.diff(t,2))
#accel_th2=sympy.solve(sympy.Eq(TH2,Q_2),th2.diff(t,2))
#accel_th3=sympy.solve(sympy.Eq(TH3,Q_3),th3.diff(t,2))


TH1=TH1.subs(l1,L1).subs(l2,L2).subs(l3,L3).subs(m1,M1).subs(m2,M2).subs(m3,M3).subs(K,k).subs(B,b)
TH2=TH2.subs(l1,L1).subs(l2,L2).subs(l3,L3).subs(m1,M1).subs(m2,M2).subs(m3,M3).subs(K,k).subs(B,b)
TH3=TH3.subs(l1,L1).subs(l2,L2).subs(l3,L3).subs(m1,M1).subs(m2,M2).subs(m3,M3).subs(K,k).subs(B,b)


def eqn(t):
    return {TH1==Q_1,TH2==Q_2, TH3==Q_3}
init_cond=np.array([0,0,0,0,0,0])
t_span=np.array([0,1])

sol=RK45(eqn(t),t_span[0],[init_cond[0],init_cond[1]],t_span[1])
plot(sol)
show()
