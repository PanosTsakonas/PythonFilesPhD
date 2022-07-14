import matplotlib.pyplot as plt
import numpy as np
from numpy import pi,cos,sin,exp
from scipy.signal import butter, filtfilt
from matplotlib.pyplot import plot,show,legend,title,xlabel,ylabel
from scipy.optimize import curve_fit as cf
import pandas as pd

b=["thumb","index","middle","ring","little"]
a=int(input("What digit are you working on: "))
#m="C:/Users/panos/Desktop/"+b[a-1]+"_abd_spring.csv"
m="C:/Users/panos/Desktop/little_abd_spring.csv"
print(m)
th=pd.read_csv(m,delimiter=",")

wn=5
print(th)
