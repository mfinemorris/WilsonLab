# -*- coding: utf-8 -*-
"""
The following modified by Morgan Fine-Morris, October 2014. 
Dendritic component removed and variable values of somatic component adjusted to match model 1 of
Butera et al. , J. Neurophysiol, July 1999 Vol. 82 no. 1 382-397

Vleaks values from Butera paper are -60,-57.5,-54 (all mV)

Model of somato-dendritic burster
From Toporikova and Butera, J. Comp neuroscience, DOI 10.1007/s10827-010-0274-z
Model with following parameters reproduce Figure 1B of the paper above
setting gcan=0, convert model to Pruvis and Butera 2007
units: V = mV; Kca = uM; Cm = uF/cm^2; g = uS/cm^2; phi = 1/s;
Parameter fi is normalized with respect to intracellular volume, so it's units fi=1/pL; 
"""
import numpy as np
import matplotlib.pyplot as plt
import math
import scipy
from scipy.integrate import odeint

import time


#init V and vleaks
v = -60.0

# INITIAL CONDITIONS
Vs0=v
ns0=0.004
hs0=0.33

# variables
tEnd=20*1000.
dt=0.3
tLen=int(tEnd/dt)

# PARAMETERS

# conductances of persistent sodium current
gnaps=2.8

# Somatic Parameters (same as Butera and Pruvis 2007)
Cms=21

vna=50
Vk=-85

gk=11.2
gna=28

vm=-34
vn=-29
vmp=-40
vh=-48

sm=-5
sn=-4
smp=-6
sh=6

taunb=10
tauhb=10000

# leak current variables
gls=2.8
vleaks=v

Iaps=0

def df(y,t): 
    Vs=y[0]
    ns=y[1]
    hs=y[2]
    
    # SOMATIC FUNCTIONS
    minfs=1/(1+math.exp((Vs-vm) /sm))
    ninfs=1/(1+math.exp((Vs-vn) /sn))
    minfps=1/(1+math.exp((Vs-vmp)/smp))
    hinfs=1/(1+math.exp((Vs-vh) /sh))
    
    tauns=taunb/math.cosh((Vs-vn)/(2*sn))
    tauhs=tauhb/math.cosh((Vs-vh)/(2*sh))
    
    I_nas=gna*math.pow(minfs,3)*(1-ns)*(Vs-vna)
    I_ks=gk*math.pow(ns,4)*(Vs-Vk)
    I_naps=gnaps*minfps*hs*(Vs-vna)
    
    I_ls =gls*(Vs-vleaks)

    # SOMATIC EQUATIONS

    dVs= (-I_ks - I_nas-I_naps-I_ls-Iaps)/Cms
    dns= (ninfs-ns)/tauns
    dhs= (hinfs-hs)/tauhs
    
    return [dVs,dns,dhs]


#start timer
start = time.clock()
#solve eq
y0 = [Vs0,ns0,hs0]
t = scipy.linspace(0,tEnd,tLen)
soln = odeint(df, y0, t)
Vs = soln[:, 0]
ns=soln[:, 0]
hs=soln[:, 0]
#end timer
end = time.clock()
print "Duration: ", end - start, "seconds"
 

#plot voltage vs time
plt.figure()
plt.plot(t/1000.,Vs,'k')
print Vs
plt.ylabel('V (mV)', fontsize=16)
plt.xlabel('time (sec)', fontsize=16)
plt.ylim(-80,10)
#plt.savefig('ca_urst.png')
plt.show()

