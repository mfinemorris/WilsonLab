# -*- coding: utf-8 -*-
"""
Model of somato-dendritic burster
From Toporikova and Butera, J. Comp neuroscience, DOI 10.1007/s10827-010-0274-z
Model with following parameters reproduce Figure 1B of the paper above
setting gcan=0, convert model to Pruvis and Butera 2007
units: V = mV; Kca = uM; Cm = uF/cm^2; g = uS/cm^2; phi = 1/s;
Parameter fi is normalized with respect to intracellular volume, so it's units fi=1/pL; 

Modified by Morgan Fine-Morris, October 2014


"""
import numpy as np
import matplotlib.pyplot as plt
import math
import scipy
from scipy.integrate import odeint
#
# INITIAL CONDITIONS
Vs0=-54
ns0=0.004
hs0=0.33
Vd0=-50
C0=0.03
l0=0.93

# variables
tEnd = 10000.
dt=0.3
tLen = int(tEnd/dt)
# PARAMETERS

# conductance of CaN current
gcan=1.

# conductances of persistent sodium current
gnaps=2.

#IP3 concentration
I=1

#Total Ca
Ct=1.25

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
sl=6
sn=-4
smp=-6
sh=5
taunb=10
tauhb=10000
gls=2.3
vleaks=-57
Iaps=0

# The ER parameters
fi=0.0001	
Vi=4
LL=0.37
P=31000
Ki=1.0
Ka=0.4
Ve=400
Ke=0.2
A=0.005 
Kd=0.4
sigma=0.185

# Other parameters 
Cmd=5
gc=1
k=0.3
Kcan=0.74
ncan=0.97
def df(y,t): 
    Vs=y[0]
    ns=y[1]
    hs=y[2]
    Vd=y[3]
    C=y[4]
    l=y[5]
    
    #Calculate ER Ca

    Ce = (Ct - C)/sigma
    
    # SOMATIC FUNCTIONS
    minfs=1/(1+math.exp((Vs-vm) /sm))
    ninfs=1/(1+math.exp((Vs-vn) /sn))
    minfps=1/(1+math.exp((Vs-vmp)/smp))
    hinfs =1/(1+math.exp((Vs-vh) /sh))
    
    tauns=taunb/math.cosh((Vs-vn)/(2*sn))
    tauhs=tauhb/math.cosh((Vs-vh)/(2*sh))
    
    I_nas=gna*math.pow(minfs,3)*(1-ns)*(Vs-vna)
    I_ks=gk*math.pow(ns,4)*(Vs-Vk)
    I_naps=gnaps*minfps*hs*(Vs-vna)
    
    I_ls =gls*(Vs-vleaks)
    I_ds =gc*(Vd-Vs)/k
    
    #DENDRITIC FUNCTIONS
    
    J_ER_in=(LL + P*math.pow(((I*C*l)/((I+Ki)*(C+Ka))),3))*(Ce - C)
    J_ER_out=Ve*math.pow(C,2)/(math.pow(Ke,2)+math.pow(C,2))
    
    caninf =1/(1+math.pow((Kcan/C),ncan))
    I_can=gcan*caninf*(Vd-vna)
    I_sd =gc*(Vs-Vd)/(1-k)

    # SOMATIC EQUATIONS

    dVs= (-I_ks - I_nas-I_naps-I_ls-Iaps-I_sd)/Cms
    dns= (ninfs-ns)/tauns
    dhs= (hinfs-hs)/tauhs

    # DENDRITIC EQUATIONS

    dVd= (-I_can-I_ds)/Cmd
    dC = fi/Vi*( J_ER_in- J_ER_out)
    dl = A*( Kd - (C + Kd)*l )


    
    return [dVs,dns,dhs,dVd,dC,dl]

y0 = [Vs0,ns0,hs0,Vd0,C0,l0]
t = scipy.linspace(0,tEnd,tLen)
soln = odeint(df, y0, t)
Vs = soln[:, 0]
ns=soln[:, 0]
hs=soln[:, 0]
Vd=soln[:, 0]
C0=soln[:, 0]
l0=soln[:, 0]
 
plt.figure()
plt.plot(t/1000.,Vs,'k')
plt.ylabel('V (mV)', fontsize=16)
plt.xlabel('time (sec)', fontsize=16)
plt.ylim(-60,10)
#plt.savefig('ca_urst.png')
plt.show()

