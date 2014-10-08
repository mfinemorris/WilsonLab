ti'''
In the following, m, h, and n are gating variables controlling the activation or inactivation of certain channels.
M is approximated as starting instantaneously, so we never need a dm/dt equation-we can just use m_inf in place of m.
Below are three strings containing different models, different in only the following way: 
one has both I_NaPh and I_NaP, one has the first only, and one has the second only. 
I've tested all three and found that only the equation with I_NaPh works properly. The only current issue for the code
is that the voltages are too small. Otherwise, the shape of the curve seems correct.
'''

from brian import *
from brian.library.ionic_currents import *
import time

reinit_default_clock()
print defaultclock.dt
#number of seconds to run simulation (include units)
runtime = 5.0*second


# temperature dependence, borrowed from Cochlear neuron model of Rothman & Manis in Brian Examples @
# http://www.briansimulator.org/docs/examples-frompapers_Rothman_Manis_2003.html
# I've used q10 term in gateing variable eqns of eqs_noNaP.
celsius = 36.
q10 = 1.0#3.**((celsius-22)/10.)#commented out bc q10 is related to temperature change not absolute temp.

#The nominal params for leak current are gL=2.8 nS and EL=-65 mV.
EL = -60*mvolt
EK = -85.0*mvolt
ENa = 50.0*mvolt
C = 21 * pF
gL = 2.8*nsiemens
gK = 11.2*nsiemens
gNa = 28*nsiemens
gNaP = 2.8*nsiemens

#---Model 1 of butera99a
#---time constants are consistent for all currents in Model 1 of butera99a
tauMaxH = 10000*ms
tauMaxN = 10*ms

#constants for K and Na
thetaN = -29*mvolt
thetaM = -34*mvolt
sigmaN = -4*mvolt
sigmaM = -5*mvolt

#constants for NaP and Nap-h
thetaH = -48*mvolt
sigmaH = 6*mvolt
thetaM_nap = -40*mvolt
sigmaM_nap = -6*mvolt

#eqns for Na and K
mInf = lambda v: (1+exp((v-thetaM)/sigmaM))**-1
nInf = lambda v: (1+exp((v-thetaN)/sigmaN))**-1
tauN = lambda v: tauMaxN/cosh((v-thetaN)/(2*sigmaN))

#eqns for NaP and NaP-h
hInf = lambda v: (1+exp((v-thetaH)/sigmaH))**-1
mInf_nap = lambda v: (1+exp((v-thetaM_nap)/sigmaM_nap))**-1
tauH = lambda v: tauMaxH/cosh((v-thetaH)/(2*sigmaH))

#no nap - this one seems to be correct
eqs_noNaP = '''
dV/dt = -((gNa*(1-n)*(V-ENa)*mInf_1**3)+ (gNaP*(V-ENa)*h*mInf_2) + \
(gK*(V-EK)*n**4) + gL*(V-EL))/C:mvolt
#
#Sodium and Potassium Activation
mInf_1 = mInf(V):1
n_inf = nInf(V):1
tau_n=tauN(V):ms
dn/dt =  q10*(n_inf - n)/tau_n:1
#
#sodium inactivation (nap or naph)
mInf_2 = mInf_nap(V):1
h_inf = hInf(V):1
tau_h=tauH(V):ms
dh/dt = q10*(h_inf - h)/tau_h:1
I:amp
'''

start = time.clock()

#neurons=NeuronGroup(1, eqs_noNaP, clock = defaultclock)
neurons=NeuronGroup(1, eqs_noNaP, method = 'RK', clock = defaultclock)

neurons.V = EL
neurons.h = hInf(EL)
neurons.n = nInf(EL)
M = StateMonitor(neurons, 'V', record=True)

run(runtime)
print "duration = ",(time.clock() - start)
print "time ", time.strftime("%H:%M:%S")

p=plot(M.times, M[0]/mV)
xlabel('Time (ms)')
ylabel('Membrane potential (mV)')
plt.show()
