from numpy import *
import matplotlib.pyplot as plt
import matplotlib
import pandas
from scipy import integrate
from time import clock
import os, time
set_printoptions(4)
import NMResources
import utility


#Define BRS and TB Models
class BRSModel(object):

    #### MEMBRANE CAPACITANCE (pF) ####
    Cms=21

    #### REVERSAL POTENTIALS (mV) ####
    # sodium potential (for normal sodium current and persistent sodium current)
    ena=50
    # potassium potential
    ek=-85
    # leak reversal potential
    eL=-60.0 # suggested value for vleak from Butera paper are -60,-57.5,-54
    # reversal potential of non-NMDA glutamatergic synaptic currents.
    #eSyn=0.0

    #### CONDUCTANCE (nS) ####
    # potassium current conductance
    gk=11.2
    # sodium current conductance
    gna=28
    # persistent sodium current conductance
    gnaps=2.8 #should this be 1.5 or 2.0 as in Toparikova and Butera 2011?
    # leak current conductance
    gL=2.8
    # tonic conductance
    #gtonic=0.0

    #### APPLIED CURRENT (pA) ####
    Iaps=0

    #### HALF (IN)ACTIVATION VOLTAGE (mV) ####
    vm=-34
    vn=-29
    vmp=-40
    vh=-48

    #### INVERSE SLOPE MULTIPLIER (mV) ####
    sm=-5
    sn=-4
    smp=-6
    sh=5

    #### TIME CONSTANT MAX VALUES (ms) #####
    taunb=10 #tau max for n
    tauhb=10000 #tau max for h
    
    # In milliseconds
    dt = .1
    
    Vs0=-60.0 #initial somatic memberane voltage
    ns0=0.004 #initial value of gating var n
    hs0=0.33 #initial value of gating var h
    
    #set here so that only the variables are in param_names
    _param_names = [i for i in dir() if "__" not in i]
    
    #let user know if parameter dictionary passed in to setVars is incomplete
    warn_for_missing_vars = False
        
    def __init__(self, param_dict={}):          
        self.setVars(param_dict)
    
    def setVars(self, param_dict, warn = warn_for_missing_vars):
        """
        Set model parameter values.
        param_dict should be a dictionary where the keys are names for the class variables, 
        and the values are whatever you want those variables set too.
        """
        for i in self._param_names:
            try:
                self.__dict__[i] = float(param_dict[i])
            except (KeyError, AttributeError) as e:
                #if is missing from var dictionary, let user know or not, depending on boolean val
                self.__dict__[i] = self.__class__.__dict__[i]
                if  warn:
                    print ('Param '+i+" is either missing or invalid. Using default.")
    
    def printVarsHelper(self, cols=1, exclude=[]):
        '''
        prints parameters with less formatting
        Cols is how many columns to print the output in, 
        exclude is a list of strings, which should be parameter keys to exclude from the output
        '''
        print "{:<9s}: {:<10} \t ".format('Model', self.__class__.__name__)
        
        i = 1        
        s = sorted(self._param_names)
        for key in self._param_names:
            if key in exclude:
                continue
            try:
                item = vars(self)[key]
            except (KeyError, AttributeError):
                item = self.__class__.__dict__[key]

            print "{:<9s}: {:<10} \t ".format(key, item),
            if i%cols == 0:
                print
            i +=1
    
    def printVars(self,cols = 3):
        '''
        Print the model parameters in an easy to read format. 
        Set cols to the number of columns in which you would like the variables to appear.
        '''
        
        line = "-----------------------"
        print line + "\n"+"Model Parameters: "+"\n"+line+"\n"
        print "Initial state: Vs = {}, ns = {}, hs = {} ".format(self.Vs0, self.ns0, self.hs0,)
        print
        self.printVarsHelper(3,['Vs0', 'ns0', 'hs0'])
        print "\n"+line + "\n"+"End Parameters "+"\n"+line+'\n'
    
    def get_parameter_keys(self):
        '''
        Returns a list of the model parameters that can be set with setVars. 
        You can use printVars to see a formatted list of parameter keys and values.
        '''
        return deepcopy(self._param_names)
    
    def model(self, y, t):
        '''
        Primarily for internal use by the simulate function. y is the cell state from the previous time step 
        (consisting of an array with Vs, n, and h) and t is the time at which this timestep occurs. 
        '''
        Vs = y[0]
        ns = y[1]
        hs = y[2]
        
        # SOMATIC FUNCTIONS
        minfs = 1/(1+exp((Vs-self.vm) /self.sm))
        ninfs = 1/(1+exp((Vs-self.vn) /self.sn))
        minfps = 1/(1+exp((Vs-self.vmp)/self.smp))
        hinfs = 1/(1+exp((Vs-self.vh) /self.sh))

        tauns = self.taunb/cosh((Vs-self.vn)/(2*self.sn))
        tauhs = self.tauhb/cosh((Vs-self.vh)/(2*self.sh))

        # CURRENT EXPRESSIONS
        # currents in Soma
        I_nas = self.gna*(minfs**3)*(1-ns)*(Vs-self.ena) #sodium current
        I_ks = self.gk*(ns**4)*(Vs-self.ek) #potassium current
        I_naps = self.gnaps*minfps*hs*(Vs-self.ena) #persistent sodium current
        I_L = self.gL*(Vs-self.eL) #pA
        #I_tonic = self.gtonic*(Vs-eSyn)

        #### DIFFERENTIAL EQUATIONS
        # SOMATIC EQUATIONS
        dVs = (-I_ks-I_nas-I_naps-I_L+self.Iaps)/self.Cms
        dns = (ninfs-ns)/tauns
        dhs = (hinfs-hs)/tauhs        
        
        return [dVs, dns, dhs]

    def jacobian(self,y,t):
        '''
        Primarily for internal use by the simulate function. y is the cell state from the previous time step 
        (consisting of an array with Vs, n, and h) and t is the time at which this timestep occurs. 
        '''
        Vs = y[0]
        n = y[1]
        h = y[2]
        
        j11 = (-h*self.gnaps/(exp((Vs - self.vmp)/self.smp) + 1) + h*self.gnaps*(Vs - self.ena)*exp((Vs - self.vmp)/self.smp)/(self.smp*(exp((Vs - self.vmp)/self.smp) + 1)**2) - n**4*self.gk - self.gL - self.gna*(-n + 1)/(exp((Vs - self.vm)/self.sm) + 1)**3 + 3*self.gna*(Vs - self.ena)*(-n + 1)*exp((Vs - self.vm)/self.sm)/(self.sm*(exp((Vs - self.vm)/self.sm) + 1)**4))/self.Cms 
        j12 =  (-4*n**3*self.gk*(Vs - self.ek) + self.gna*(Vs - self.ena)/(exp((Vs - self.vm)/self.sm) + 1)**3)/self.Cms 
        j13 = -self.gnaps*(Vs - self.ena)/(self.Cms*(exp((Vs - self.vmp)/self.smp) + 1))
        j21 = (-n + 1/(exp((Vs - self.vn)/self.sn) + 1))*sinh((Vs - self.vn)/(2*self.sn))/(2*self.sn*self.taunb) - exp((Vs - self.vn)/self.sn)*cosh((Vs - self.vn)/(2*self.sn))/(self.sn*self.taunb*(exp((Vs - self.vn)/self.sn) + 1)**2)
        j22 = -cosh((Vs - self.vn)/(2*self.sn))/self.taunb
        j23 = 0
        j31 = (-h + 1/(exp((Vs - self.vh)/self.sh) + 1))*sinh((Vs - self.vh)/(2*self.sh))/(2*self.sh*self.tauhb) - exp((Vs - self.vh)/self.sh)*cosh((Vs - self.vh)/(2*self.sh))/(self.sh*self.tauhb*(exp((Vs - self.vh)/self.sh) + 1)**2)
        j32 = 0
        j33 = -cosh((Vs - self.vh)/(2*self.sh))/self.tauhb
        
        return [[j11, j12, j13], [j21, j22, j23], [j31,j32,j33]] 
    
    def _autosave(self, autosave_dir, sim_data):
        utility.mkdir_p(autosave_dir)
        #filename
        temp = ['%s',self.__class__.__name__]+list(map(lambda x: str(x), time.localtime()[0:6]))
        file_name = os.path.join(autosave_dir, "-".join(temp)+'.csv')
        
        #make and save dataframe
        str_formatted_data = sim_data.to_csv()#file_name%'all_sim_traces')
        #also save all model info to txt file
        with open(file_name%'simulation_autosave', "w") as f:
            with utility.stdout_redirected(f):
                self.printVarsHelper()
                print '--------------------'
                print str_formatted_data

                
    def simulate(self, simulationTime, use_jacobian=True, autosave_dir = 'autosaved_sim_data'):
        '''
        Simulates the model for number of milliseconds indicated by simulationTime. 
        Returns two arrays, V and t: t is the array of time points at which the system was calculated
        and V is the membrane voltage found at each of the time points
        Create dir specified in autosave_dir (do not overwrite if it already exists),
        autosaves all traces and model parameters to uniquely named, time-stamped files in that dir.
        Note that a value for autosave_dir is already specified. To prevent autosaving set
        autosave_dir to '' or None.
        '''
        #make time array and initial state array for odeint function
        t = linspace(0,simulationTime, simulationTime/self.dt)
        initial_state = [self.Vs0, self.ns0, self.hs0]

        try:
            if use_jacobian:
                y = integrate.odeint(self.model, initial_state, t, Dfun=self.jacobian)
            else:
                y = integrate.odeint(self.model, initial_state, t)
        except KeyboardInterrupt:
            if autosave_dir:
                sim_data = pandas.DataFrame(data=y, index=t, columns=['Vs0', 'ns0', 'hs0'])
                self._autosave(autosave_dir, sim_data)
            raise
        except Exception as e:
                raise("Could not run simulation. Exceptions: %s"%e)
        '''
        try:
            if use_jacobian:
                y = integrate.odeint(self.model, initial_state, t, Dfun=self.jacobian)
            else:
                raise Exception() #to be caught by except, so that integration is tried w/out jacobian
        except Exception as e:
            try:
                print "Running simulation without jacobian.", e
                y = integrate.odeint(self.model, initial_state, t)
            except Exception as a:
                raise("Could not run simulation. Exceptions: %s"%a)
        '''
        #extract membrane voltage
        V = y.T[0] 

        #save data to filename based on date and time
        if autosave_dir:
            sim_data = pandas.DataFrame(data=y, index=t, columns=['Vs0', 'ns0', 'hs0'])
            self._autosave(autosave_dir, sim_data)

        return V, t  


class TBModel(BRSModel):
        
    #### MEMBRANE CAPACITANCE (pF) ####
    Cms=21.
    #dendrite membrane capacitance
    Cmd=5.

    #### REVERSAL POTENTIALS (mV) ####
    # sodium potential (for normal sodium current and persistent sodium current)
    ena=50.
    # potassium potential
    ek=-85.
    # leak reversal potential
    eL=-60.0 # suggested value for vleak from Butera paper are -60,-57.5,-54
    # reversal potential of non-NMDA glutamatergic synaptic currents.
    #eSyn=0.0
    
    #### CONDUCTANCE (nS) ####
    # potassium current conductance
    gk=11.2
    # sodium current conductance
    gna=28.
    # persistent sodium current conductance
    gnaps=2.8#1.5#2.8 #should this be 1.5 or 2.0 as in Toparikova and Butera 2011?
    # leak current conductance
    gL=2.3 #in Toporikova_Butera_2010_code.py
    # tonic conductance
    #gtonic=0.0
    # calcium channel conductance
    gcan = 1.5
    # gc (the conductance for the link terms)
    gc = 1.0

    #### APPLIED CURRENT (pA) ####
    Iaps=0

    #### HALF (IN)ACTIVATION VOLTAGE (mV) ####
    vm=-34.
    vn=-29.
    vmp=-40.
    vh=-48.

    #### INVERSE SLOPE MULTIPLIER (mV) ####
    sm=-5.
    sn=-4.
    smp=-6.
    sh=5.

    #### TIME CONSTANT MAX VALUES (ms) #####
    taunb=10. #tau max for n
    tauhb=10000. #tau max for h

    #### Constants for calculating Ca Flux: ER --> Cytosol ####
    IP=1. #IP3 concentration
    LL=0.37 #ER leak permeability
    P=31000. #maximum total permeability of IP3 channels
    Ki=1.0 #dissociation consts for IP3 receptor activity by IP3
    Ka=0.4 #dissociation consts for IP3 receptor activity by Ca

    #### Constants for calculating Ca Flux: Cytosol --> ER ####
    Ve=400. #Maximal SERCA pump rate
    Ke=0.2 #coefficient for SERCA pumps

    #### ER Ca CONCENTRATION ####
    Ct=1.25 #Total Ca
    sigma=0.185 #ratio of cytosolic to ER volume

    #The ER parameters
    fi=0.0001 #bound Ca concentration in cytosol
    Vi=4. #free Ca concentration in cytosol
    A=0.005 #scaling const.
    Kd=0.4 #dissociation constant for IP3 receptor inactivation by Ca

    #### Calcium current activation ####
    Kcan=0.74 # microM 
    ncan=0.97

    #### Ratio of somatic to total area ####
    k=0.3
    
    # In milliseconds
    dt = 0.1
    
    Vs0=-60.0 #initial somatic memberane voltage
    ns0=0.004 #initial value of gating var n
    hs0=0.33 #initial value of gating var h
    Vd0=-50. #initial dendritic membrane voltage
    C0=0.03 #initial calcium 2+ balance
    l0=0.93 # initial value of IP3 channel gating variable
    
    #set here so that only the variables are in param_names
    _param_names = sorted([i for i in dir() if "__" not in i])
    
    #let user know if parameter dictionary passed in to setVars is incomplete
    warn_for_missing_vars = False
    
    def __init__(self, param_dict={}):
        self.setVars(param_dict)
    
    def setVars(self,var, warn = warn_for_missing_vars):
        
        for i in self._param_names:
            try:
                self.__dict__[i] = float(var[i])
            except Exception as e:
                self.__dict__[i] = self.__class__.__dict__[i]

                #if is missing from var dictionary, let user know or not, depending on boolean val
                if  warn:
                    print ('Param '+i+" is either missing or invalid. Using .")
                 

    def printVars(self,cols = 3):
        line = "-----------------------"
        print line + "\n"+"Model Parameters: "+"\n"+line+"\n"
        print "initial state: Vs = {}, ns = {}, hs = {}, Vd = {}, C = {}, l = {} ".format(self.Vs0, self.ns0, self.hs0, self.Vd0, self.C0, self.l0)
        print
        self.printVarsHelper(3,['Vs0', 'ns0', 'hs0', 'Vd0', 'C0', 'l0'])
        print "\n"+line + "\n"+"End Parameters "+"\n"+line+'\n'
        
            

    def model(self, y, t):

        #Vs,ns,hs,Vd,C,l = y
        Vs=y[0] #initial somatic memberane voltage
        ns=y[1] #initial value of gating var n
        hs=y[2] #initial value of gating var h
        Vd=y[3] #initial dendritic membrane voltage
        C=y[4] #initial calcium 2+ balance
        l=y[5] # initial value of IP3 channel gating variable

        # SOMATIC FUNCTIONS
        minfs = 1/(1+exp((Vs-self.vm) /self.sm))
        ninfs = 1/(1+exp((Vs-self.vn) /self.sn))
        minfps = 1/(1+exp((Vs-self.vmp)/self.smp))
        hinfs = 1/(1+exp((Vs-self.vh) /self.sh))

        tauns = self.taunb/cosh((Vs-self.vn)/(2*self.sn))
        tauhs = self.tauhb/cosh((Vs-self.vh)/(2*self.sh))

        # DENDRITIC FUNCTIONS
        #Calculate ER Ca
        Ce = (self.Ct - C)/self.sigma
        # Flux of Ca from ER to cytosol(regulated by IP3 receptors)
        J_ER_in=(self.LL + self.P*((self.IP*C*l/((self.IP+self.Ki)*(C+self.Ka)))**3))*(Ce - C) 
        # Flux from cytosol back to ER (controlled by SERCA pumps)
        J_ER_out=self.Ve*(C**2)/((self.Ke**2)+(C**2))
        # Activation of calcium current (I_can)
        caninf = 1/(1+((self.Kcan/C)**self.ncan))

        # CURRENT EXPRESSIONS
        # currents in Soma
        I_nas = self.gna*(minfs**3)*(1-ns)*(Vs-self.ena) #sodium current
        I_ks = self.gk*(ns**4)*(Vs-self.ek) #potassium current
        I_naps = self.gnaps*minfps*hs*(Vs-self.ena) #persistent sodium current
        I_L = self.gL*(Vs-self.eL) #pA
        I_sd = self.gc*(Vs-Vd)/(1-self.k) # modification of dendritic current due to somatic current
        #I_tonic = self.gtonic*(Vs-eSyn)
        # currents in Dendrite
        I_can = self.gcan*caninf*(Vd-self.ena) # calcium current
        I_ds = self.gc*(Vd-Vs)/self.k # modification of somatic current due to dendritic current

        #### DIFFERENTIAL EQUATIONS
        # SOMATIC EQUATIONS
        dVs = (-I_ks-I_nas-I_naps-I_L-I_sd+self.Iaps)/self.Cms
        dns = (ninfs-ns)/tauns
        dhs = (hinfs-hs)/tauhs
        
        # DENDRITIC EQUATIONS
        dVd = (-I_can-I_ds)/self.Cmd
        dC = (self.fi/self.Vi)*( J_ER_in - J_ER_out)
        dl = self.A*( self.Kd - (C + self.Kd)*l )

        dy = [dVs, dns, dhs, dVd, dC, dl]
        #print dy
        return dy

    def jacobian(self, y, t):
        #Vs,n,h,Vd,C,l = y
        Vs=y[0] #initial somatic memberane voltage
        n=y[1] #initial value of gating var n
        h=y[2] #initial value of gating var h
        Vd=y[3] #initial dendritic membrane voltage
        C=y[4] #initial calcium 2+ balance
        l=y[5] # initial value of IP3 channel gating variable

        return [[(-h*self.gnaps/(exp((Vs - self.vmp)/self.smp) + 1) + h*self.gnaps*(Vs - self.ena)*exp((Vs - self.vmp)/self.smp)/(self.smp*(exp((Vs - self.vmp)/self.smp) + 1)**2) - n**4*self.gk - self.gL - self.gc/(-self.k + 1) - self.gna*(-n + 1)/(exp((Vs - self.vm)/self.sm) + 1)**3 + 3*self.gna*(Vs - self.ena)*(-n + 1)*exp((Vs - self.vm)/self.sm)/(self.sm*(exp((Vs - self.vm)/self.sm) + 1)**4))/self.Cms, (-4*n**3*self.gk*(Vs - self.ek) + self.gna*(Vs - self.ena)/(exp((Vs - self.vm)/self.sm) + 1)**3)/self.Cms, -self.gnaps*(Vs - self.ena)/(self.Cms*(exp((Vs - self.vmp)/self.smp) + 1)), self.gc/(self.Cms*(-self.k + 1)), 0, 0], [(-n + 1/(exp((Vs - self.vn)/self.sn) + 1))*sinh((Vs - self.vn)/(2*self.sn))/(2*self.sn*self.taunb) - exp((Vs - self.vn)/self.sn)*cosh((Vs - self.vn)/(2*self.sn))/(self.sn*self.taunb*(exp((Vs - self.vn)/self.sn) + 1)**2), -cosh((Vs - self.vn)/(2*self.sn))/self.taunb, 0, 0, 0, 0], [(-h + 1/(exp((Vs - self.vh)/self.sh) + 1))*sinh((Vs - self.vh)/(2*self.sh))/(2*self.sh*self.tauhb) - exp((Vs - self.vh)/self.sh)*cosh((Vs - self.vh)/(2*self.sh))/(self.sh*self.tauhb*(exp((Vs - self.vh)/self.sh) + 1)**2), 0, -cosh((Vs - self.vh)/(2*self.sh))/self.tauhb, 0, 0, 0], [self.gc/(self.Cmd*self.k), 0, 0, (-self.gc/self.k - self.gcan/((self.Kcan/C)**self.ncan + 1))/self.Cmd, -self.gcan*self.ncan*(self.Kcan/C)**self.ncan*(Vd - self.ena)/(C*self.Cmd*((self.Kcan/C)**self.ncan + 1)**2), 0], [0, 0, 0, 0, self.fi*(2*C**3*self.Ve/(C**2 + self.Ke**2)**2 - 2*C*self.Ve/(C**2 + self.Ke**2) + (-1 - 1/self.sigma)*(C**3*l**3*self.IP**3*self.P/((C + self.Ka)**3*(self.IP + self.Ki)**3) + self.LL) + (-C + (-C + self.Ct)/self.sigma)*(-3*C**3*l**3*self.IP**3*self.P/((C + self.Ka)**4*(self.IP + self.Ki)**3) + 3*C**2*l**3*self.IP**3*self.P/((C + self.Ka)**3*(self.IP + self.Ki)**3)))/self.Vi, 3*C**3*l**2*self.IP**3*self.P*self.fi*(-C + (-C + self.Ct)/self.sigma)/(self.Vi*(C + self.Ka)**3*(self.IP + self.Ki)**3)], [0, 0, 0, 0, -l*self.A, self.A*(-C - self.Kd)]]
   
    def simulate(self, simulationTime, use_jacobian=True, autosave_dir = 'autosaved_sim_data'):
        #make time array and initial state array for odeint function
        t = linspace(0,simulationTime, simulationTime/self.dt)
        initial_state = [self.Vs0, self.ns0, self.hs0, self.Vd0, self.C0, self.l0]
        try:
            if use_jacobian:
                y = integrate.odeint(self.model, initial_state, t, Dfun=self.jacobian)
            else:
                raise Exception() #to be caught be except, so that integration is tried w/out jacobian
        except Exception as e:
            try:
                y = integrate.odeint(self.model, initial_state, t)
                print "Running simulation without jacobian.", e
            except Exception as a:
                Raise("Could not run simulation. Exceptions: %s"%a)
                
        V = y.T[0] #extract membrane voltage

        #save data to filename based on date and time
        if autosave_dir:
            sim_data = pandas.DataFrame(data=y, index=t, columns=['Vs0', 'ns0', 'hs0', 'Vd0', 'C0', 'l0'])
            self._autosave(autosave_dir, sim_data)

        return V, t
