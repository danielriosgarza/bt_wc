# -*- coding: utf-8 -*-
"""
Created on Fri Jul 16 08:57:06 2021

@author: u0139894
"""

import numpy as np
import scipy.integrate as solver
import pandas as pd




class Strain:
    '''
    Parameters of a single strain.
    
    Parameters
    ---------
    
    name : str
    strain id
    
    x_0 : float
    initial abundance
    
    q_0 : float
    initial lag phase term
    
    mu : float
    max growth rate
    
    dr: float
    death rate. Assuming a fixed death rate
    
    metab : list
    ordered list of metab ids that are imported or secreted by the strain.
    
    metab_v : np.array/list
    consumption/production of metabolites in the metab list.
    (positive - consumption; negative - production)
             
    metab_k : float
    Monod constants for the metabolites in the metab list
    
    growth_model : list of tuples
    composition of the growth rate term. In each tuple in the list, the first term
    must be a constant followed by a list of metabolites. The constant and the monod
    equations for each of the metabolites in a tuple are multipled. 
    the final growth rate consists of the sum of these products. An arbitrary number 
    of terms can be used.
    E.g. [(0.1, 'Fructose', 'Pyruvate'), (0.8, 'acetate')]
    
    gr = 0.1*((s_fruc/(k_fruc + s_fruc) * (s_pyr/(k_pyr + s_pyr)) + 
               0.8*(s_ac/(k_ac + s_ac)
                  
    feeding : list of tuples
    terms used in the importing/secreting of metabs. Binary indices that
    indicate which growth rate terms to use for each metabolite in the 
    metab list.
    
    E.g. In the example of the growth model, the strain has two grwoth 
    terms ([(0.1, 'Fructose', 'Pyruvate'), (0.8, 'acetate')])
    For the feeding/secreting we need to determine for each of the
    three metabolites which growth rate term to take into account.
    So, if the first term is used for fructose and pyruvate and both
    terms are used for acetate, we have:
    
        [(1,0), (1,0), (1,1)]
    
    
    parametrizations: list
    A list with strain objects and alternative
    parametrizations (e.g. monoculture parameters) (optional)
    
    parametrization_in_use : int
    index of the parametrization to use
    
   
    '''
    def __init__(self, name, x_0, q_0, mu, dr, metab, metab_v, metab_k, growth_model, feeding):
        
        self.name = name
        self.x_0 = x_0
        self.q_0 = q_0
        self.mu = mu
        self.dr = dr
        
        self.metab = np.array(metab)
        self.metab_v = np.array(metab_v)
        self.metab_k = np.array(metab_k)
        self.growth_model = growth_model
        self.feeding = np.array([np.array(feeding[i])*self.metab_v[i] for i in range(len(self.metab))]) 
        
        
        self.parametrizations=None
        self.parametrization_in_use = None




class Culture:
    
    '''
    Parameters
    ----------
    
    strains : list
    list with defined 'Strain' objs
    
    metabolome : list
    list with the names of external metabolites
    
    metabolome_c : list
    list with the concentration of external metabolites
    
    feed_c : np.array
    same dimension as 'metabolome_c'
    the concentration of metabolites in the feed
    
    
    dilution : float
    dilution parameter for chemostat dynamics
    
    
    '''
    def __init__(self, strains, metabolome, metabolome_c, feed_c = None, metab_deg = None, dilution = 0):
        
        
        self.strain_obj = self.__take_strains(strains)
        self.nstrains = len(self.strain_obj)
        self.strains = np.array([i.name for i in strains])
        self.metabolome = np.array(metabolome)
        self.metabolome_c = np.array(metabolome_c)
        self.nmets = len(self.metabolome_c)
        self.ndims = self.nstrains + self.nmets
        self.metab_v, self.metab_k, self.mus, self.qs, self.x_0, self.drs = self.__parse_strains(self.strain_obj)
        self.feeding_m = self.__get_feeding_matrices(self.strain_obj)
        self.growth_rate_functions = [self.__make_growth_function(i) for i in self.strain_obj]
        self.community_dyn = None
        self.environment_dyn = None
        self.lag_dyn = None
        self.system_time = None
        
        self.feeding_c = np.array(feed_c)
        self.metab_deg = np.array(metab_deg)
        
        self.dilution = 0
        
        if dilution is not None:
            self.dilution = dilution
    
    def __take_strains(self, strains):
        '''
        Check for the strain parametrization to use.

        Parameters
        ----------
        strains : list
        List of Strain objs
        
        Returns
        -------
        strain_objs : list
        returns a list of Strain objs with their parametrization
        indicated by their 'parametrization_in_use' in use parameter.

        '''
        strain_objs = []
        
        for i in strains:
            if i.parametrization_in_use is not None:
                strain_objs.append(i.parametrizations[i.i.parametrization_in_use])
            else:
                strain_objs.append(i)
                
        return strain_objs
    
        
    
    def __parse_strains(self, strains_obj):
        '''
        Retrieve the parameters for the strains in the culture.

        Parameters
        ----------
        strains_obj : list
        List of Strain objs

        Returns
        -------
        metab_v : dict
        [strain] : {v for each metabolite in self.metabolome}

        metab_k : dict
        [strain] : {k for each metabolite in self.metabolome}


        mus : np.array
        max growth rates in the same order as the strains list
        
        
        qs : np.array
        qs in the same order as the strains list

        x_0 : np.array
        x_0 in the same order as the strains list

        '''
        
        
        
        metab_v= {i: np.zeros(self.nstrains) for i in self.metabolome}
        metab_k= {i: np.zeros(self.nstrains) for i in self.metabolome}
        
        mus = np.zeros(self.nstrains)
        qs = np.zeros(self.nstrains)
        x_0 = np.zeros(self.nstrains)
        drs = np.zeros(self.nstrains)
        
        
        for i,v in enumerate(strains_obj):
            
            mus[i] = v.mu
            qs[i] = v.q_0
            x_0[i] = v.x_0
            drs[i] = v.dr
            
            for z,v2 in enumerate(v.metab):
                if v2 in self.metabolome:
                    metab_v[v2][i],metab_k[v2][i]  = v.metab_v[z],v.metab_k[z]
        
        metab_v = np.array([metab_v[i] for i in self.metabolome])
        metab_k = np.array([metab_k[i] for i in self.metabolome])
        
        return metab_v, metab_k, mus, qs, x_0, drs
    
    
    def __get_feeding_matrices (self, strains_obj):
        '''
        Build feeding matrices per strain.

        Parameters
        ----------
        strains_obj : list
        List of Strain objs

        Returns
        -------
        matrices : list of np.arrays
        order list of feeding matrices. Each matrix indicates
        the coeficient v applied to the growth rate term 
        for each metabolite

        '''
        
        matrices=[]
        for i,v in enumerate(strains_obj):
            fm = np.zeros((self.nmets, len(v.feeding[0])))
            
            for z,m in enumerate(self.metabolome):
                if m in v.metab:
                    fm[z] = v.feeding[list(v.metab).index(m)]
            matrices.append(fm)
        
        return matrices
    
    def __make_growth_function(self, strain_obj):
        '''
        Make a growth rate function for each strain
        

        Parameters
        ----------
        strains_obj : list
        List of Strain objs

        Returns
        -------
        list of grwoth rate functions

        '''
        
        multipliers = []
        
        for i in strain_obj.growth_model:
            d=[i[0]]
            for z in i[1::]:
                if z in self.metabolome:
                    d.append((strain_obj.metab_k[list(strain_obj.metab).index(z)], list(self.metabolome).index(z)))
            multipliers.append(d)
    
        
        
        def growth_function(s):
            '''
            

            Parameters
            ----------
            s : np.array
            concentration of external metabolits as in self.metabolome_c

            Returns
            -------
            np.array
            value of the growth rate terms.
            
    
            '''
            
            l=[]
            
            for i in multipliers:
                a = i[0]
                
                for z in i[1::]:
                    a*= s[z[1]]/(z[0] + s[z[1]])
                l +=[a]
            
            return np.array(l)
        return growth_function
    
            
                
        
        
        
        
    def lagPhase(self, q):
        """ 
        function for the lag phase.
        
        lag_term = lag_variable/(1+lag_variable)
        
        The time derivative of the lag variable is given by:
            d(lag_variable)/dt = u_max*lag_variable
        
        """
        return q/(1+q)
    
    def growthRates(self, s,q):
        """ 
        computes the growth rate term for the strain in the culture.
        Uses growth rate functions and the lagphase function
        
        gr = lag * mu_max * growth_rate terms
        
        """
        lag = self.lagPhase(q)
        gr = []
        for i,v in enumerate(self.growth_rate_functions):
            gr.append(self.mus[i]*lag[i]*v(s))
            
            
        return gr
    
    def dynamics(self, t, vars_t):
        '''
        

        Parameters
        ----------
        t : 
        integration time
        
        vars_t : np.array
        array with the current state of the system
        
        var_t[0:nstrains] = liveCell cocentrations (same order as self.strains)
        
        var_t[nstrains: nstrains + nmets] = metabolite concentrations 
        (same order as self.metabolome)
        
        var_t[nstrains + nmets::] = lagphase_variable (same order as self.strains)
        
        Returns
        -------
        None.

        '''
        
        x = vars_t[0:self.nstrains]
        s = vars_t[self.nstrains : self.nstrains + self.nmets]
        q = vars_t[self.nstrains + self.nmets:]
        
        x = x*(x>0)
        s = s*(s>0)
        
        growth = self.growthRates(s,q)
        
        #derivative for bacterial growth
        dxdt = np.array([sum(i) for i in growth])*x - self.dilution*x - self.drs*x
        
        #derivative for the change in substrate concentration
        dsdt =  np.zeros(self.nmets)
        
        #consumption/secretion by microbes
        for i in range(self.nstrains):
            feedk=self.feeding_m[i].dot(growth[i]*x[i])
            dsdt-=feedk
        
        #dilution
        
        
        
        dsdt += self.dilution * (self.feeding_c-s)
        
        if self.metab_deg is not None:
            dsdt-=self.metab_deg * s
        #derivative of the lagphase
        dqdt = self.mus * q
        
        ddt = np.append(dxdt, np.append(dsdt,dqdt))
        
        return ddt
        
    
    def simulate(self, t_start=0, t_end = 50, nsteps = 1000, method = 'bdf'):
        """ 
        Simulate the ODE's 
        
        """
        # Solve ODE
        
        ode = solver.ode(self.dynamics)
        
        ode.set_integrator('vode',nsteps=nsteps,method= method)
        
        y_init = np.append(self.x_0, np.append(self.metabolome_c, self.qs))
        
        # Time
       
        t_step = (t_end - t_start)/nsteps
        
        ode.set_initial_value(y_init,t_start)
        
        ts = []
        ys = []
        
        while ode.successful() and ode.t < t_end:
            ode.integrate(ode.t + t_step)
            ts.append(ode.t)
            ys.append(ode.y)
        
        
        #integration finished, store the results
        time = np.array(ts)
        self.system_time=time
        y_total = np.vstack(ys).T    #y = (x,s,q)
        
        community_dyn = {}
        environment_dyn = {}
        lag_dyn={}
        
        for i,v in enumerate(self.strains):
            community_dyn[v] = y_total[i]
        
        for i in range(self.nmets):
            environment_dyn[self.metabolome[i]]=y_total[i+self.nstrains]
        for i,v in enumerate(self.strains):
            lag_dyn[v] = y_total[i+self.nstrains + self.nmets]
        
        self.community_dyn = pd.DataFrame.from_dict(community_dyn)
        self.community_dyn.index=time
        self.environment_dyn = pd.DataFrame.from_dict(environment_dyn)
        self.environment_dyn.index=time
        self.lag_dyn = pd.DataFrame.from_dict(lag_dyn)
        self.lag_dyn.index = time
        #update sharing of resources
        





