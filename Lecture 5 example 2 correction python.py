# -*- coding: utf-8 -*-
"""
Created on Sun Nov  4 12:57:50 2018

@author: lemitri
"""

import os
import pandas as pd
import scipy.stats as sp
import gurobipy as gb
import numpy as np

class expando(object):
    '''    
    A small class which can have attributes set
    '''
    pass


class example_2:
    def __init__(self):
        self.data = expando()
        self.variables = expando()
        self.constraints = expando()
        self._build_model()

    
    def optimize(self):
        self.model.optimize()
    
    def _build_model(self):
        
        self.model = gb.Model()
        self._build_variables()
        self._build_objective()
        self._build_constraints()
    
    def _build_variables(self):
        
        #indexes shortcuts 
        m = self.model
            

        self.variables.x =  m.addVar(lb=2,name='x')

        self.variables.y =  m.addVar(lb=-gb.GRB.INFINITY,name='y')

        self.variables.mu_1 =  m.addVar(name='mu 1') # positive dual variable
        self.variables.mu_2 =  m.addVar(name='mu 2') # positive dual variable
        self.variables.mu_3 =  m.addVar(name='mu 3') # positive dual variable
        self.variables.mu_4 =  m.addVar(name='mu 4') # positive dual variable

        self.variables.l_1 =  m.addVar(name='l 1') # positive auxiliary variable
        self.variables.l_2 =  m.addVar(name='l 2') # positive auxiliary variable
        self.variables.l_3 =  m.addVar(name='l 3') # positive auxiliary variable
        self.variables.l_4 =  m.addVar(name='l 4') # positive auxiliary variable
                   
        m.update()

    def _build_objective(self): # building the objective function of the leader!

        #indexes shortcuts    
        m = self.model      
        
        m.setObjective(-self.variables.y,   
            gb.GRB.MINIMIZE)         

    
    def _build_constraints(self):

        #indexes shortcuts     
        m = self.model

        # lower level inequality constraints 
        
        self.constraints.c_1=m.addConstr(
                    - 7 + self.variables.y + self.variables.x,
                    gb.GRB.EQUAL,
                    self.variables.l_1)

        self.constraints.c_2=m.addConstr(
                    9 - self.variables.y,
                    gb.GRB.EQUAL,
                    self.variables.l_2)
        
        
        self.constraints.c_3=m.addConstr(
                    1 + 2*self.variables.y - self.variables.x,
                    gb.GRB.EQUAL,
                    self.variables.l_3)
        
        
        self.constraints.c_4=m.addConstr(
                    45 - self.variables.y - 6*self.variables.x,
                    gb.GRB.EQUAL,
                    self.variables.l_4)

        # Stationarity condition
        
        self.constraints.L=m.addConstr(1-self.variables.mu_1+self.variables.mu_2-2*self.variables.mu_3+self.variables.mu_4,
                                       gb.GRB.EQUAL,
                                       0)
        
        # complementarity constraints

        self.constraints.SOS_1=m.addSOS(gb.GRB.SOS_TYPE1,[self.variables.l_1,self.variables.mu_1])
        self.constraints.SOS_2=m.addSOS(gb.GRB.SOS_TYPE1,[self.variables.l_2,self.variables.mu_2])
        self.constraints.SOS_3=m.addSOS(gb.GRB.SOS_TYPE1,[self.variables.l_3,self.variables.mu_3])        
        self.constraints.SOS_4=m.addSOS(gb.GRB.SOS_TYPE1,[self.variables.l_4,self.variables.mu_4])


solution = example_2()
solution.optimize()

print('x=',solution.variables.x.x)
print('y=',solution.variables.y.x)