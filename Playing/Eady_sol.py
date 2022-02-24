# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
import numpy as np

mu= 1.6
k= mu

U= 1

c_i = 0.2*U
c_r= 0.5*U


z= np.linspace(-1,1, 20)
x= np.linspace(0, 2* np.pi, 30)


phi= np.cosh(mu * z)  - U* c_r * np.sinh(mu * z) /(mu * (c_i**2 + c_r**2))


