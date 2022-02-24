#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar  8 12:54:09 2021

@author: pst019
"""


import numpy as np
from mpi4py import MPI
import time
comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()
print(rank)
  
nx=10000024
Lx = 2*np.pi
x_b = np.linspace(0,Lx,nx)
nx_r = int(nx/size)
x_s = np.empty(nx_r)
comm.Scatter(x_b,x_s,root = 0)

a= np.sin(x_s)
b = np.cos(x_s)
c_s = rank*(a**2 + b**2)
time.sleep(5)
c_b = np.empty(x_b.shape)
comm.Gather(c_s,c_b,root=0)
if rank ==0:
    print(c_b)