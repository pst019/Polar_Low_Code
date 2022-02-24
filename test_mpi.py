#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar  8 12:16:38 2021

@author: pst019
"""

from mpi4py import MPI

comm = MPI.COMM_WORLD
size = comm.Get_size()
print('size', size)
rank = comm.Get_rank()
print('rank', rank)