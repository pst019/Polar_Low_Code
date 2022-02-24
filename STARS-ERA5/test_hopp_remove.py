#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 11 16:46:14 2020

@author: pst019
"""

S= Stoll.loc[ind]
Snode_rr= [int(k) for k,_g in itertools.groupby(S.node.fillna(-1) )] #dr - repetitive remove , removes if several instances after each other are the same    


for i,node in enumerate(Snode_rr):
    if i == 0:#the first
        Snode_rr_rem_hopp= [node]
    elif i < len(Snode_rr)-1: #all middle ones
        if Snode_rr[i-1]!= Snode_rr[i+1]:
            if Snode_rr_rem_hopp[-1] != node:
                Snode_rr_rem_hopp += [node]
        elif Snode_rr[i-1] == Snode_rr[i+1] and Snode_rr[i-1]== -1:
            #nodes between nans should not be removed
            if Snode_rr_rem_hopp[-1] != node:
                Snode_rr_rem_hopp += [node]
    elif i == len(Snode_rr)-1: #the last
        if Snode_rr_rem_hopp[-1] != node:
            Snode_rr_rem_hopp += [node]        
        
print(Snode_rr)
print(Snode_rr_rem_hopp)