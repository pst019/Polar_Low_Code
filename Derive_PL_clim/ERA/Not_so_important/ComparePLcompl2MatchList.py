#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct  3 11:36:27 2017

@author: pst019
"""
import os
user = os.getcwd().split('/')[2]

import numpy as np
import sys  #to import the functions from a different directory

sys.path.insert(0, '/home/'+user+'/polar_low_code/Functions/')
from f_imp_ERA2 import *
from f_impLists import *    


model= 'ERA'
model= 'ASR'

Tname='6' #version of the complet PL list

yearT= 2003 # year in which the number of TRACK PLs is counted

""" import TRACKmatchSTARS matchlist
 all the matches, tPLstart gives the starting point in time of the match"""

if model== 'ERA':
    filename= Mediadir+'PL/PLclim/MatchList/TRACKmatchSTARS.csv' #_Sdurlim3.csv'
elif model== 'ASR':
    filename= Mediadir+'PL/PLclim/MatchList/ASR_TRACKmatchSTARS.csv' #_Sdurlim3.csv'
    
    
import pandas as pd
PLlist= np.array(pd.read_table(filename, sep="\t"))

"""combine double PLs in TRACKmatchSTARS"""
doubleTRACKS= remove_dublicate([x for x in PLlist[:, 1] if list(PLlist[:, 1]).count(x) > 1])
print('double TRACKs', doubleTRACKS)

removeind= []
for n in range(len(doubleTRACKS)): #remove the double tracks
    dind=  np.where(PLlist[:, 1] == doubleTRACKS[n])[0]
    tPLstartfirst= PLlist[:, 4][dind[0]]
    tPLstartsecond= PLlist[:, 4][dind[-1]]
    tPLendfirst= PLlist[:, 5][dind[0]]
    if tPLstartsecond >= tPLstartfirst  and tPLstartsecond <=  tPLendfirst:
        print('merge, SPLs:', PLlist[:,0][dind[0]], PLlist[:,0][dind[-1]], 'at times:', 'startPL1', tPLstartfirst, 'StartPL2', tPLstartsecond, 'endPL1', tPLendfirst)
        
        PLlist[dind[0], 5] = PLlist[dind[-1], 5]  #put the end of the second polar low for the end of the first pl
        removeind += [dind[-1]]

PLlist= np.delete(PLlist, removeind, 0) #remove the second PL

print('nr of PLs:', len(PLlist[:, 0]))
                 
SPLnr= PLlist[:, 0]
TPLnr= PLlist[:, 1]
yearL= PLlist[:, 2]
monthL= PLlist[:, 3]

if model == 'ERA':
    tPLstart= (PLlist[:, 4]+ 1)//2
    tPLend= PLlist[:, 5]//2
elif model == 'ASR':
    tPLstart= PLlist[:, 4]
    tPLend= PLlist[:, 5]

nrSTARSPLs= 0       
inclSPLlist= []       
              
for i in range(len(SPLnr)):
    """ import the complete PL list"""
    TPL=TRACK_PLlist_2(yearL[i], monthL[i], hemi= '', model= model, name=Tname)
    
#    nrtTPLi= len(np.where(TPL.PLli1st[0] == TPLnr[i])[0]) #amount of time steps which the TPL is included in the TPL list
    nrtTPLi= len(np.where(np.logical_and.reduce((TPL.PLlist[1] >= tPLstart[i], TPL.PLlist[1] <= tPLend[i], TPL.PLlist[0] == TPLnr[i] , TPL.PLlist[-1] == 1)))[0]) #amount of time steps which the TPL is included in the TPL list

    if nrtTPLi > 0:
        nrSTARSPLs += 1
        inclSPLlist+= [SPLnr[i]]
        
#    print( SPLnr[i], TPLnr[i], nrtTPLi, nrSTARSPLs)
#    print('STARSPLnr', SPLnr[i], 'matching TRACK', TPLnr[i], 'length of TRACK in derived PL list', nrtTPLi, nrSTARSPLs)


nTPL= 0 #number of TRACK PLs
for month_i in [1,2,3,4,10,11,12]:
    """ import the complete PL list"""
    TPL=TRACK_PLlist_2(yearT, month_i, hemi= '', model= model, name=Tname)
    nTPL += len(remove_dublicate(TPL.PLlist[0]))
    
    print(month_i, nTPL)
    
    
    
print('number of STARS PLs in complete PL list', nrSTARSPLs)
print('number of TRACK PLs in complete PL list in that year', nTPL)