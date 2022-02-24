#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 21 12:40:31 2017

@author: pst019

Write data into a csv file
"""

import csv


file='/home/'+user+'/Documents/SIC-PCs.csv'

column_name= ['Year', 'PC_1', 'PC_2', 'PC_3', 'PC_4']

data= np.hstack((year.reshape(-1, 1), pcs))  #shape (36, 5)
#for hstack shape(year)= (36, 1) and shape(pcs)= (36,4)


with open(file, 'w') as csvfile:
    writer = csv.writer(csvfile, delimiter="\t")
    writer.writerow(column_name)
    [writer.writerow(r) for r in data]  