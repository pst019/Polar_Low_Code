#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Mar 13 10:44:11 2021

@author: pst019


Run in shell: python Passing_argument.py argument
"""

import sys

data = int(sys.argv[1])
print("In python code")
print(data)

print(data +1)
print(type(data))

if data == 1:
    print('hello')