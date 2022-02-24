#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 15 16:15:40 2021

@author: pst019
"""

import time
import multiprocessing

start = time.perf_counter()

def do_something():
    print('Sleep 1 sec...')
    time.sleep(1)
    print('Done Sleeping..')
    
"""synochroniesly"""
#do_something()
#do_something()


"""manually"""
# p1= multiprocessing.Process(target=do_something)
# p2= multiprocessing.Process(target=do_something)


# p1.start()
# p2.start()


# p1.join()
# p2.join()


"""loop"""
# processes = []

# for _ in range(10):
#     p = multiprocessing.Process(target=do_something)
#     p.start()
#     processes.append(p)
    
 
# for process in processes:
#     process.join()       


"""put arguments"""

# def do_something(seconds):
#     print(f'Sleep {seconds} sec...')
#     time.sleep(seconds)
#     print('Done Sleeping..')

# processes = []

# for _ in range(10):
#     p = multiprocessing.Process(target=do_something, args=[1.5])
#     p.start()
#     processes.append(p)
    
 
# for process in processes:
#     process.join()   


"""process pool executer - do not need join"""

def do_something(seconds):
    print(f'Sleep {seconds} sec...')
    time.sleep(seconds)
    return f'Done Sleeping..{seconds}'

import concurrent.futures #instead multiprocessing

"""manual"""
# with concurrent.futures.ProcessPoolExecutor() as executor:
#     f1= executor.submit(do_something, 1)
#     f2= executor.submit(do_something, 1)

#     print(f1.result())
#     print(f2.result())

"""list comprehension"""
# with concurrent.futures.ProcessPoolExecutor() as executor:
#     secs= [5, 4, 3, 2, 1]
#     results= [executor.submit(do_something, sec ) for sec in secs]

#     for f in concurrent.futures.as_completed(results):
#         print(f.result())

"""print in starting order"""
with concurrent.futures.ProcessPoolExecutor() as executor:
    secs= [5, 4, 3, 2, 1]
    results= executor.map(do_something, secs)
    
    # for result in results:
    #     print(result)


finish= time.perf_counter()

print(f'Finished in {round(finish - start,2)} seconds')




"""zip to loop over two lists at the same time
for name, hero in zip(names, heroes):
    print(f'{name} us actually {hero}')

""" 