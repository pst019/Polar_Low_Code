#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec  2 13:48:03 2019

@author: pst019
"""
import numpy as np
import os
from shutil import copyfile


file= "txt_files/track_list_PMC.txt" #used until version 3 - non-python length: 173
file= "txt_files/track_list_PMC_appent.txt" #used from version 4 - non-python length: 170

trackid, startyear, startmonth, startday, endyear, endmonth, endday= np.loadtxt(file).T.astype(int)
#trackid is actually tracktimeid, so a count for the time block in which PLs occurred

PMC_dir= "/media/pst019/1692A00D929FEF8B/ERA5_STARS/pmctrack-master/"
os.chdir(PMC_dir)


"""run through the indexes"""
i= 7 #this is python counting
outname=''
#outname='ec_incl_'
#outname='smth'
#outname='h_'
#test='Denis_1'
#test='Denis_1h'

#test='test_nolsm'
#test='test_smth40'
#test='version2'
#test='version2_gamma01'
#test='version2_gamma01smth30'
#test='test_1'
#test='version3'
test='version4'


#test='version1_nosmth'
#test='version1_no2merge'


for i in range(1, 170): #173 #non-python counting starting with 1
#for i in range(i, i+1):
    i -= 1
    conf_file= "conf_files/settings_ERA5_"+test+"_"+str(trackid[i]).zfill(3)+".conf"

    if 'Denis' in test:   copyfile("settings_ERA5_Denis.conf", conf_file )
    else:  copyfile("settings_ERA5.conf", conf_file )
    
    """rewrite the config file"""
    with open(conf_file, 'r') as file :
      filedata = file.read()
    
    # Replace the target string
    filedata = filedata.replace('1999_12_18',
                                str(startyear[i])+ '_'+str(startmonth[i]).zfill(2)+ '_' 
                                +str(startday[i]).zfill(2) )
    
    filedata = filedata.replace('1999_12_19',
                                str(endyear[i])+ '_'+str(endmonth[i]).zfill(2)+ '_' 
                                +str(endday[i]).zfill(2) )
    
    filedata = filedata.replace('outdir="../tracks/conf"', 'outdir="../tracks/'+test+'/tracks_'+str(trackid[i]).zfill(3)+'"')
    
    # Write the file out again
    with open(conf_file, 'w') as file:
      file.write(filedata)
      
    
    """run the PMC tracking"""
    os.chdir(PMC_dir)
    os.system('./track.x '+conf_file)


#orig: vo:850, zeta_max (min): 3 (2.5), del_psea_min=.5
#t1 :  vo:925, zeta_max (min)= 2 (1.5), del_psea_min=50