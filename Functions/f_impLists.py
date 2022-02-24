#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 13 15:07:08 2017

@author: pst019

This is from plot_comparison to import data
"""

from netCDF4 import Dataset
import numpy as np
import time
import pandas as pd
import csv
import datetime

import os
user = os.getcwd().split('/')[2]

homedir= '/home/'+user+'/home/'
Mediadir= '/media/'+user+'/PatsOrange/'


import sys
sys.path.insert(0, homedir +'Polar_Low/Code2/Derive_PL_clim/ERA/')
from f_imp_ERA2 import *


class GunnarsList:
    def __init__(self, filename=homedir +'Polar_Low/Code2/PolarLowList_Gunnar_adapted.csv',
                 windlim= -1):
        """import the exel list made from Gunnars STARS dataset
        windlim specifies at which wind limit the polar lows should be taken
        used before STARS database with hourly resolution was available"""
        with open(filename, 'r') as f:
            reader = csv.reader(f)
            PLlist = np.array(list(reader))[1:]
            
        len_G= len(PLlist[:,0])
        year_G, month_G, day_G, hour_G, lat_G, lon_G, MSLP_G, wind_G= np.zeros((len_G)), np.zeros((len_G)), np.zeros((len_G)), np.zeros((len_G)), np.zeros((len_G)), np.zeros((len_G)), np.zeros((len_G)), np.zeros((len_G))
        for n in range(len_G):
        #    print(PLlist[n,2][:2])
            if PLlist[n,0]!= '':   #only take coloumns which are not empty
                year_G[n]=PLlist[n,0][-4:]
                month_G[n]=PLlist[n,0][-7:-5]
                day_G[n]=PLlist[n,0][:2]
                
                if len(PLlist[n,1]) >=4: hour_G[n]=PLlist[n,1][:2]
                elif len(PLlist[n,1]) >=3: hour_G[n]=PLlist[n,1][:1]
                else: hour_G[n]=PLlist[n,1]
                    
                lat_G[n]= PLlist[n,2][:2]
                
                if PLlist[n,3][-1]== 'E':
                    lon_G[n]= PLlist[n,3][:2]
                else:
                    lon_G[n]= -1*int(PLlist[n,3][:2])
                    
                if PLlist[n,6] not in ['','NA']: wind_G[n]= float(PLlist[n,6][:2])*0.514
                if PLlist[n,5] not in ['','NA', '**']: MSLP_G[n]= float(PLlist[n,5])
            
        
        wind_G= wind_G[year_G!=0]
        self.month= month_G[year_G!=0][wind_G >= windlim]
        self.day= day_G[year_G!=0][wind_G >= windlim]
        self.hour= hour_G[year_G!=0][wind_G >= windlim]
        self.lat= lat_G[year_G!=0][wind_G >= windlim]
        self.lon= lon_G[year_G!=0][wind_G >= windlim]
        self.MSLP= MSLP_G[year_G!=0][wind_G >= windlim]
        self.year= year_G[year_G!=0][wind_G >= windlim]
        self.wind=wind_G[wind_G >= windlim]
        
        self.tPL= (self.day-1)*4 + np.round(self.hour/6)


class STARSList:
    def __init__(self, syear=None, smonth=None, filename='STARS_TRACKS.csv', durationlim= 6, timestep= 'sixhours'):
        """import the STARS database"""
        filedir= Mediadir +'/PL/PLclim/STARS/'
        
        with open(filedir+ filename, 'r') as f:
            reader = csv.reader(f)
            PLlist = np.array(list(reader)[2:])
        
        len_G= len(PLlist)
        PLnumber, year, month, day, hour, lat, lon= np.zeros((len_G)), [0]*len_G, [0]*len_G, [0]*len_G, [0]*len_G, np.zeros((len_G)), np.zeros((len_G)) 
        
        if filename=='STARS_TRACKS.csv':
            """for the normal list"""
            for n in range(len_G):
#                print(PLlist[n][0])
                PLnumber[n]= int(PLlist[n][0][:3])
                year[n]=int(PLlist[n][0][6:10])
                month[n]=int(PLlist[n][0][13:15])
                day[n]=int(PLlist[n][0][18:20])
                hour[n]=int(PLlist[n][0][22:24])
                lat[n]=float(PLlist[n][0][31:37])
                lon[n]=float(PLlist[n][0][39:46]) 
#                if PLnumber[n] == 136:
#                    print(PLnumber)
            
        elif filename=='STARS_TRACKS_south.csv':
            for n in range(len_G):
                """for the south list"""
                PLnumber[n]= int(PLlist[n][0][:2])+200
                year[n]=int(PLlist[n][0][4:10])
                month[n]=int(PLlist[n][0][11:15])
                day[n]=int(PLlist[n][0][16:20])
                hour[n]=int(PLlist[n][0][21:24])
                lat[n]=float(PLlist[n][0][31:37])
                lon[n]=float(PLlist[n][0][39:46])
                
        PLlist= np.vstack((PLnumber, year, month, day, hour, lat, lon))

        """calculate duration"""
        if durationlim >1:
            duration= np.zeros((len_G))
            #for n, PLnr in enumerate(PLnumber):
            for PLnr in remove_dublicate(PLnumber):
    #            print(PLnr)
                ns= np.where(PLnumber == PLnr)[0]
                starttime= datetime.datetime(year[ns[0]], month[ns[0]], day[ns[0]], hour[ns[0]])
                endtime= datetime.datetime(year[ns[-1]], month[ns[-1]], day[ns[-1]], hour[ns[-1]])
                duration[ns]= (endtime -starttime).seconds//3600 +(endtime -starttime).days* 24
#                if PLnr == 136: print('found1', starttime, endtime, (endtime -starttime).seconds//3600)
                        
            PLlist= PLlist[:, duration >= durationlim]

        """calculate tPL and keep only the one closest to 6hourly or 3hourly timestep"""
        if timestep== 'sixhours':
            tPL= (PLlist[3]-1)*4 + PLlist[4]/6
        elif timestep== 'threehours':
            tPL=(PLlist[3]-1)*8 + PLlist[4]/3
                      
        
        sixhourPLs=[]
             
        for PLnr in remove_dublicate(PLlist[0]):
            ns= np.where(PLlist[0] == PLnr)[0]  
            roundtPL=np.round(tPL[ns])
#            if PLnr== 136: print('found')
            for rtPL in remove_dublicate(roundtPL):
                diff= np.abs((tPL[ns] -roundtPL)[roundtPL == rtPL])
                if np.min(diff) <= .2:
                    index= np.where(diff== np.min(diff))[0]
                    sixhourPLs += list(np.where(PLlist[0] == PLnr)[0][roundtPL == rtPL][index])
        
        PLlist= np.vstack((PLlist, np.round(tPL)))
        PLlist= PLlist[:, sixhourPLs]

        """just for a certain months and year"""
        if syear !=None: PLlist= PLlist[:, PLlist[1]== syear]
        if smonth != None: PLlist= PLlist[:, PLlist[2]== smonth]
        
#        print(list(PLlist[0]))
        
#        """reorganize the PLnumber"""                                      
#        for n, PLn in enumerate(remove_dublicate(PLlist[0])):#rearange the list
#            PLlist[0][PLlist[0]== PLn]= n+1
        
        [self.PLnumber, self.lat, self.lon, self.tPL]= np.vstack((PLlist[0], PLlist[5:8]))
#        [self.PLnumber, self.year, self.month, self.day, self.hour, self.lat, self.lon, self.tPL]= PLlist



class JuliarsList:
    def __init__(self, filename=homedir +'Polar_Low/Code2/TablePL_2000-2004-1_JuliaSmirnova_adapted.csv'):
        """import the PL list from Julia"""
        with open(filename, 'r') as f:
            reader = csv.reader(f)
            PLlist = np.array(list(reader)[1:])
            
        PLlist[PLlist == ''] = -1 #fill value
        PLlist= np.array(PLlist, dtype= 'float')
        
        """make the polar low number"""
        listlen= np.shape(PLlist)[0]
        PLlist= np.hstack((np.zeros((listlen,1)), PLlist)) #add a column of zeros , which will become the polar low number
        
        PLnrcounter= 1
        for n in range(listlen):
            if PLlist[n, 1] == -1:
                PLnrcounter += 1
            else:        
                PLlist[n, 0]= PLnrcounter
        
        PLlist= PLlist[PLlist[:,0] != 0]
        
        self.PLnumber= PLlist[:, 0]
        self.year= PLlist[:, 1]
        self.month= PLlist[:, 2]
        self.day= PLlist[:, 3]
        self.lat= PLlist[:, 4]
        self.lon= PLlist[:, 5]
        self.hour= PLlist[:, 6]
        
        self.tPL= (self.day -1)*4 + np.round(self.hour/6)
        

        



class TRACK_PLlist_2:
    def __init__(self, year, month, hemi='', name='', model='ERA'):
        """include the list created by TRACK2PL_4 and by WritePLlistComp
        this is often used
        month= 'all' gives out the full list for all months including also information about the year"""
        
        filename= Mediadir +'/PL/PLclim/PL_file_'+model+'/PL_file_complete'+hemi+name+str(year)+'.csv'
        
        PLlist= np.array(pd.read_table(filename, sep="\t"))
        
        if month=='all':
            self.PLlist= PLlist.T
        else:
            yearlist= PLlist[:,0]
            monthlist=PLlist[yearlist== year,1]
            PLlist= PLlist[yearlist==year, 2:]
           
            self.PLlist=PLlist[monthlist==month, :].T



class TRACK_PLlist_Year:
    def __init__(self, year, name='', model='ERA'):
        """include the list created by WritePLlistComp2S"""
        filename= Mediadir +'/PL/PLclim/PL_file_'+model+'/PL_file_complete'+name+str(year)+'.csv'
        
        PLlist= np.array(pd.read_table(filename, sep="\t"))       
        self.PLlist=PLlist.T



                          
                          
class TRACK_list:
    def __init__(self, year, season):
        """import the TRACK list
        season should be given in 'A', 'B', 'C', or 'D'
        PLlist = (PLnr, duration, timestep in season from start period, longitude, latitude, vorticity) """

        print('should use TRACK_list month')
        filedir= "/media/"+user+"/1692A00D929FEF8B/PL/PLclim/ERA_NH_TRACKoutput/ERA_VOR850_"

        PLlist= []
        with open(filedir+str(year)+'.'+season+'_NH_tracks', 'r') as f:
            reader = csv.reader(f)
            for rn, row in enumerate(reader):
                if rn >2: #to exclude the first three lines
                    r= row[0]
                    if 'TRACK_ID' in r:
                        plnr= [int(r[10:])]
        
                    elif 'POINT_NUM' in r:
                        duration= [int(r[11:])]
        
                    else:
                        out= [float(x) for x in r[:-1].split(' ')]
                        PLlist += plnr + duration + out
         
            self.PLlist= np.reshape(np.array(PLlist), (-1, 6)).transpose()
                        



class TRACK_list_month:
    def __init__(self, year, month, duration_limit= 2, approved= True, timestep='sixhours'):
        """import the TRACK list for one given months - this is better than the two versions before
        the timestep specifies in which interval the time tPL is taken
        mainly version 2 is used"""
        import datetime

        if approved == True:
            filedir= "/media/"+user+"/1692A00D929FEF8B/PL/PLclim/ERA_NH_TRACKoutput/ERA_VOR850_"
        else:
            filedir= "/media/"+user+"/1692A00D929FEF8B/PL/TRACKoutput_notapproved/ERA_VOR850_"
            print('not approved files')
            
            
        if month in [1,2,3]: season, startdate= 'A', datetime.date(year, 1, 1)
        elif month in [4,5,6]: season, startdate= 'B', datetime.date(year, 4, 1)
        elif month in [7,8,9]: season, startdate= 'C', datetime.date(year, 7, 1)
        elif month in [10,11,12]: season, startdate= 'D', datetime.date(year, 10, 1)

        self.filename= filedir+str(year)+'.'+season+'_NH_tracks'
        PLlist= []
        with open(self.filename, 'r') as f:
            reader = csv.reader(f)
            for rn, row in enumerate(reader):
                if rn >2: #to exclude the first three lines
                    r= row[0]
                    if 'TRACK_ID' in r:
                        plnr= [int(r[10:])]
        
                    elif 'POINT_NUM' in r:
                        duration= [int(r[11:])]
        
                    else:
                        if duration[0] >= duration_limit:
                            out= [float(x) for x in r[:-1].split(' ')]
                            PLlist += plnr + duration + out
                            
        PLlist= np.reshape(np.array(PLlist), (-1, 6)).transpose()
        
        """exclude time steps"""       
        if timestep== 'sixhours':
            PLlist= PLlist[:, PLlist[2]%2 == 1]
            PLlist[2] //=2
            tdiffend= (datetime.date(year+month//12, month%12 +1, 1)-startdate).days*4
            tdiffbegin= (datetime.date(year, month, 1)-startdate).days*4
        
        elif timestep== 'threehours':
            PLlist[2] -= 1 #to let the time count start with 0
            tdiffend= (datetime.date(year+month//12, month%12 +1, 1)-startdate).days*8
            tdiffbegin= (datetime.date(year, month, 1)-startdate).days*8
        
        else: print('wrong timestep')
    
#        print(np.min(PLlist[2]), np.max(PLlist[2]), tdiffbegin, tdiffend)
        
        """exclude data from other months"""
        PLlist= PLlist[:, PLlist[2] < tdiffend]
        PLlist= PLlist[:, PLlist[2] >= tdiffbegin]
        PLlist[2] -= tdiffbegin
              
        [self.PLnumber, duration, self.tPL, self.lon, self.lat, self.vort]= PLlist
        self.lon[self.lon > 180] -= 360 #rotate the longitude data
        self.PLlist= np.array([self.PLnumber, self.tPL, self.lon, self.lat, self.vort])


    def local(self, minlon, maxlon, minlat, maxlat):
        """include only the polar lows in a certain area (box), defined by the coordinates"""
                
        if os.path.isfile(self.filename):                                                                                                                     
            lonminPL=remove_dublicate(self.PLnumber[self.lon> minlon])
            lonmaxPL=remove_dublicate(self.PLnumber[self.lon< maxlon])
            latminPL=remove_dublicate(self.PLnumber[self.lat> minlat])
            latmaxPL=remove_dublicate(self.PLnumber[self.lat< maxlat])            
            intersectPL= list(set(lonminPL).intersection(lonmaxPL).intersection(latminPL).intersection(latmaxPL))                                         
            
            self.PLlist= np.array([self.PLlist[:, x] for x, PLn in enumerate(self.PLnumber) if PLn in intersectPL]).T                                        
#            self.lon= np.array([self.lon[x] for x, PLn in enumerate(self.PLnumber) if PLn in intersectPL])                                        
#            self.lat= np.array([self.lat[x] for x, PLn in enumerate(self.PLnumber) if PLn in intersectPL])                                        
#            self.PLnumber= np.array([self.PLnumber[x] for x, PLn in enumerate(self.PLnumber) if PLn in intersectPL])                                        
            [self.PLnumber, self.tPL, self.lon, self.lat, self.vort]= self.PLlist




class TRACK_list_month2:
    def __init__(self, year, month, duration_limit= 2, approved= True, timestep='sixhours', hemisp='NH', model='ERA'):
        """import the TRACK list for one given months - this is better than the two versions before
        the timestep specifies in which interval the time tPL is taken
        gives the output only in a list: PLlist= [PLnumber, tPL, lon, lat, vort]
        hemisp = ['NH', 'SH']"""
        import datetime

        if model== 'ASR':
            filedir= "/media/"+user+"/1692A00D929FEF8B/PL/PLclim/ASR_TRACKout/ASR_VOR850_"
        elif model== 'ERA':
            if approved == True:
                filedir= "/media/"+user+"/1692A00D929FEF8B/PL/PLclim/ERA_NH_TRACKoutput/ERA_VOR850_"
            else:
                if hemisp == 'NH':
                    filedir= "/media/"+user+"/1692A00D929FEF8B/PL/TRACKoutput/ERA_VOR850_"
                    print('not approved files')
    
                elif hemisp == 'SH':
                    filedir= "/media/"+user+"/1692A00D929FEF8B/PL/PLclim/ERA_SH_TRACKout/ERA_VOR850_"
                    #approving? the SH files were regarded and no real mistake was observed
            
        if month in [1,2,3]: season, startdate= 'A', datetime.date(year, 1, 1)
        elif month in [4,5,6]: season, startdate= 'B', datetime.date(year, 4, 1)
        elif month in [7,8,9]: season, startdate= 'C', datetime.date(year, 7, 1)
        elif month in [10,11,12]: season, startdate= 'D', datetime.date(year, 10, 1)

        if model=='ERA': self.filename= filedir+str(year)+'.'+season+'_'+hemisp+'_tracks'
        elif model=='ASR': self.filename= filedir+str(year)+'_'+season+'_tracks'
        
        PLlist= []
        with open(self.filename, 'r') as f:
            reader = csv.reader(f)
            for rn, row in enumerate(reader):
                if rn >2: #to exclude the first three lines
                    r= row[0]
                    if 'TRACK_ID' in r:
                        plnr= [int(r[10:])]
        
                    elif 'POINT_NUM' in r:
                        duration= [int(r[11:])]
        
                    else:
                        if duration[0] >= duration_limit:
                            out= [float(x) for x in r[:-1].split(' ')]
                            PLlist += plnr + duration + out
                            
        PLlist= np.reshape(np.array(PLlist), (-1, 6)).transpose()
        
        """exclude time steps"""       
        if timestep== 'sixhours':
            PLlist= PLlist[:, PLlist[2]%2 == 1]
            PLlist[2] //=2
            tdiffend= (datetime.date(year+month//12, month%12 +1, 1)-startdate).days*4
            tdiffbegin= (datetime.date(year, month, 1)-startdate).days*4
        
        elif timestep== 'threehours':
            PLlist[2] -= 1 #to let the time count start with 0
            tdiffend= (datetime.date(year+month//12, month%12 +1, 1)-startdate).days*8
            tdiffbegin= (datetime.date(year, month, 1)-startdate).days*8
        
        else: print('wrong timestep')
    
#        print(np.min(PLlist[2]), np.max(PLlist[2]), tdiffbegin, tdiffend)
        
        """exclude data from other months"""
        PLlist= PLlist[:, PLlist[2] < tdiffend]
        PLlist= PLlist[:, PLlist[2] >= tdiffbegin]
        PLlist[2] -= tdiffbegin
        
        PLlist[3][PLlist[3] > 180] -= 360 #rotate the longitude data

        #PLlist= [PLnumber, duration, tPL, lon, lat, vort]
        self.PLlist= np.delete(np.array(PLlist), 1, axis= 0)  #remove the duration
        #self.PLlist= [PLnumber, tPL, lon, lat, vort]


    def local(self, minlon, maxlon, minlat, maxlat):
        """include only the polar lows in a certain area (box), defined by the coordinates"""
                
        if os.path.isfile(self.filename):
            self.PLlist= self.PLlist[:, self.PLlist[2]> minlon]
            self.PLlist= self.PLlist[:, self.PLlist[2]< maxlon]
            self.PLlist= self.PLlist[:, self.PLlist[3]> minlat]
            self.PLlist= self.PLlist[:, self.PLlist[3]< maxlat]



def MatchingTRACKtoSTARS(SPLnrs, model='ERA'):
    """Return a list of the corresponding TRACKcyclone to the STARS PLs with SPLnrs"""
    if model=='ERA':
        filename= "/media/"+user+"/1692A00D929FEF8B/PL/PLclim/MatchList/TRACKmatchSTARS.csv"
    elif model=='ASR':
        filename= "/media/"+user+"/1692A00D929FEF8B/PL/PLclim/MatchList/ASR_TRACKmatchSTARS.csv"
    else: print('wrong model')
    PLlist= np.array(pd.read_table(filename, sep="\t"))
    
    return [PLlist[PLlist[:,0]== SPLnr , 1][0] for SPLnr in SPLnrs]
        
def daynr_month(year, month):
    """this gives how many days the month has"""
    if month in [1,3,5, 7, 8, 10, 12]:
        return 31
    elif month in [4, 6, 9, 11]:
        return 30
    elif month==2 and year%4==0:
        return 29
    elif month==2 and year%4 !=0:
        return 28
