#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 25 22:28:10 2018

@author: pst019
"""

import numpy as np


def flat_list(l):
    """l is list of list which is flattened to a list"""
    return [item for sublist in l for item in sublist]

 
def common_member(a, b): 
    """ check if lists a and b have at-least one element common member"""
    a_set = set(a) 
    b_set = set(b) 
    if (a_set & b_set): 
        return True 
    else: 
        return False


def common_list(a, b): 
    """returns a list of elements that are in lists a and b"""
    a_set = set(a) 
    b_set = set(b) 
    if (a_set & b_set): 
        return list(a_set & b_set)
    else: 
        print("No common elements")  


def FourierFilter2D(lon, lat, vort, T_low= 40, T_up= 100):
    """Filter of waves bigger than T_low and smaller than T_up. 
    T represents how often the wave would go around the equator
    lon and lat have to be equally gridded
    this is tested for 2D data, lon, lat are 1D"""
    
    from scipy.fftpack import fftfreq, fft2, ifft2
    transf= fft2(vort)

    dlat, dlon= lat[0]-lat[1], lon[1]-lon[0]
    
    wny= fftfreq(len(lat), dlat/360) #wavenr in y
    wny= np.tile(wny, (len(lon), 1)).T #dublication of the wavenr in y
    wnx= np.array([list(fftfreq(len(lon), dlon*np.cos(np.deg2rad(la))/360)) for la in range(len(lat))]) # wavenr in x
    
    transf[np.sqrt(wnx**2+wny**2) <T_low] = 0
    transf[np.sqrt(wnx**2+wny**2) >T_up] = 0
    
    """advanced filter. maybe one should take it even more serious with the Gibbs oscillations"""
#    transf[np.sqrt(wnx**2+wny**2) <0.8* T_low] = 0
#    transf[np.sqrt(wnx**2+wny**2) <0.9* T_low] *= .32
#    transf[np.sqrt(wnx**2+wny**2) <1* T_low] *= .625
#    transf[np.sqrt(wnx**2+wny**2) <1.2* T_low] *= .8
#           
#    transf[np.sqrt(wnx**2+wny**2) >1.2* T_up] = 0
#    transf[np.sqrt(wnx**2+wny**2) >1.1* T_up] *= 0.32
#    transf[np.sqrt(wnx**2+wny**2) >1* T_up] *= 0.625
#    transf[np.sqrt(wnx**2+wny**2) >0.8* T_up] *= 0.8 
    
    return np.real(ifft2(transf))
    

def FourierFilter3D(lon, lat, vort, T_low= 40, T_up= 100):
    """Filter of waves bigger than T_low and smaller than T_up. 
    T represents how often the wave would go around the equator
    lon and lat have to be equally gridded and being 1D
    vort is 3D data (time, lat, lon)"""
    
    from scipy.fftpack import fftfreq, fft2, ifft2
    transf= fft2(vort)

    dlat, dlon= lat[0]-lat[1], lon[1]-lon[0]
    ltim= np.shape(vort)[0]
    
    wny= fftfreq(len(lat), dlat/360) #wavenr in y
    wny= np.tile(wny, (len(lon), 1)).T #dublication of the wavenr in y
    wny= np.tile(wny, (ltim, 1, 1))
    wnx= np.array([list(fftfreq(len(lon), dlon*np.cos(np.deg2rad(la))/360)) for la in range(len(lat))]) # wavenr in x
    wnx= np.tile(wnx, (ltim, 1, 1))
    
    transf[np.sqrt(wnx**2+wny**2) <T_low] = 0
    transf[np.sqrt(wnx**2+wny**2) >T_up] = 0

    """advanced filter. maybe one should take it even more serious with the Gibbs oscillations"""
#    transf[np.sqrt(wnx**2+wny**2) <0.8* T_low] = 0
#    transf[np.sqrt(wnx**2+wny**2) <0.9* T_low] *= .32
#    transf[np.sqrt(wnx**2+wny**2) <1* T_low] *= .625
#    transf[np.sqrt(wnx**2+wny**2) <1.2* T_low] *= .8
#           
#    transf[np.sqrt(wnx**2+wny**2) >1.2* T_up] = 0
#    transf[np.sqrt(wnx**2+wny**2) >1.1* T_up] *= 0.32
#    transf[np.sqrt(wnx**2+wny**2) >1* T_up] *= 0.625
#    transf[np.sqrt(wnx**2+wny**2) >0.8* T_up] *= 0.8 
           
    return np.real(ifft2(transf))
    


def FourierFilter2d_equaldist(vort, dist= 2.5, T_low= 40, T_up= 100):
    """filter of waves bigger than T_low and smaller than T_up. 
    T represents how often the wave would go around the equator
    the data has equal grid resolution of dist [km].
    vort is 2d data (xgrid, ygrid)"""

                     
    from scipy.fftpack import fftfreq, fft2, ifft2
    transf= fft2(vort)
           
    wny= fftfreq(np.shape(vort)[1], 2.5/40000) #fftfreq(len(lat), dlat/360) #wavenr in y
    wny= np.tile(wny, (np.shape(vort)[0], 1)) #dublication of the wavenr in y
    wnx= fftfreq(np.shape(vort)[0], 2.5/40000) #np.array([list(fftfreq(len(lon), dlon*np.cos(np.deg2rad(la))/360)) for la in range(len(lat))]) # wavenr in x
    wnx= np.tile(wnx, (np.shape(vort)[1], 1)).T
    
    transf[np.sqrt(wnx**2+wny**2) <T_low] = 0
    transf[np.sqrt(wnx**2+wny**2) >T_up] = 0

    return np.real(ifft2(transf))
    
    
def LocalMax(data, neighborhood_size= (0,12,25), threshold= 6):
    """find local maxima in data in the distance of neighborhoodsize (a constant box),
    neighborhoodsize can be a tuple= (ntim, nlat, nlon) or an integer (then one has the same distance in all directions). A matrix could be chosen, to make it round, since it is a box now.
    above the value of threshold
    returns array with same shape as data, with True for maxima and False else"""
    
    import scipy.ndimage.filters as filters
    data_max = filters.maximum_filter(data, neighborhood_size)
    maxima= (data==data_max)
    maxima[data_max < threshold]=0
    return maxima

    
def LocalMax_dist(data, lat, distlat, threshold= 6):
    """find local maxima in data in a box with distance of distlat (in deg)
    so the neighborhood_size is varying with latitude
    lat - latitude of the dataset
    neighborhoodsize can be a tuple= (ntim, nlat, nlon) or an integer (then one has the same distance in all directions). A matrix could be chosen, to make it round, since it is a box now.
    above the value of threshold
    returns array with same shape as data, with True for maxima and False else
    it gets less efficient with distlat getting big"""    
    import scipy.ndimage.filters as filters

    data_max= np.zeros(np.shape(data))
    distlat= distlat*2 #distance Between maxima gridpoints (since it has a resolution of 0.5 degree)
    for n in range(distlat, len(lat)- distlat):
        distlon= int(distlat/np.cos(np.deg2rad(lat[n])))
        neighborhood_size= (1, distlat, distlon)
    #    print(neighborhood_size)
    #    print(data[:, n-distlat: n+distlat+1, :].shape)
        data_max[:, n] = filters.maximum_filter(data[:, n-distlat: n+distlat+1, :], neighborhood_size, mode= "wrap")[:, distlat]
    
    maxima= (data==data_max)
    maxima[data_max < threshold]=0
    return maxima



        
def remove_dublicate(seq):
    """remove double numbers from list (seq) without changing the order
	seq can also be a list of tuples, e.g [(a,b), (c,d), (c,d)] that can be created by seq= [(x[i],y[i]) for i in range(len(x))]"""
    seen = set()
    seen_add = seen.add
    return [x for x in seq if not (x in seen or seen_add(x))]

    
def remove_dublicate2D(seq):
    seq= [item for sublist in seq for item in sublist]
    return remove_dublicate(seq)



def remove_from_list(item_list, remove_list):
    remove_list= set(remove_list)
    return [e for e in item_list if e not in remove_list]

import itertools
def remove_repeater(inputlist):
    """remove item from list if it occurs several times in a row"""
    return [k for k,_g in itertools.groupby(inputlist)]

import operator
import functools
def split_items(inputlist, splitexpression):
    """splitexpression ' - ' splits ['C', 'C - U', 'U', 'C', '-'] to ['C', 'C', 'U', 'U', 'C', '-']"""
    listoflists= [a.split(' - ') for a in inputlist]
    return functools.reduce(operator.iconcat, listoflists, [])

#def replaceMultiple(mainString, toBeReplaces, newString):
#    """Replace a set of multiple sub strings with a new string in main string."""
#    # Iterate over the strings to be replaced
#    for elem in toBeReplaces :
#        # Check if string is in the main string
#        if elem in mainString :
#            # Replace the string
#            mainString = mainString.replace(elem, newString)
#    
#    return  mainString
    



def EasyTrack(data, lon, lat, distance, threshold, propdistance, latbound=None, lonbound=None):
    """track local maxima in of data
    data - 3D numpy array (time, lon, lat)
    lon, lat - 2D numpy arrays of the longitude and latitude
    threshold - threshold the local maxima have to exceed    
    distance - (integer in gridcells) minimum distance between two local maxima at a given time step
    propdistance - (distance in km) maximum distance in which consequitive points in time are matched
    make distance reasonable larger than propdistance to avoid "doublematches" twice is on the safe side

    latbound [latmin, latmax] - area outside the data is set to zero - so only local max within this are detected
    
    output: local maxima at every time step with their lon and lat and tracked systems are labeled the same
    type list of lists (one per time step) - a bit ugly
    """
    
    import scipy.ndimage.filters as filters
    
#    datad= data[:]
#    if latbound is not None:
#        latdestroy = lat[:] #the np.tile "destroys" latdestroy and therefore a copy is made
#        datad[np.tile(latdestroy, (datad.shape[0],1,1) ) < latbound[0]] = 0        
#        datad[np.tile(latdestroy, (datad.shape[0],1,1) ) > latbound[1]] = 0
#
#    if lonbound is not None:
#        londestroy = lon[:] 
#        datad[np.tile(londestroy, (datad.shape[0],1,1) ) < lonbound[0]] = 0        
#        datad[np.tile(londestroy, (datad.shape[0],1,1) ) > lonbound[1]] = 0
 
    data_max = filters.maximum_filter1d(data, distance, axis= -1)
    data_max = filters.maximum_filter1d(data_max, distance, axis= -2)
    maxima= (data == data_max)
    maxima[data_max < threshold]=0
          
    
    #start outputlists      
    ti= 0
#    print(maxima[ti])
#    print(lat)
    datalist= [list(data[ti][maxima[ti]])]
    lonlist= [list(lon[maxima[ti]] )]
    latlist= [list(lat[maxima[ti]] )]
    maxpointlist= [np.where(maxima[ti])]
    labellist = [list(np.arange(1, len(lonlist[0])+1) )]
    
    maxlabelnr= len(lonlist[0]) +1
    
    for ti in range(1, data.shape[0]):
    
        datalist.append(list(data[ti][maxima[ti]] ))
        lonlist.append(list(lon[maxima[ti]] ))
        latlist.append(list(lat[maxima[ti]] ))
        maxpointlist.append(list(np.where(maxima[ti])))
        
#        print(maxpointlist)
        
        dist=  np.array([[110E3* np.sqrt(np.cos(np.deg2rad(latlist[ti][m])) * (lonlist[ti][m]- lonlist[ti-1][n])**2 + (latlist[ti][m]- latlist[ti-1][n])**2  ) for n in range(len(lonlist[ti-1]))] for m in range(len(lonlist[ti])) ]) #in m
#        print('ti', ti, 'dist', dist)
        mindist= np.min(dist, axis=1)
                                
        labellist1= np.array([ np.array(labellist[ti-1])[dist[n] == mindist[n] ][0] for n in range(len(dist))])
        
        newmaxs= len(labellist1[mindist > propdistance])
        labellist1[mindist > propdistance]= np.arange(maxlabelnr, maxlabelnr+ newmaxs)
        
        if len(list(labellist1)) > len(remove_dublicate(list(labellist1))) :
            print('two matches to same system at time ', ti, labellist1)
        
        labellist.append(list(labellist1))
        maxlabelnr += newmaxs 

    return datalist, lonlist, latlist, labellist #possibly maxpointlist
 

def Data_One_Track(datalist, lonlist, latlist, labellist, cyclnr):
    """ uses the output of EasyTrack to give the data for the local max with cyclnr"""
    timemax= len(datalist) #the amount of timesteps
    
    tcycl= np.array([ti for ti in range(timemax) if len(np.array(datalist[ti])[np.array(labellist[ti])== cyclnr]) == 1])
    datacycl= np.array([ np.array(datalist[ti])[np.array(labellist[ti])== cyclnr][0] for ti in range(timemax) if len(np.array(datalist[ti])[np.array(labellist[ti])== cyclnr]) == 1])
#    maxpointcycl= [ (maxpointlist[ti][0][labellist[ti]== cyclnr], maxpointlist[ti][1][labellist[ti]== cyclnr]) for ti in range(48) if len(datalist[ti][labellist[ti]== cyclnr]) == 1]
    loncycl= np.array([ np.array(lonlist[ti])[np.array(labellist[ti])== cyclnr][0] for ti in range(timemax) if len(np.array(datalist[ti])[np.array(labellist[ti])== cyclnr]) == 1])
    latcycl= np.array([ np.array(latlist[ti])[np.array(labellist[ti])== cyclnr][0] for ti in range(timemax) if len(np.array(datalist[ti])[np.array(labellist[ti])== cyclnr]) == 1])
    return tcycl, datacycl, loncycl, latcycl #possibly also maxpointcycl


def Tracks_in_lat_lon(datalist, lonlist, latlist, labellist, latbound, lonbound, tbound=None):
    """ excludes all tracks that do not pass through a certain area defined by latbound, lonbound
    and writes a new output list of the same style as the input list
    - to be used after EasyTrack"""
    timemax= len(datalist)
    if tbound== None: tbound= [0, timemax]
#    maxlabelnr= np.max(np.array([np.max(labellist[t]) for t in range(timemax)]))
    
    newlabellist= [[] for i in range(timemax)] #(data.shape[0] ) 
    newdatalist= [[] for i in range(timemax)]
    newlatlist= [[] for i in range(timemax)]
    newlonlist= [[] for i in range(timemax)]
    
    for cyclnr in remove_dublicate2D(labellist): # range(1, maxlabelnr+1):
        tcycl, datacycl, loncycl, latcycl= Data_One_Track(datalist, lonlist, latlist, labellist, cyclnr)
        
        #condition for a cyclone to be included - it has to be True for at least one point in time
        condition = np.max( np.logical_and.reduce((latcycl > latbound[0], latcycl < latbound[1], 
                                                   loncycl >lonbound[0], loncycl < lonbound[1],
                                                   tcycl >=tbound[0], tcycl <= tbound[1]                 )) )
    
        if condition == True:
            for ni, ti in enumerate(tcycl):
                newlabellist[ti].append(cyclnr)
                newdatalist[ti].append(datacycl[ni])
                newlatlist[ti].append(latcycl[ni])
                newlonlist[ti].append(loncycl[ni])

    return newdatalist, newlonlist, newlatlist, newlabellist



def Track_OtherMax(data, dist, tcycl, latcycl, datalat, datalon, local='max', diff_grid_data=False, loncycl= False):
    """gets a numpy array of the the local maximum in the data within dist for one track,
    e.g. the maximum wind speed connected to the cyclone center
    local = ('max','min') - specifies if local maximum or minumum is calculated
    data - 3d array
    datalat, the latitude of the data (2d)
    tcycl, loncycl, latcycl - the 1d np.arrays of the the variables of the cyclone
    output: maxdata - the value of the maximum value associated to the cyclone
    maxlon, maxlat - the position of maxdata
    diff_grid_data [False,True] - True if datalat, datalon uses different grid than was used for construction of the TRACK, loncycl must be given"""
    import scipy.ndimage.filters as filters

    maxdata, maxlat, maxlon= np.array([]), np.array([]), np.array([]) 
    for i in range(len(tcycl)):
        if local == 'max':
            if diff_grid_data== False:
                newmax= filters.maximum_filter(data[tcycl[i]], dist)[datalat == latcycl[i]]
            else:
                distfromgrid= (latcycl[i] - datalat)**2 + (np.cos(np.deg2rad(latcycl[i])) * (loncycl[i] - datalon)**2)
                newmax= filters.maximum_filter(data[tcycl[i]], dist)[distfromgrid == np.min(distfromgrid)]
        elif local == 'min':
            if diff_grid_data== False: newmax= filters.minimum_filter(data[tcycl[i]], dist)[datalat == latcycl[i]] 
            else:
                distfromgrid= (latcycl[i] - datalat)**2 + (np.cos(np.deg2rad(latcycl[i])) * (loncycl[i] - datalon)**2)
                newmax= filters.minimum_filter(data[tcycl[i]], dist)[distfromgrid == np.min(distfromgrid)]
        if len(newmax >= 1):
            maxdata= np.append(maxdata, newmax )
            maxlon= np.append(maxlon, datalon[data[tcycl[i]] == newmax][0])
            maxlat= np.append(maxlat, datalat[data[tcycl[i]] == newmax][0])

        else: #there is no local max/min - possibly because latcycl is outside latdata
            maxdata= np.append(maxdata, np.nan )
            maxlon= np.append(maxlon, np.nan)
            maxlat= np.append(maxlat, np.nan)

    return maxdata, maxlon, maxlat



#def radMask3D(data_array, centerpoint , radiusx, radiusy):
#    """uses radMask2D and takes data_array[centerpoint[0]] """
#    return radMask2D(data_array[centerpoint[0]], centerpoint[1:], radiusx, radiusy)
    
def radMask2D(data_shape, centerpoint, radiusx, radiusy,  periody=True):
    """makes a mask in dimensions of the 2D - data_shape
    around the 2D centerpoint with radiusx and radiusy
    it is periodic in the last dimension if continous_y =True"""
    a,b = centerpoint
    nx,ny = data_shape
    
#    print(radiusx, radiusy, nx, ny)
    x,y = np.ogrid[-a:nx-a,-b:ny-b]
#    print(x,y)
    if periody==True:
        y[y > ny/2] -= ny
#    print(y)

    mask = x**2/(radiusx**2) + y**2/radiusy**2 <= 1
    return mask   

    
def radMaskPixel(centerpoint, radius, array, periody= True):
    """makes a radial mask in the shape of array around the centerpoint with radius (in pixel)
    can be periodic in second dimension"""
    radiusx= radius
    radiusy= radius
    a,b = centerpoint
    nx,ny = array.shape
    
    x,y = np.ogrid[-a:nx-a,-b:ny-b]

    if periody == True:
        y[y > ny/2] -= ny #periodic in second dimension

    mask = x**2/(radiusx**2) + y**2/radiusy**2 <= 1
    return mask



def MeanAroundPoint(centerpoint, radius, data, periody=False):
    """calculates the mean of data [(t,) x,y] around centerpoint [(t), x,y] within radius [pixels]
    data can be of shape x,y or t,x,y"""
    
    if len(np.shape(data)) == 2:
        return np.mean(data[radMaskPixel(centerpoint, radius, data, periody= periody)])
    
    else:
        lent= np.shape(data)[0]
        mean= np.zeros(lent)
        
        for t in range(lent):
#            print(t)
#            print(data[t])
            mean[t]= np.mean(data[t][radMaskPixel(centerpoint[t], radius, data[t], periody= periody)])
        
        return mean
    
        

def radMask(centerpoint,radius,array, dx, dy):
    """makes a radial (eliptic) mask in the shape of array around centerpoint with radius where
    dx - distance between latitudes
    dy - distance between longitudes
    can be periodic in second dimension"""
    radiusx= radius/dx
    radiusy= radius/dy
    a,b = centerpoint
    nx,ny = array.shape
    
#    print(radiusx, radiusy, nx, ny)
    x,y = np.ogrid[-a:nx-a,-b:ny-b]
#    print(x,y)
    y[y > ny/2] -= ny
#    print(y)

    mask = x**2/(radiusx**2) + y**2/radiusy**2 <= 1
    return mask    






def CircularFilter(data, radius, dx, lat):
    """calculate the mean of data= array((lat, lon))
    of the cells within radius
    dx= distance between latitudes, lat= array of the latitudes to calculate dy"""
    
    meanAll=np.zeros(np.shape(data))
    
    print('start (the Circular filter is quite unefficient)')
    for x in range(np.shape(data)[0]):
        for y in range(np.shape(data)[1]):
            centerMask=(x,y)
            dy=dx*np.cos(np.deg2rad(lat[x]))
            mask=radMask(centerMask,radius,data, dx, dy)
    
            #get the mean
            meanAll[x,y]=np.mean(data[mask])
    
#    print('end')
    return meanAll


def CircularFilter_samedist(data, radius, periody=False):
    """calculate the mean of data= array((lat, lon))
    of the cells within radius
    dx= distance between latitudes, lat= array of the latitudes to calculate dy"""
    
    meanAll=np.zeros(np.shape(data))
    
#    print('start (the Circular filter is quite unefficient)')
    for x in range(np.shape(data)[0]):
        for y in range(np.shape(data)[1]):
            centerpoint=(x,y)
            mask=radMaskPixel(centerpoint,radius, data, periody= periody)
    
            #get the mean
            meanAll[x,y]=np.mean(data[mask])
    
#    print('end')
    return meanAll





def unit_vector(vector):
    """ Returns the unit vector of the vector.  """
    return vector / np.linalg.norm(vector)

def angle_between(v1, v2):
    """ Returns the angle in degree between vectors 'v1' and 'v2'::
    """
    v1_u = unit_vector(v1)
    v2_u = unit_vector(v2)
    return np.rad2deg(np.arccos(np.clip(np.dot(v1_u, v2_u), -1.0, 1.0)))


def distance(latlon1, latlon2):
    """returns the distance between the two points latlon1, latlon2, each a tuple of (lat, lon)
    one of the "points" can also be an array"""
        
    R = 6373.0 # approximate radius of earth in km
    
    lat1 = np.radians(latlon1[0])
    lon1 = np.radians(latlon1[1])
    lat2 = np.radians(latlon2[0])
    lon2 = np.radians(latlon2[1])
    
    dlon = lon2 - lon1
    dlat = lat2 - lat1
    
    a = np.sin(dlat / 2)**2 + np.cos(lat1) * np.cos(lat2) * np.sin(dlon / 2)**2
    c = 2 * np.arctan2(np.sqrt(a), np.sqrt(1 - a))
    
    return R * c



def grad_x(var, lat, lon):
    """calculates the gradient in x-direction of the variable with dimensions(lat, lon) with lat starting from the pole. (The ERA-5 data that I have downloaded)"""
    import xarray as xr
    if type(lat) == xr.core.dataarray.DataArray:
        lat= lat.values
        lon= lon.values
    
    diff_deg= lon[1] - lon[0]
    dx= diff_deg *110E3 * np.cos(np.deg2rad(lat))
    dx= np.tile(dx, (len(lon), 1)).T
    return np.gradient(var, axis= 1)/dx
    
    
def grad_y(var, lat, lon):
    import xarray as xr
    if type(lat) == xr.core.dataarray.DataArray:
        lat= lat.values
        lon= lon.values
        
    dy= (lat[1] - lat[0]) *110E3 
    return np.gradient(var, dy, axis= 0)




def value_in_rad(field, lat, lon, center_lat, center_lon, distance, type= 'mean', intype='pd'):
    """calculuates the 'type'=[mean, max, min, median]' around the center_lat, center_lon, with radius= distance for the field with coordinates (lat, lon)
    intype gives the type of field which can be 'pd' or 'np'"""
    
    if len(np.shape(lat))== 1: #create 2d lat, lon
        lon, lat= np.meshgrid(lon, lat)
    
#    print( center_lat )
#    if type(center_lat) == pd.core.series.Series: #convert the center location from pandas series to a value
#        center_lat, center_lon= center_lat.values[0], center_lon.values[0]
        
    
    dist_now= 110* np.sqrt( (center_lat- lat)**2+ (np.cos(np.deg2rad(center_lat))* (center_lon- lon))**2)
    
   
#    if type(field) is not np.ndarray: field= field.values
#    if not isinstance(type(field), np.ndarray): field= field.values
    if intype=='pd': field= field.values
    
    if type== 'mean': value= np.mean(field[dist_now < distance])
    elif type== 'max': value= np.max(field[dist_now < distance])
    elif type== 'min': value= np.min(field[dist_now < distance])
    elif type== 'med': value= np.median(field[dist_now < distance])

    return value


