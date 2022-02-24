#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 11 11:11:16 2017

@author: pst019
"""

"""first run TRACKmatchSTARS_stats3 and TRACKsatisfyCrit
then makes histograms out of the statistics"""

import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator

fignr= 8

plt.figure(fignr)
plt.clf()

ax= plt.figure(fignr).gca()

lim= None
var='large' #whether a variable is larger or smaller for PLs than for all cyclones

#data=np.array(stattheta_e700_max)
#data_track=np.array(stat_track_theta_e700_max)
#name= r'$\theta_{e, SST}$ - $\theta_{e, 700}$ [K]' #'stabtheta_e700'
#lim= [-25, 25]
##
#data=np.array(stattheta_e500)
#data_track=np.array(stat_track_theta_e500)
#name= '$\theta_{e, SST}$ - $\theta_{e, 500}$ [K]' #stabtheta_e500'

#data=np.array(stattheta700)
#data_track=np.array(stat_track_theta700)
#name= '$\theta_{SST}$ - $\theta_{700}$ [K]' #'stabtheta700'

#data=np.array(stattheta500)
#data_track=np.array(stat_track_theta500)
#name= r'$\theta_{SST}$ - $\theta_{500}$ [K]' #'stabtheta500'
#lim = [-40, 7]

#data=np.array(stattheta500_max)
#data_track=np.array(stat_track_theta500_max)
#name= r'$\theta_{SST}$ - $\theta_{500}$ [K]' #'stabtheta500'
#lim = [-40, 7]
             
#data=np.array(stattheta850_2)
#data_track=np.array(stat_track_theta850_2)
#name= r'$\theta_{SST}$ - $\theta_{850}$ [K]'  
#lim= [-22, 15]
                
#data=np.array(stattheta_e850_max_1)
#data_track=np.array(stat_track_theta_e850_max_1)
#name= r'$\theta_{e,SST}$ - $\theta_{e,850}$ [K]'                
#lim=[-25,26]
                
data=np.array(statstab500)
data_track=np.array(stat_track_stab500)
name= r'SST - T$_{500}$ [K]'
#
#data=np.array(statstab700)
#data_track=np.array(stat_track_stab700)
#name= 'stab700'
##
#data=np.array(statt_500)
#data_track=np.array(stat_track_t_500)
#name= 'temp_500'
##
#data=np.array(statt_700)
#data_track=np.array(stat_track_t_700)
#name= 'temp_700'
#
#data=np.array(statattheta_e700)
#data_track=np.array(stat_track_attheta_e700)
#name= 'theta_e700'
#
#data=np.array(statattheta_e850)
#data_track=np.array(stat_track_attheta_e850)
#name= 'theta_e850'
#
#data=np.array(statSST)
#data_track=np.array(stat_track_SST)
#name= 'temp_SST'
##
#data=np.array(statMCAO_Kolstad)
#data_track=np.array(stat_track_MCAO_Kolstad)
#name= 'MCAO index' #'stabtheta500'
#
#data=np.array(statMCAO_Kolstad_0)*1e3
#data_track=np.array(stat_track_MCAO_Kolstad_0)*1e3
#name= 'MCAO index [K/1000\,hPa]' #K/hPa'stabtheta500'
##
#data=np.array(statMCAO_Kolstad_0_500)*1e3
#data_track=np.array(stat_track_MCAO_Kolstad_0_500)*1e3
#name= 'MCAO index 500 [K/1000\,hPa]' #K/hPa'stabtheta500'
# 
#data=np.array(statp_PV_3)/100
#data_track=np.array(stat_track_p_PV_3)/100
#name='$p_{tr}$ [hPa]'
#lim= [100, 800]            
##
#data=np.array(statU10_2)
#data_track=np.array(stat_track_U10_2)
#name=r'$U_{10m}$ [m s$^{-1}$]'
#lim=[0, 40]

#data=np.array(statvort)
#data_track=np.array(stat_track_vort)
#name= r'$\xi_{f, 850}$ [$10^{-5} $s$^{-1}$]'
#lim= [2, 17]
#
#data=np.array(statWater)
#data_track=np.array(stat_track_Water)
#name= 'Total column water [kg m$^{-2}$]'
#lim = [0, 50]
#var= 'small'
#
#data=np.array(statUPVnorth)
#data_track=np.array(stat_track_UPVnorth)
#name=r'$U_{tr,p}$ [m s$^{-1}$]'
#var= 'small'
#
#data=np.array(statUPVnorthM2p5)
#data_track=np.array(stat_track_UPVnorthM2p5)
#name='UPVnorth minus 2.5'
#
#data=np.array(statUPL500north)
#data_track=np.array(stat_track_UPL500north)
#name='UPL 500 north'

#data= np.array(statgradtheta_e850north)
#data_track= np.array(stat_track_gradtheta_e850north)
#name= 'grad theta_e north'

#data=np.array(statUPV)
#data_track=np.array(stat_track_UPV)
#name='UPV'

#data=np.array(statthetaPV_3)
#data_track= np.array(stat_track_thetaPV_3)
#name=r'$\theta_{tr}$ [K]'
#var='small'
##lim = [275, 385]
#
#data=np.array(statthetaSSTmPV_3)
#data_track= np.array(stat_track_thetaSSTmPV_3)
#name=r'$\theta_{SST} - \theta_{tr}$[K]'
#var='large'
#lim=[-100, 5]

#data=np.array(statPBH)
#data_track=np.array(stat_track_PBH)
#name='PBH'
#
#data=np.array(statMSLP)
#data_track=np.array(stat_track_MSLP)
#name='MSLP'

#data= np.array(statgradtheta_e5)
#data_track=np.array(stat_track_gradtheta_e5)
#name=r'$\nabla \theta_{e,850}$ [$10^{-2}$ K km$^{-1}$]'
#lim= [0,25]
#var= 'small'

#data= np.array(statpdiff1)*-100
#data_track=np.array(stat_track_pdiff1)*-100
#name=r'$\overline{SLP}$ - SLP [Pa]'
#lim= [-50, 250]

data=np.array(statCAPE_max_2)
data_track=np.array(stat_track_CAPE_max_2)
name=r'CAPE [J kg$^{-1}$]'
lim= [0, 1200]
             
SPL = np.array(statSPL)
TPL= np.array(stat_track_TPL)

nSPL= len(remove_dublicate(SPL))
nTPL= len(remove_dublicate(TPL))

print('Before exclusion:  STARS PLs: ', nSPL, 'TRACK cyclones: ', nTPL)


exclusion= False
exclusion3=False

exclusion=True
exclusion3=True

percentile= 10

"""other exclusion"""
"""a)criteria"""
vortcrit= np.float64(5.04)
vortlist= statvort > vortcrit 
vortlist_track= stat_track_vort > vortcrit 

UPVnorthcrit= np.float64(31.29) #minimum value in the lifetime of a polar low                        
UPVnorthlist= statUPVnorth < UPVnorthcrit 
UPVnorthlist_track= stat_track_UPVnorth < UPVnorthcrit

stabcrit= np.float64(43.0) #(39.6)
stablist= statstab500 > stabcrit 
stablist_track= stat_track_stab500 > stabcrit 

stabtheta500crit= np.float64(-9.403) #(-4.422) # #using the max
stablist= stattheta500_max > stabtheta500crit 
stablist_track= stat_track_theta500_max > stabtheta500crit 

pdiffcrit =np.float64(-0.409)
pdifflist= statpdiff1 < pdiffcrit 
pdifflist_track= stat_track_pdiff1 <pdiffcrit

#stabtheta500crit= np.float64(-11.2) #using the mean
#stablist= stattheta500 > stabtheta500crit 
#stablist_track= stat_track_theta500 > stabtheta500crit 

#stabtheta_e500crit= np.float64(-12.4)
#stablist= stattheta_e500 > stabtheta_e500crit 
#stablist_track= stat_track_theta_e500 > stabtheta_e500crit 

U10crit =np.float64(15.) #(13.317)
U10list= statU10_2 > U10crit 
U10list_track= stat_track_U10_2 > U10crit

"""b) exclusion lists"""
#exclist= pdifflist
#exclist_track= pdifflist_track

#exclist= stablist
#exclist_track= stablist_track

#exclist= UPVnorthlist
#exclist_track= UPVnorthlist_track

#2 crits - UPVnorth, stab
#exclist_track= np.array([a and b for a, b in zip(UPVnorthlist_track, stablist_track)])
#exclist= np.array([a and b for a, b in zip(UPVnorthlist, stablist)])
#
##2 crits - pdiff, stab
#exclist_track= np.array([a and b for a, b in zip(pdifflist_track, stablist_track)])
#exclist= np.array([a and b for a, b in zip(pdifflist, stablist)])

#2 crits - pdiff, UPVnorth
exclist_track= np.array([a and b for a, b in zip(pdifflist_track, UPVnorthlist_track)])
exclist= np.array([a and b for a, b in zip(pdifflist, UPVnorthlist)])

#2crits - vort, stab
#exclist=np.array([a and b for a, b in zip(vortlist, stablist)])
#exclist_track= np.array([a and b for a, b in zip(vortlist_track, stablist_track)])

#3 crits - UPVnorth, stab, wind
#exclist_track= np.array([a and b and c for a, b, c in zip(UPVnorthlist_track, U10list_track, stablist_track)])
#exclist= np.array([a and b and c for a, b, c in zip(UPVnorthlist, U10list, stablist)])


##2 crits - vort, UPVnorth
#exclist_track= np.array([a and b for a, b in zip(vortlist_track, UPVnorthlist_track)])
#exclist= np.array([a and b for a, b in zip(vortlist, UPVnorthlist)])

#2 crits - vort, UPVnorth
#exclist_track= np.array([a and b for a, b in zip(vortlist_track, UPVnorthlist_track)])
#exclist= np.array([a and b for a, b in zip(vortlist, UPVnorthlist)])

#3 crits - vort, UPVnorth, wind
#exclist_track= np.array([a and b and c for a, b, c in zip(vortlist_track, UPVnorthlist_track, U10list_track)])
#exclist= np.array([a and b and c for a, b, c in zip(vortlist, UPVnorthlist, U10list)])
#
##3 crits - vort, stab, wind
#exclist_track= np.array([a and b and c for a, b, c in zip(vortlist_track, U10list_track, stablist_track)])
#exclist= np.array([a and b and c for a, b, c in zip(vortlist, U10list, stablist)])

#4 crits - vort, UPVnorth, U10, stab
#exclist_track= np.array([a and b and c and d for a, b, c, d in zip(vortlist_track, UPVnorthlist_track, U10list_track, stablist_track)])
#exclist= np.array([a and b and c and d for a, b, c, d in zip(vortlist, UPVnorthlist, U10list, stablist)])

"""3crits """
if exclusion3== True:
    #vort, UPVnorth, stab
#    exclist=np.array([a and b and c for a, b, c in zip(vortlist, UPVnorthlist, stablist)])
#    exclist_track= np.array([a and b and c for a, b, c in zip(vortlist_track, UPVnorthlist_track, stablist_track)])
    
    #pdiff, UPVnorth, stab
    exclist=np.array([a and b and c for a, b, c in zip(pdifflist, UPVnorthlist, stablist)])
    exclist_track= np.array([a and b and c for a, b, c in zip(pdifflist_track, UPVnorthlist_track, stablist_track)])


"""c) exclusion of data"""
if exclusion == True:
    data_track= data_track[exclist_track]
    data= data[exclist]
    SPL= SPL[exclist]
    TPL= TPL[exclist_track]

nSPL= len(remove_dublicate(SPL))
nTPL= len(remove_dublicate(TPL))
print('After exclusion:  STARS PLs: ', nSPL, 'TRACK cyclones: ', nTPL)




"""get the maximum or minimum of the cyclone
this should often be excluded """
#print('long to max/min: the statistic is for the maximum/minimum value during the lifetime of the cyclone/ polar low')
data_max, data_track_max = [], []
for n in remove_dublicate(SPL):
    if var=='small' or 'UPV' in name or 'temp_' in name or 'theta_e' in name or 'MSLP' in name or 'UPL 500 north' in name: data_max += [np.min(data[SPL== n])]
    else: data_max += [np.max(data[SPL== n])]

for n in remove_dublicate(TPL):    
    if var=='small' or'UPV' in name or 'temp_' in name or 'theta_e' in name or 'MSLP' in name  or 'UPL 500 north' in name: data_track_max += [np.min(data_track[TPL==n])]
    else: data_track_max += [np.max(data_track[TPL==n])]

data= data_max
data_track= data_track_max
SPL= np.array(remove_dublicate(SPL))
TPL= np.array(remove_dublicate(TPL))




"""criteria of reduced data"""
if exclusion == False:
    datamean= np.mean(data)
    stdfactor= 1
    if var=='small' or 'UPV' in name or 'temp_' in name or 'theta_e' in name or 'MSLP' in name or 'UPL 500 north' in name:
        crit= np.percentile(data, 100-percentile)
        stdcrit=datamean+ stdfactor* np.std(data)
        
    else:
        crit= np.percentile(data, percentile)
        stdcrit=datamean- stdfactor* np.std(data)
        
#    print('STARS mean', datamean, 'std crit', stdcrit)


else:
    print('old criteria are used')


#crit=np.float64(43.0)


"""plot the histograms"""
if exclusion == False:
    bindist= (np.max(data_track)-np.min(data_track))//80 +1 #80 was 40 before
    if '$\overline{SLP}$' in name: bindist = 15       
    mindata= np.round(np.min([np.min(data), np.min(data_track)]) -0.5)
    offset= mindata%bindist #this is done such that the first bin would always start at 0
    mindata -= offset

bins= np.arange(mindata, np.round(np.max(data)+0.5), bindist)
bins_track= np.arange(mindata, np.round(np.max(data_track)+0.5), bindist)

weights = np.ones_like(data_track)/(len(data_track)*bindist)
plt.hist(np.array(data_track), bins=bins_track, alpha= 0.7, weights=weights*100, label='Cyclones')

weights= np.ones_like(data)/(len(data)*bindist)
plt.hist(np.array(data), bins=bins, alpha= 0.7, weights=weights*100,  label='Polar lows')

plt.legend(fontsize= 16)

line= plt.plot([datamean], [0], 'o')[0]
line.set_clip_on(False) #this way the circle is going over the axis
line= plt.plot([crit], [0], 'o')[0]
line.set_clip_on(False)

plt.ylabel('Normalized frequency [%]', size= 16)

if name == 'stab500':
    plt.xlabel('SST - T$_{500}$ [K]')

elif name == 'stab700':
    plt.xlabel('SST - T$_{700}$ [K]')
#    plt.xlim([0, 40])
    
elif name == 'stabtheta_e500':
    plt.xlabel('SST - $\theta_{e, 500}$ [K]')

elif name == 'stabmax2theta_e500':
    plt.xlabel('max (SST - $\theta_{e, 500})$ [K]')

elif name == 'stabtheta_e700':
    plt.xlabel('SST - $\theta_{e, 700}$ [K]')

elif name == 'stabtheta500':
    plt.xlabel('SST - $\theta_{500}$ [K]')

elif name == 'stabtheta700':
    plt.xlabel('SST - $\theta_{700}$ [K]')

elif name == 'deltheta_e500':
    plt.xlabel('$\theta_{e, 1000}$ - $\theta_{e, 500}$ [K]')
    
elif name == 'temp_500':
    plt.xlabel('500 hPa temperature [K]')

elif name == 'temp_700':
    plt.xlabel('700 hPa temperature [K]')
    
elif name == 'temp_SST':
    plt.xlabel('SST')
    
elif name == 'UPV':
    plt.xlabel('Maximum tropopause wind in 100km radius [m/s]')

    
else:
    plt.xlabel(name, size= 16)
if lim != None: plt.xlim(lim)
plt.tick_params(labelsize= 16)

ax.yaxis.set_major_locator(MaxNLocator(integer=True))

plt.tight_layout()


#plt.savefig('/home/'+user+'/home/Polar_Low/Documents/ERA-PLs/Graphs/'+name+'.png', dpi= 70, bbox_inches='tight')

print(name, 'crit', crit)


#print('STARS PLs: ', nSPL, 'TRACK cyclones: ', nTPL)

if var=='small' or 'UPV' in name or 'temp_' in name or 'theta_e' in name or 'MSLP' in name or 'water' in name or 'UPL 500 north' in name:
#    print('exclusion of TRACK points: ', len(np.where(data_track > crit)[0])/len(data_track))
#    print('exclusion of STARS points: ', len(np.where(data > crit)[0])/len(data))
    
    trackexcl= (nTPL- len(remove_dublicate(TPL[data_track < crit])) )/nTPL
    print('exclusion of TRACK PLs: ', trackexcl, 'nr removed:', nTPL- len(remove_dublicate(TPL[data_track < crit])))    
    print('exclusion of STARS PLs: ', (nSPL- len(remove_dublicate(SPL[data < crit])) )/nSPL, 'nr removed:',nSPL- len(remove_dublicate(SPL[data < crit])) )    
#    , SPLnr[data > crit]
    
else:
#    print('exclusion of TRACK points: ', len(np.where(data_track < crit)[0])/len(data_track))
#    print('exclusion of STARS points: ', len(np.where(data < crit)[0])/len(data))

    trackexcl= (nTPL- len(remove_dublicate(TPL[data_track > crit])) )/nTPL
    print('exclusion of TRACK PLs: ', trackexcl, 'nr removed:',nTPL- len(remove_dublicate(TPL[data_track > crit])))    
    print('exclusion of STARS PLs: ', (nSPL- len(remove_dublicate(SPL[data > crit])) )/nSPL,'nr removed:', nSPL- len(remove_dublicate(SPL[data > crit])))
#    , SPLnr[data < crit]

#plt.title(str(np.round(trackexcl, 3)))