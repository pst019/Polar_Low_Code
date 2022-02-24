#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 11 11:11:16 2017

@author: pst019
"""

"""first run TRACKmatchSTARS_stats3 and TRACKsatisfyCrit
then makes histograms out of the statistics"""

from matplotlib.ticker import MaxNLocator

fignr= 6

plt.figure(fignr)
plt.clf()
ax= plt.figure(fignr).gca()

lim= None
var='large'


data=np.array(stattheta700)
data_track=np.array(stat_track_theta700)
name= '$\Theta_{SST}$ - $\Theta_{700}$ [K]' #'stabtheta700'

#data=np.array(stattheta500_3)
#data_track=np.array(stat_track_theta500_3)
#name= 'stabtheta500'
#
data=np.array(stattheta500_max)
data_track=np.array(stat_track_theta500_max)
name= r'$\theta_{SST}$ - $\theta_{500}$ [K]' #'stabtheta500'
lim = [-40, 7]

#data=np.array(statstab500_max)
#data_track=np.array(stat_track_stab500_max)
#name= 'stab500'
#
#data=np.array(statstab700_max)
#data_track=np.array(stat_track_stab700_max)
#name= 'stab700'
###
#data=np.array(statt_500)
#data_track=np.array(stat_track_t_500)
#name= 'temp_500'
#var='small'
###
#data=np.array(statt_700)
#data_track=np.array(stat_track_t_700)
#name= 'temp_700'
#var='small'

##
#data=np.array(statSST)
#data_track=np.array(stat_track_SST)
#name= 'temp_SST'
#var='small'


#data=np.array(statU10_2)
#data_track=np.array(stat_track_U10_2)
#name='U10'
#
#data=np.array(statvort)
#data_track=np.array( stat_track_vort)
#name= 'vort'
###
#data=np.array(statU500north)
#data_track=np.array(stat_track_U500north)
#name=r'$U_{500,p}$ [m s$^{-1}$]'
#var='small'
#lim= [0, 70]

#data=np.array(statU500)
#data_track=np.array(stat_track_U500)
#name='U500'
#var='small'

#data= np.array(statpdiff3)*10
#data_track= np.array(stat_track_pdiff3) *10
#name='$\overline{SLP}$ - SLP [Pa]'
#lim= [-500, 1500]

percentile = 50                       
                        
exclusion= False
exclusion3= False
#
exclusion=True
exclusion3= True

SPL = np.array(statSPL)
TPL= np.array(stat_track_TPL)

nSPL= len(remove_dublicate(SPL))
nTPL= len(remove_dublicate(TPL))

print('Before exclusion:  STARS PLs: ', nSPL, 'TRACK cyclones: ', nTPL)


"""other exclusion"""
"""a)criteria"""
vortcrit= np.float64(4.0)# (3.61)
vortlist= statvort > vortcrit 
vortlist_track= stat_track_vort > vortcrit 

U500northcrit= np.float64(29.58) #(31.46) #(29.8) # #minimum value in the lifetime of a polar low                     
U500northlist= statU500north < U500northcrit 
U500northlist_track= stat_track_U500north < U500northcrit

#stabcrit= np.float64(39.6)
#stablist= statstab500 > stabcrit 
#stablist_track= stat_track_stab500 > stabcrit 


stabtheta500crit= np.float64(-8.514) # for the normal systems (10th percentile)
stabtheta500crit= np.float64(-3.978) # for the most intense systems (50th percentile)

stablist= stattheta500_max > stabtheta500crit 
stablist_track= stat_track_theta500_max > stabtheta500crit 

#stabtheta500crit= np.float64(-13.0)
#stablist= stattheta500 > stabtheta500crit 
#stablist_track= stat_track_theta500 > stabtheta500crit 

U10crit =np.float64(17.5)
U10list= statU10_2 > U10crit 
U10list_track= stat_track_U10_2 > U10crit

#this should be positive now!!
pdiffcrit= np.float64(23.82) #for the normal systems(10th percentile)
pdiffcrit= np.float64(50.71) #for the most intense systems (50th percentile)
pdifflist= statpdiff3 > pdiffcrit
pdifflist_track= stat_track_pdiff3 > pdiffcrit

#"""b) exclusion lists"""
#1crit
#exclist= pdifflist
#exclist_track= pdifflist_track

#1crit
exclist= stablist
exclist_track= stablist_track

#2crits - pdiff, stab
exclist=np.array([a and b for a, b in zip(pdifflist, stablist)])
exclist_track= np.array([a and b for a, b in zip(pdifflist_track, stablist_track)])

#2crits - pdiff, U500north
exclist=np.array([a and b for a, b in zip(pdifflist, U500northlist)])
exclist_track= np.array([a and b for a, b in zip(pdifflist_track, U500northlist_track)])

#2crits - stab, vort
#exclist=np.array([a and b for a, b in zip(vortlist, stablist)])
#exclist_track= np.array([a and b for a, b in zip(vortlist_track, stablist_track)])

#2 crits - U500north, stab
#exclist_track= np.array([a and b for a, b in zip(U500northlist_track, stablist_track)])
#exclist= np.array([a and b for a, b in zip(U500northlist, stablist)])


#"""3crits - vort, U500north, stab
#exclist=np.array([a and b and c for a, b, c in zip(vortlist, U500northlist, stablist)])
#exclist_track= np.array([a and b and c for a, b, c in zip(vortlist_track, U500northlist_track, stablist_track)])

if exclusion3== True:
    #3 crits - U500north, stab, pdiff
    exclist_track= np.array([a and b and c for a, b, c in zip(U500northlist_track, pdifflist_track, stablist_track)])
    exclist= np.array([a and b and c for a, b, c in zip(U500northlist, pdifflist, stablist)])

#3 crits - U500north, stab, wind
#exclist_track= np.array([a and b and c for a, b, c in zip(U500northlist_track, U10list_track, stablist_track)])
#exclist= np.array([a and b and c for a, b, c in zip(U500northlist, U10list, stablist)])


#2 crits - vort, U500north
#exclist_track= np.array([a and b for a, b in zip(vortlist_track, U500northlist_track)])
#exclist= np.array([a and b for a, b in zip(vortlist, U500northlist)])

#2 crits - vort, U500north
#exclist_track= np.array([a and b for a, b in zip(vortlist_track, U500northlist_track)])
#exclist= np.array([a and b for a, b in zip(vortlist, U500northlist)])

#3 crits - vort, U500north, wind
#exclist_track= np.array([a and b and c for a, b, c in zip(vortlist_track, U500northlist_track, U10list_track)])
#exclist= np.array([a and b and c for a, b, c in zip(vortlist, U500northlist, U10list)])
#
##3 crits - vort, stab, wind
#exclist_track= np.array([a and b and c for a, b, c in zip(vortlist_track, U10list_track, stablist_track)])
#exclist= np.array([a and b and c for a, b, c in zip(vortlist, U10list, stablist)])

#4 crits
#exclist_track= np.array([a and b and c and d for a, b, c, d in zip(vortlist_track, U500northlist_track, U10list_track, stablist_track)])
#exclist= np.array([a and b and c and d for a, b, c, d in zip(vortlist, U500northlist, U10list, stablist)])

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
#print('long to max/min: the statistic is for the maximum value during the cyclone/ polar low')
data_max, data_track_max = [], []
for n in remove_dublicate(SPL):
    if var=='small': data_max += [np.min(data[SPL== n])]
    else: data_max += [np.max(data[SPL== n])]

for n in remove_dublicate(TPL):    
    if var=='small': data_track_max += [np.min(data_track[TPL==n])]
    else: data_track_max += [np.max(data_track[TPL==n])]

data= data_max
data_track= data_track_max
SPL= np.array(remove_dublicate(SPL))
TPL= np.array(remove_dublicate(TPL))



"""criteria of reduced data"""
if exclusion == False:
    stdfactor= 1
    if var=='small':
        crit= np.percentile(data, 100 - percentile)
#        crit= np.percentile(data, 95)
        meandata=np.mean(data)    
        stdcrit= np.mean(data) + stdfactor * np.std(data)
        
    else:
        crit= np.percentile(data, percentile)
#        crit= np.percentile(data, 5)
        meandata= np.mean(data)    
        stdcrit= np.mean(data) - stdfactor * np.std(data)

#print('STARS mean', meandata, 'STARS std', np.std(data), '10% perc', crit2,'5% perc', crit,  'std crit', stdcrit)

#crit= np.float64(261.6) #temp_700
#crit = np.float64(243.9)
#crit=np.float64(U10crit)

"""plot the histograms"""
if exclusion == False:
    bindist= (np.max(data)-np.min(data))//20 +1
    if '$\overline{SLP}$' in name: bindist = 100             
    mindata= np.round(np.min([np.min(data), np.min(data_track)]) -0.5)
    bins= np.arange(mindata, np.round(np.max(data)+0.5), bindist)
    bins_track= np.arange(mindata, np.round(np.max(data_track)+0.5), bindist)

weights = np.ones_like(data_track)/(len(data_track)*bindist)
plt.hist(data_track, bins=bins_track, alpha= 0.7, weights=weights*100, label='Cyclones')


weights= np.ones_like(data)/(len(data)*bindist)
plt.hist(data, bins=bins, alpha= 0.7, weights=weights*100,  label='Polar lows')

plt.legend(fontsize= 16)

line= plt.plot([meandata], [0], 'o')[0]
line.set_clip_on(False)

#plt.plot([crit2], [0], 'o')
line= plt.plot([crit], [0], 'o')[0]
line.set_clip_on(False)

plt.ylabel('Normalized frequency [%]', size= 16)

#crit = crit2


if name == 'vort':
    plt.xlabel('Vorticity [$10^{-5} $s$^{-1}$]')
    plt.xlim([2, 18])

elif name == 'stab500':
    plt.xlabel('SST - T$_{500}$ [K]')

elif name == 'stab700':
    plt.xlabel('SST - T$_{700}$ [K]')
#    plt.xlim([0, 40])
    
elif name == 'stabtheta_e500':
    plt.xlabel('SST - $\Theta_{e, 500}$ [K]')

elif name == 'stabmax2theta_e500':
    plt.xlabel('max (SST - $\Theta_{e, 500})$ [K]')

elif name == 'stabtheta_e700':
    plt.xlabel('SST - $\Theta_{e, 700}$ [K]')

elif name == 'stabtheta500':
    plt.xlabel('$\Theta_{SST}$ - $\Theta_{500}$ [K]')
    plt.xlim([-45, 5])

elif name == 'stabtheta700':
    plt.xlabel('SST - $\Theta_{700}$ [K]')

elif name == 'deltheta_e500':
    plt.xlabel('$\Theta_{e, 1000}$ - $\Theta_{e, 500}$ [K]')
    
elif name == 'temp_500':
    plt.xlabel('500 hPa temperature [K]')

elif name == 'temp_700':
    plt.xlabel('700 hPa temperature [K]')
    
elif name == 'U10':
    plt.xlabel('10m wind speed in 250km radius [m/s]')
    plt.xlim([0, 40])
    
elif name == 'U500':
    plt.xlabel('Maximum tropopause wind in 100km radius [m/s]')

elif name == 'water':
    plt.xlabel('Total column water [kg/m$^2$]')

else:
    plt.xlabel(name, size= 16)


if lim != None: plt.xlim(lim)
plt.tick_params(labelsize= 16)
ax.yaxis.set_major_locator(MaxNLocator(integer=True))

#ax.yaxis.set_major_locator(MaxNLocator(4))  #if the yaxis is smaller 1


plt.tight_layout()
#plt.savefig('/home/'+user+'/home/Polar_Low/Documents/ERA-PLs/Graphs/'+name+'.png', dpi= 70, bbox_inches='tight')

print(name, 'crit', crit)


#print('STARS PLs: ', nSPL, 'TRACK cyclones: ', nTPL)

if var=='small':
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