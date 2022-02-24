#!/usr/bin/env python
"""
Regrid data in a regular lat-lon grid to an equal area EASE grid.
Plot an example of Arctic SLP (sea-level pressure).
using files such as:
    slp_erai_0.75deg_6h_20021201-20180131.nc
    where the variable is SLP(time,lat,lon)
Output: slp_erai_0.75deg_6h_20021201-20180131_eag.nc
Change input variables: dates and netcdf variable name in the end of this script.
"""

__author__  = "<petteri.uotila@helsinki.fi>"

import sys
import os
import re
import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from mpl_toolkits.basemap import Basemap, interp
from datetime import datetime, timedelta

dcmap = LinearSegmentedColormap.from_list(name='custom_div_cmap',colors=['blue','white','red'],N=11)

class NCVar(object):
    def __init__(self,sdate='19980101',edate='20141231',freq='6h',\
                 ncvar='MSL',resolution='0.75deg',dset='erai',\
                 latrange=[55,90],FillValue=0.):
        self.dset, self.sdate, self.edate = dset, sdate, edate
        self.freq, self.resolution = freq, resolution
        self.bn = "%s_%s_%s_%s_%s-%s" % \
                   (ncvar.lower(),dset,resolution,freq,sdate,edate)
        self.ncvar = ncvar
        self.fn = "%s.nc" % self.bn
        self.latrange = latrange
        self.FillValue = FillValue

    def formDateArray(self,fp):
        time = fp.variables['time']
        self.cdftime = nc.netcdftime.utime(time.units,calendar=time.calendar)
        self.timeunits = time.units
        self.calendar  = time.calendar
        return [self.cdftime.num2date(t) for t in time[:]]

    def readData(self):
        try:
            fp = nc.Dataset(self.fn)
        except:
            print("Cant read %s!" % self.fn)
            sys.exit(0)
        self.dates = self.formDateArray(fp)
        lat = np.array(fp.variables['lat'][:])
        latidx = np.where((lat>=self.latrange[0])&(lat<self.latrange[1]))
        #latidx = np.where(lat>=self.latlim)
        self.lat = lat[latidx]
        self.lon = np.array(fp.variables['lon'][:])
        self.d = np.ma.array(fp.variables[self.ncvar][:,latidx[0],:]).squeeze()
        fp.close()

    def getMapParams(self):
        parallels = np.arange(-90.,90.,10.)
        meridians = [0,45,90,135,180,225,270,315]
        m = Basemap(width=3100000,height=3100000, rsphere=(6378137.00,6356752.3142), resolution='l', area_thresh=1000.,projection='lcc',lat_0=67.,lon_0=17.)
        x,y = m(*np.meshgrid(self.lon,self.lat))
        return m, x, y, parallels, meridians

    def vectorTransform2map(self,u,v,m=None):
        """ Rotate and scale lat, lon directed vectors to
            the map projection
        """
        if m is None:
            m, x, y, parallels, meridians = self.getMapParams()
        lon2d, lat2d = np.meshgrid(self.lon, self.lat)
        x1, y1 = m(lon2d+u, lat2d+v)
        u_map, v_map = x1-x, y1-y
        # Rescale the magnitudes of the vectors...
        mag_scale = np.ma.hypot(u_map, v_map) / np.ma.hypot(u, v)
        i = np.ma.where(mag_scale>0)
        u_map[i] /= mag_scale[i]
        v_map[i] /= mag_scale[i]
        return u_map, v_map

    def interpolate2eag(self,fldin,m):
        """ Interpolate field to equal area in a npstereo
            from a regular lat-lon
            Note: basemap's interp seems only work when
            lon ranges from -180 to 180
        """
        # curvilinear (EASE type?) grid for plotting
        r = np.arange(0.1,3.9,.05)*1e6
        xout,yout = np.meshgrid(r,r)
        lonout,latout = m(xout,yout,inverse=True)
        if self.lon.max()>180:
            lidx = np.where(self.lon>180)[0][0]
            lonin = np.hstack((self.lon[lidx:]-360,self.lon[:lidx]))
            fldin = np.hstack((fldin[:,lidx:],fldin[:,:lidx]))
        else:
            lonin=self.lon
        if self.lat[0]>self.lat[1]:
            fldout = interp(fldin[::-1,:],lonin,self.lat[::-1],lonout,latout,masked=0)
        else:
            fldout = interp(fldin,lonin,self.lat,lonout,latout,masked=0)
        #fldout[np.where(latout<self.latlim)] = self.FillValue
        print(xout)
        return fldout, xout, yout

    def outputNetCDF(self):
        m, x, y, parallels, meridians = self.getMapParams()
        for didx, date in enumerate(self.dates):
            field, xeag, yeag = self.interpolate2eag(self.d[didx],m)
            if didx==0:
                fp = nc.Dataset(self.bn+'_eag_sca.nc','w',format='NETCDF3_CLASSIC')
                lon, lat = m(xeag,yeag,inverse=True)
        #        fp = nc.Dataset(self.bn+'_eag_sca.nc','w') 
                fp.createDimension('x',field.shape[0])
                fp.createDimension('y',field.shape[1])
                #fp.createDimension('xy',field.shape[0]*field.shape[1])
                fp.createDimension('time',None)
                nav_lon = fp.createVariable('nav_lon','f',('y','x',))
                nav_lon.setncattr('name','longitude')
                nav_lon.setncattr('units','degrees_east')
                nav_lon.setncattr('axis','X')
                nav_lon[:,:] = lon
                nav_lat = fp.createVariable('nav_lat','f',('y','x',))
                nav_lat.setncattr('name','latitude')
                nav_lat.setncattr('units','degrees_north')
                nav_lat.setncattr('axis','Y')
                nav_lat[:,:] = lat
                time = fp.createVariable('time','f',('time',))
                time.setncattr('name','time')
                time.setncattr('calendar',self.calendar)
                time.setncattr('axis','T')
                time.setncattr('units',self.timeunits)
                vari = fp.createVariable(self.ncvar,'f',('time','y','x',),fill_value=self.FillValue)
                vari.setncattr('coordinates','time nav_lat nav_lon')
                vari.setncattr('standard_name','air_pressure_at_sea_level')
                vari.setncattr('long_name','Mean sea-level pressure')
                vari.setncattr('units','Pa')
                vari.setncattr('missing_value',self.FillValue)
                #vari1d = fp.createVariable(self.ncvar+'1d','f',('time','xy',),fill_value=self.FillValue)
                #vari1d.setncattr('standard_name','air_pressure_at_sea_level')
                #vari1d.setncattr('long_name','Mean sea-level pressure')
                #vari1d.setncattr('units','Pa')
                #vari1d.setncattr('missing_value',self.FillValue)
            time[didx] = self.cdftime.date2num(date)
            vari[didx,:,:] = field
            #vari1d[didx,:] = field.reshape((field.shape[0]*field.shape[1],))
        fp.setncattr('history',"Created by %s on %s. For help contact %s." % \
                     (os.path.basename(__file__),\
                      datetime.today().strftime("%Y-%m-%d"),\
                      __author__))
        fp.close()

    def plotMap(self):
        m, x, y, parallels, meridians = self.getMapParams()
        fig = plt.figure()
        ax = fig.add_subplot(1,1,1)
        field, xeag, yeag = self.interpolate2eag(self.d.mean(axis=0),m)
        cs = m.pcolormesh(xeag,yeag,field)
        cb = fig.colorbar(cs)
        m.drawcoastlines()
        m.drawparallels(parallels,labels=[0,0,0,0],fontsize=8,linewidth=0.5)
        m.drawmeridians(meridians,labels=[1,0,0,1],fontsize=8,linewidth=0.5)
        for suf in ['.png']:
            plt.savefig(self.bn+suf)
        plt.close()

if __name__=="__main__":
    slp = NCVar(sdate='19980101',\
                edate='20141231',\
                ncvar='MSL')
    #slp = NCVar(sdate='19790101',\
    #            edate='20180131',\
    #            freq='1d',\
    #            ncvar='T2M')
    # read data
    slp.readData()
    # plot data (optional)
    # slp.plotMap()
    # output data
    slp.outputNetCDF()
    print("Finnished!")
