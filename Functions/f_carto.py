# -*- coding: utf-8 -*-



import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.ticker as mticker
import numpy as np
import matplotlib.pyplot as plt

import cartopy.geodesic as cgeo





def Plot_Polar_Stereo(fig, central_longitude= 30, extent= [-10, 50, 60, 80],
                      scalebar=True, subplot= (1,1,1), hem='NH', circular= False):
    """plot a polar stereo map
    extent can be False for whole hem."""
    
    #ax = fig.add_subplot(1, 1, 1, projection= ccrs.Stereographic(central_longitude=10))
    sub_a, sub_b, sub_c= subplot
    if hem== 'NH': ax = fig.add_subplot(sub_a, sub_b, sub_c, projection= ccrs.NorthPolarStereo(central_longitude= central_longitude))
    elif hem== 'SH': ax = fig.add_subplot(sub_a, sub_b, sub_c, projection= ccrs.SouthPolarStereo(central_longitude= central_longitude))
    
    if extent: ax.set_extent(extent, crs=ccrs.PlateCarree())
    
    #ax.coastlines(resolution='50m', color='grey') 
    feature= cfeature.NaturalEarthFeature(name="coastline", category="physical", scale="50m", alpha=0.5, edgecolor="#333333", facecolor="#AAAAAA")
    ax.add_feature(feature)
    
    gl= ax.gridlines() #draw_labels=True)
    gl.xlocator = mticker.FixedLocator(np.arange(-180, 180, 30))
    gl.ylocator = mticker.FixedLocator(np.arange(-80, 80.1, 10))
    if scalebar: scale_bar(ax, 500, location= (0.06, 0.04))
    

    if circular:
        theta = np.linspace(0, 2*np.pi, 100)
        center, radius = [0.5, 0.5], 0.5
        verts = np.vstack([np.sin(theta), np.cos(theta)]).T
        
        import matplotlib.path as mpath
        circle = mpath.Path(verts * radius + center)

        ax.set_boundary(circle, transform=ax.transAxes)
    
    return ax


def Plot_PlateCarree(fig, central_longitude= 0, extent= [-15, 55, 50, 80],
                     scalebar=True, subplot= (1,1,1)):
    sub_a, sub_b, sub_c= subplot
    ax = fig.add_subplot(sub_a, sub_b, sub_c, projection= ccrs.PlateCarree(central_longitude= central_longitude))
    
    ax.set_extent(extent, crs=ccrs.PlateCarree())
    
#    ax.coastlines(resolution='50m', color='grey') 
    feature= cfeature.NaturalEarthFeature(name="coastline", category="physical", scale="50m", alpha=0.5, edgecolor="#333333", facecolor="#AAAAAA")
    ax.add_feature(feature)
    
    gl= ax.gridlines() #draw_labels=True)
    gl.xlocator = mticker.FixedLocator(np.arange(-180, 180, 10))
    gl.ylocator = mticker.FixedLocator(np.arange(-80, 80.1, 5))
    if scalebar: scale_bar(ax, 500, location= (0.06, 0.04))


    return ax






def PlotWind(lon, lat, u, v, ax, nx= 15, ny=20, arrowtype= 'quiver', alen=None, scale=False, color= 'k'):
    """ winds with barbs
    alen= arrow length: the lenght of the arrows can be specified by alen
    arrowtype= ('quiver', 'barb') - chose one type of arrows, quiver are normal arrows, barbs have a tail
    nx, ny, specifies how many arrows there are in the x and y direction
    """
#    print('plot winds')
    if len(np.shape(lon))==1: #for ERA data lon and lat are put into 2D array
        lon, lat= np.meshgrid(lon, lat)
        
    everyx= np.shape(u)[0]//nx #put a barb in every nth box, such that you have 20 barbs in the longer dim
    everyy= np.shape(u)[1]//ny #20 for a quadratical map , 50 for era
    u = u[0::everyx, 0::everyy]
    v = v[0::everyx, 0::everyy]
    lon= lon[0::everyx, 0::everyy]
    lat= lat[0::everyx, 0::everyy]
    
    
#    if arrowtype== 'barb':
#        map.barbs(Lon,Lat,u,v, np.sqrt(u**2+v**2), cmap= plt.cm.Reds)
        
    if arrowtype== 'quiver':
        if alen != None:
#            Q = ax.quiver(lon,lat,u,v, pivot= 'mid', scale=20* alen, headwidth= 4, width=0.005, transform=ccrs.PlateCarree(), color= color)           
            #this should fix for wrong cartopy arrow direction:
            Q = ax.quiver(lon,lat, u/np.cos(lat /180 * np.pi),v, pivot= 'mid', scale=20* alen,
                          headwidth= 4, width=0.005, transform=ccrs.PlateCarree(), color= color , angles = "xy")
            qk = plt.quiverkey(Q, 0.18, 0.98, alen, str(alen)+' m/s', coordinates='figure', labelpos='W')
        
#        else:        
#            Umax= np.max(np.sqrt(u**2+v**2)) #//5*5 #rounds to nearest 5
#    #        print(Umax)
#            if Umax < 40:
#
#                if scale== False: scale= 500
#                Q = map.quiver(Lon,Lat,u,v, pivot= 'mid', scale=scale, headwidth= 4, width=0.005)
#                qk = plt.quiverkey(Q, 0.18, 0.98, 20, '20 m/s', coordinates='figure', labelpos='W')
#            else:
#                if scale== False: scale= 1000
#                Q = map.quiver(Lon,Lat,u,v, pivot= 'mid', scale=scale, headwidth= 4, width=0.005)
#                qk = plt.quiverkey(Q, 0.18, 0.98, 75, '75 m/s', coordinates='figure' , labelpos='W')        
#        #    Q = map.quiver(Lon,Lat,uproj,vproj, pivot= 'mid', scale=25*Umax)
#        #    qk = plt.quiverkey(Q, 0.22, 0.97, Umax, str(int(Umax))+' m/s', labelpos='W')





def scale_bar(ax, length, location= (0.06, 0.04), metres_per_unit=1000, unit_name='km',
              tol=0.01, angle=0, color='black', linewidth=3, text_offset=0.005,
              ha='center', va='bottom', plot_kwargs=None, text_kwargs=None,
              **kwargs):
    """Add a scale bar to CartoPy axes.

    For angles between 0 and 90 the text and line may be plotted at
    slightly different angles for unknown reasons. To work around this,
    override the 'rotation' keyword argument with text_kwargs.

    Args:
        ax:              CartoPy axes.
        location:        Position of left-side of bar in axes coordinates.
        length:          Geodesic length of the scale bar.
        metres_per_unit: Number of metres in the given unit. Default: 1000
        unit_name:       Name of the given unit. Default: 'km'
        tol:             Allowed relative error in length of bar. Default: 0.01
        angle:           Anti-clockwise rotation of the bar.
        color:           Color of the bar and text. Default: 'black'
        linewidth:       Same argument as for plot.
        text_offset:     Perpendicular offset for text in axes coordinates.
                         Default: 0.005
        ha:              Horizontal alignment. Default: 'center'
        va:              Vertical alignment. Default: 'bottom'
        **plot_kwargs:   Keyword arguments for plot, overridden by **kwargs.
        **text_kwargs:   Keyword arguments for text, overridden by **kwargs.
        **kwargs:        Keyword arguments for both plot and text.
    """
    # Setup kwargs, update plot_kwargs and text_kwargs.
    if plot_kwargs is None:
        plot_kwargs = {}
    if text_kwargs is None:
        text_kwargs = {}

    plot_kwargs = {'linewidth': linewidth, 'color': color, **plot_kwargs,
                   **kwargs}
    text_kwargs = {'ha': ha, 'va': va, 'rotation': angle, 'color': color,
                   **text_kwargs, **kwargs}

    # Convert all units and types.
    location = np.asarray(location)  # For vector addition.
    length_metres = length * metres_per_unit
    angle_rad = angle * np.pi / 180

    # End-point of bar.
    end = _point_along_line(ax, location, length_metres, angle=angle_rad,
                            tol=tol)

    # Coordinates are currently in axes coordinates, so use transAxes to
    # put into data coordinates. *zip(a, b) produces a list of x-coords,
    # then a list of y-coords.
    ax.plot(*zip(location, end), transform=ax.transAxes, **plot_kwargs)

    # Push text away from bar in the perpendicular direction.
    midpoint = (location + end) / 2
    offset = text_offset * np.array([-np.sin(angle_rad), np.cos(angle_rad)])
    text_location = midpoint + offset

    # 'rotation' keyword argument is in text_kwargs.
    ax.text(*text_location, f"{length} {unit_name}", rotation_mode='anchor',
            transform=ax.transAxes, **text_kwargs)





def plot_circle(ax, lon, lat, radius, edgecolor= 'k'):
    """plot a circle around the center point (lon, lat) with radius [km]"""
    
    def compute_radius(ortho, radius_degrees, lon, lat):
        phi1 = lat + radius_degrees if lat <= 0 else lat - radius_degrees
        _, y1 = ortho.transform_point(lon, phi1, ccrs.PlateCarree())
        return abs(y1)

    #the projection is important for equal distance
    #proj= ccrs.NorthPolarStereo(central_longitude= 30)
    proj= ccrs.Orthographic(central_longitude=lon, central_latitude=lat)
    
    radius_degrees= radius/110
    r_ortho= compute_radius(proj, radius_degrees, lon, lat)
    
    import matplotlib.patches as mpatches
    ax.add_patch(mpatches.Circle(xy=[lon, lat], radius=r_ortho, edgecolor=edgecolor, fill= 0,
                                 transform= proj))
    




"""---------------- some of Mathias functions ------------------ """





def _axes_to_lonlat(ax, coords):
    """(lon, lat) from axes coordinates."""
    display = ax.transAxes.transform(coords)
    data = ax.transData.inverted().transform(display)
    lonlat = ccrs.PlateCarree().transform_point(*data, ax.projection)

    return lonlat


def _upper_bound(start, direction, distance, dist_func):
    """A point farther than distance from start, in the given direction.

    It doesn't matter which coordinate system start is given in, as long
    as dist_func takes points in that coordinate system.

    Args:
        start:     Starting point for the line.
        direction  Nonzero (2, 1)-shaped array, a direction vector.
        distance:  Positive distance to go past.
        dist_func: A two-argument function which returns distance.

    Returns:
        Coordinates of a point (a (2, 1)-shaped NumPy array).
    """
    if distance <= 0:
        raise ValueError(f"Minimum distance is not positive: {distance}")

    if np.linalg.norm(direction) == 0:
        raise ValueError("Direction vector must not be zero.")

    # Exponential search until the distance between start and end is
    # greater than the given limit.
    length = 0.1
    end = start + length * direction

    while dist_func(start, end) < distance:
        length *= 2
        end = start + length * direction

    return end


def _distance_along_line(start, end, distance, dist_func, tol):
    """Point at a distance from start on the segment  from start to end.

    It doesn't matter which coordinate system start is given in, as long
    as dist_func takes points in that coordinate system.

    Args:
        start:     Starting point for the line.
        end:       Outer bound on point's location.
        distance:  Positive distance to travel.
        dist_func: Two-argument function which returns distance.
        tol:       Relative error in distance to allow.

    Returns:
        Coordinates of a point (a (2, 1)-shaped NumPy array).
    """
    initial_distance = dist_func(start, end)
    if initial_distance < distance:
        raise ValueError(f"End is closer to start ({initial_distance}) than "
                         f"given distance ({distance}).")

    if tol <= 0:
        raise ValueError(f"Tolerance is not positive: {tol}")

    # Binary search for a point at the given distance.
    left = start
    right = end

    while not np.isclose(dist_func(start, right), distance, rtol=tol):
        midpoint = (left + right) / 2

        # If midpoint is too close, search in second half.
        if dist_func(start, midpoint) < distance:
            left = midpoint
        # Otherwise the midpoint is too far, so search in first half.
        else:
            right = midpoint

    return right


def _point_along_line(ax, start, distance, angle=0, tol=0.01):
    """Point at a given distance from start at a given angle.

    Args:
        ax:       CartoPy axes.
        start:    Starting point for the line in axes coordinates.
        distance: Positive physical distance to travel.
        angle:    Anti-clockwise angle for the bar, in radians. Default: 0
        tol:      Relative error in distance to allow. Default: 0.01

    Returns:
        Coordinates of a point (a (2, 1)-shaped NumPy array).
    """
    # Direction vector of the line in axes coordinates.
    direction = np.array([np.cos(angle), np.sin(angle)])

    geodesic = cgeo.Geodesic()

    # Physical distance between points.
    def dist_func(a_axes, b_axes):
        a_phys = _axes_to_lonlat(ax, a_axes)
        b_phys = _axes_to_lonlat(ax, b_axes)

        # Geodesic().inverse returns a NumPy MemoryView like [[distance,
        # start azimuth, end azimuth]].
        return geodesic.inverse(a_phys, b_phys).base[0, 0]

    end = _upper_bound(start, direction, distance, dist_func)

    return _distance_along_line(start, end, distance, dist_func, tol)


