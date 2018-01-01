
# coding: utf-8

# T. Lavergne - 19 june 2017 02:02

# In[98]:

import os
from netCDF4 import Dataset, date2num
import pyresample as pr
import numpy as np
from datetime import datetime, timedelta,  time
import pyproj


# In[99]:


shapes = {'nh':(177,119), 'sh':(131,125)}
def get_areadef(hemis):
    nx, ny = shapes[hemis][::-1]
    a = 62500
    if hemis == 'nh':
        xl = -3781250.0
        xL = xl+0.5*a + (nx-1)*a + 0.5*a
        yL = +5781250.0
        yl = yL-0.5*a - (ny-1)*a - 0.5*a
        pdict = {u'a': u'6378273', u'b': u'6356889.44891', u'lat_ts': u'70', u'lon_0': u'-45', u'proj': u'stere', u'lat_0': u'90', u'units':u'm'}
    elif hemis == 'sh':
        xl = -3906250.0 
        xL = xl+0.5*a + (nx-1)*a + 0.5*a
        yL = 4281250.0
        yl = yL-0.5*a - (ny-1)*a - 0.5*a
        pdict = {u'a': u'6378273', u'b': u'6356889.44891', u'lat_ts': u'-70', u'lon_0': u'0', u'proj': u'stere', u'lat_0': u'-90', u'units':u'm'}
    else:
        raise ValueError("Not a valid area ('nh' and 'sh' are valid)")
    
    area_def  = pr.geometry.AreaDefinition('osi405_{}'.format(hemis,),
                                          'osi405_{}'.format(hemis,),
                                          'osi405_{}'.format(hemis,),
                                          pdict,
                                          nx, ny,
                                          [xl, yl, xL, yL])
    if area_def.pixel_size_x != a or area_def.pixel_size_y != a:
        raise ValueError("The grid parameters do not match a grid spacing of 62.5km! (x:{} y:{})".format(area_def.pixel_size_x, area_def.pixel_size_y))
    
    return area_def

def write_icedrift(odir,meth,hemis,dt0,dt1,dx,dy,):
    ### odir = output directory
    ### meth = ice drift method (e.g. 'winds', 'ida')
    ### hemis = hemisphere ('nh' or 'sh')
    ### dt0 = a datetime object for the start day of the drift at 12 o'clock
    ### dt1 = a datetime object for the stop day of the drift at 12 o'clock
    ### dx  = x component of the drift
    ### dy  = y component of the drift
    
    # create output name
    fname = 'ice_drift_{}_polstere-625_{}_{:%Y%m%d}1200-{:%Y%m%d}1200.nc'.format(hemis,meth,dt0,dt1)
    fname = os.path.join(odir,fname)
    
    adef = get_areadef(hemis)
    ny, nx = adef.shape
    lons, lats = adef.get_lonlats()
    
    grid_x, grid_y = np.meshgrid(adef.proj_x_coords,adef.proj_y_coords,)
    endpoint_x = grid_x + dx*1000.
    endpoint_y = grid_y + dy*1000.
    pj = pyproj.Proj(adef.proj4_string)
    endpoint_lons, endpoint_lats = pj(endpoint_x, endpoint_y,inverse=True)
    
    with Dataset(fname, 'w', format='NETCDF3_CLASSIC') as dataset:

        # dimension and auxiliary datasets
        dimx  = dataset.createDimension('xc', nx)
        dimy  = dataset.createDimension('yc', ny)
        dimt  = dataset.createDimension('time', 1)
        dimnv = dataset.createDimension('nv', 2)

        crs = dataset.createVariable('Polar_Stereographic_Grid', np.int32)
        pdict = dict([el.split('=') for el in adef.proj4_string.split()])
        crs.grid_mapping_name = 'polar_stereographic'
        crs.straight_vertical_longitude_from_pole = np.float32(pdict['+lon_0'])
        crs.latitude_of_projection_origin = np.float32(pdict['+lat_0'])
        crs.standard_parallel = np.float32(pdict['+lat_ts'])
        crs.false_easting = 0.0; 
        crs.false_northing = 0.0;
        crs.semi_major_axis = np.float32(pdict['+a'])
        crs.semi_minor_axis = np.float32(pdict['+b'])
        crs.proj4_string = str(adef.proj4_string)
        
        xc = dataset.createVariable('xc', np.float64, ('xc',))
        xc.axis = "X" 
        xc.units = "km" 
        xc.long_name = "x coordinate of projection (eastings)" 
        xc.standard_name = "projection_x_coordinate" 
        xc[:] = adef.proj_x_coords*0.001
        
        yc = dataset.createVariable('yc', np.float64, ('yc',))
        yc.axis = "Y" 
        yc.units = "km" 
        yc.long_name = "y coordinate of projection (northings)" 
        yc.standard_name = "projection_y_coordinate" 
        yc[:] = adef.proj_y_coords*0.001
        
        tunits = "seconds since 1978-01-01 00:00:00"
        time = dataset.createVariable('time', np.float64, ('time',))
        time.axis = "T" 
        time.long_name = "reference time of product" 
        time.standard_name = "time" 
        time.units =  tunits
        time.calendar = "standard" 
        time.bounds = "time_bnds"
        time[:] = date2num(dt1,tunits)
        
        tbounds = dataset.createVariable('time_bnds', np.float64, ('time','nv',))
        tbounds.units = tunits
        tbounds[:] = [date2num(dt0,tunits), date2num(dt1,tunits)]
        
        lat = dataset.createVariable('lat', np.float32, ('yc','xc'))
        lat.long_name = "latitude coordinate"
        lat.standard_name = "latitude" 
        lat.units = "degrees_north"
        lat[:] = lats
        
        lon = dataset.createVariable('lon', np.float32, ('yc','xc'))
        lon.long_name = "longitude coordinate"
        lon.standard_name = "longitude" 
        lon.units = "degrees_east"
        lon[:] = lons
        
        # data
        d0 = dataset.createVariable('dt0', np.int32, ('time','yc','xc'),)
        d0.long_name = "delta time for start of displacement" 
        d0.units = "seconds since {}".format(dt0)
        d0.valid_min = -43200 
        d0.valid_max = 43200 
        d0.grid_mapping = "Polar_Stereographic_Grid" 
        d0.coordinates = "lat lon"
        d0[0,:] = np.ma.array(np.zeros(dx.shape),mask=dx.mask)
        
        d1 = dataset.createVariable('dt1', np.int32, ('time','yc','xc'),)
        d1.long_name = "delta time for end of displacement" 
        d1.units = "seconds since {}".format(dt1)
        d1.valid_min = -43200 
        d1.valid_max = 43200 
        d1.grid_mapping = "Polar_Stereographic_Grid" 
        d1.coordinates = "lat lon"
        d1[0,:] = np.ma.array(np.zeros(dx.shape),mask=dx.mask)
        
        dX = dataset.createVariable('dX', np.float32, ('time','yc','xc'),)
        dX.long_name = "component of the displacement along the x axis of the grid"
        dX.standard_name = "sea_ice_x_displacement"
        dX.units = "km"
        dX.grid_mapping = "Polar_Stereographic_Grid"
        dX.coordinates = "lat lon"
        dX.ancillary_variables = "uncert_dX_and_dY status_flag"
        dX[0,:] = dx
        
        dY = dataset.createVariable('dY', np.float32, ('time','yc','xc'),)
        dY.long_name = "component of the displacement along the y axis of the grid"
        dY.standard_name = "sea_ice_y_displacement"
        dY.units = "km"
        dY.grid_mapping = "Polar_Stereographic_Grid"
        dY.coordinates = "lat lon"
        dY.ancillary_variables = "uncert_dX_and_dY status_flag"
        dY[0,:] = dy
        
        lat1 = dataset.createVariable('lat1', np.float32, ('time','yc','xc'),)
        lat1.long_name = "latitude at end of displacement"
        lat1.units = "degrees_north"
        lat1.grid_mapping = "Polar_Stereographic_Grid"
        lat1.coordinates = "lat lon"
        lat1[0,:] = endpoint_lats
        
        lon1 = dataset.createVariable('lon1', np.float32, ('time','yc','xc'),)
        lon1.long_name = "longitude at end of displacement"
        lon1.units = "degrees_east"
        lon1.grid_mapping = "Polar_Stereographic_Grid"
        lon1.coordinates = "lat lon"
        lon1[0,:] = endpoint_lons
        
    return fname
        


# In[101]:

if __name__ == '__main__':

    dt1 = datetime.today() - timedelta(days=3)
    dt1 = datetime.combine(dt1.date(),time(12))
    dt0 = dt1 - timedelta(days=2)

    meth = 'ida'
    # meth = 'wind'

    for hemis in ('nh','sh'):
        f="http://thredds.met.no/thredds/dodsC/osisaf/met.no/ice/drift_lr/merged/{d1:%Y/%m}/ice_drift_{a:}_polstere-625_multi-oi_{d0:%Y%m%d}1200-{d1:%Y%m%d}1200.nc".format(a=hemis,d0=dt0,d1=dt1)
        print (f)
        with Dataset(f) as inp:
            dx = inp.variables['dX'][0,:]
            dy = inp.variables['dY'][0,:]
        odir = '/tmp'
        oname = write_icedrift(odir,meth,hemis,dt0,dt1,dx,dy)
        print (oname)


# In[ ]:




# In[ ]:



