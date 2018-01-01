# -*- coding: utf-8 -*-
"""
Created on Sun Oct 29 23:08:07 2017

@author: admin
"""
#!/users/bin/env python

import os
from netCDF4 import Dataset, date2num
import pyresample as pr
import numpy as np
from datetime import datetime, timedelta,  time
import pyproj

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


def write_params(start_yr, end_yr, season, hemis, real_A, imag_A, real_C, imag_C, real_c, imag_c, absA, ThA, RMS_ResR, RMS_ResI, out_dir):
    """
    start_ yr = initial year of period 
    end_yr = final year of period
    season = 'summer' or 'winter'
    hemis = 'nh' or 'sh'
    real_A = real component of the complex A parameter
    imag_A = imaginary component of the complex A parameter
    real_C = real component of the complex C parameter
    imag_C = imaginary component of the complex C parameter
    real_c = normalized real component of the complex C parameter
    imag_c = normalized imaginary component of the complex C parameter
    absA = absolute value of complex A parameter
    ThA = angle value of complex A parameter
    RMS_ResR = yearly-averaged Root Mean Squared of real component of residuals
    RMS_ResI = yearly-averaged Root Mean Squared of imaginary component of residuals
    odir = output directory
    """
    
    fname = 'inversion_parameters_{}-{}_{}.nc'.format(start_yr, end_yr, season)
    fname = os.path.join(out_dir,fname)
    
    adef = get_areadef(hemis)
    ny, nx = adef.shape
    lons, lats = adef.get_lonlats()
    
#    grid_x, grid_y = np.meshgrid(adef.proj_x_coords,adef.proj_y_coords,)
#    endpoint_x = grid_x + dx*1000.
#    endpoint_y = grid_y + dy*1000.
#    pj = pyproj.Proj(adef.proj4_string)
#    endpoint_lons, endpoint_lats = pj(endpoint_x, endpoint_y,inverse=True)
    
    with Dataset(fname, 'w', format='NETCDF3_CLASSIC') as dataset:

        # dimensions 
        dimx  = dataset.createDimension('xc', nx)
        dimy  = dataset.createDimension('yc', ny)
#        dimt  = dataset.createDimension('time', 1)
#        dimnv = dataset.createDimension('nv', 2)

        # auxiliary (coordinate) variables
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
        crs.proj4_string = str(adef.proj4_string)       # ??? #
        
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
        A_real = dataset.createVariable('A_real', np.float32, ('yc','xc'),)
        A_real.long_name = "real component of the parameter A"
        A_real.standard_name = "A real"
        A_real.units = " "
        A_real.grid_mapping = "Polar_Stereographic_Grid"
        A_real.coordinates = "lat lon"
        A_real[:] = real_A
    
        A_imag = dataset.createVariable('A_imag', np.float32, ('yc','xc'),)
        A_imag.long_name = "imaginary component of the parameter A"
        A_imag.standard_name = "A imag"
        A_imag.units = " "
        A_imag.grid_mapping = "Polar_Stereographic_Grid"
        A_imag.coordinates = "lat lon"
        A_imag[:] = imag_A
        
        C_real = dataset.createVariable('C_real', np.float32, ('yc','xc'),)
        C_real.long_name = "real component of the parameter C"
        C_real.standard_name = "C real"
        C_real.units = " "
        C_real.grid_mapping = "Polar_Stereographic_Grid"
        C_real.coordinates = "lat lon"
        C_real[:] = real_C
    
        C_imag = dataset.createVariable('C_imag', np.float32, ('yc','xc'),)
        C_imag.long_name = "imaginary component of the parameter C"
        C_imag.standard_name = "C imag"
        C_imag.units = " "
        C_imag.grid_mapping = "Polar_Stereographic_Grid"
        C_imag.coordinates = "lat lon"
        C_imag[:] = imag_C
        
        absA = dataset.createVariable('absA', np.float32, ('yc','xc'),)
        absA.long_name = "absolute value of the parameter A"
        absA.standard_name = "abs A"
        absA.units = " "
        absA.grid_mapping = "Polar_Stereographic_Grid"
        absA.coordinates = "lat lon"
        absA[:] = absA
        
        ThetaA = dataset.createVariable('ThetaA', np.float32, ('yc','xc'),)
        ThetaA.long_name = "angle value of the parameter A"
        ThetaA.standard_name = "Theta A"
        ThetaA.units = "degrees"
        ThetaA.grid_mapping = "Polar_Stereographic_Grid"
        ThetaA.coordinates = "lat lon"
        ThetaA[:] = ThA

        c_real = dataset.createVariable('c_real', np.float32, ('yc','xc'),)
        c_real.long_name = "normalized real component of the parameter C"
        c_real.standard_name = "c real"
        c_real.units = " "
        c_real.grid_mapping = "Polar_Stereographic_Grid"
        c_real.coordinates = "lat lon"
        c_real[:] = real_c
    
        c_imag = dataset.createVariable('c_imag', np.float32, ('yc','xc'),)
        c_imag.long_name = "normalized imaginary component of the parameter C"
        c_imag.standard_name = "c imag"
        c_imag.units = " "
        c_imag.grid_mapping = "Polar_Stereographic_Grid"
        c_imag.coordinates = "lat lon"
        c_imag[:] = imag_c
        
        RMS_Res_real = dataset.createVariable('RMS_Res_real', np.float32, ('yc','xc'),)
        RMS_Res_real.long_name = "Root mean square of the real component of residuals"
        RMS_Res_real.standard_name = "RMS Res real"
        RMS_Res_real.units = " "
        RMS_Res_real.grid_mapping = "Polar_Stereographic_Grid"
        RMS_Res_real.coordinates = "lat lon"
        RMS_Res_real[:] = RMS_ResR
    
        RMS_Res_imag = dataset.createVariable('RMS_Res_imag', np.float32, ('yc','xc'),)
        RMS_Res_imag.long_name = "Root mean square of the imaginary component of residuals"
        RMS_Res_imag.standard_name = "RMS Res imag"
        RMS_Res_imag.units = " "
        RMS_Res_imag.grid_mapping = "Polar_Stereographic_Grid"
        RMS_Res_imag.coordinates = "lat lon"
        RMS_Res_imag[:] = RMS_ResI
        
    return fname
    
