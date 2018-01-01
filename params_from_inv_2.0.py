# -*- coding: utf-8 -*-
"""
Created on Sat Nov 25 01:58:17 2017

@author: javimozo
"""
#!/users/bin/env python

import os
import sys
import numpy as np
import numpy.ma as ma
from numpy.linalg import inv
from netCDF4 import Dataset
from write_params_Athc import write_params
from datetime import datetime, timedelta, time, date


def load_files(i, season, wind_dir, ice_dir):

    if season == 'summer':
        start_date = date(i,6,1)
        end_date = date(i,9,30)
        wind_files = []
        ice_files = []

        st_d = start_date                
        while st_d <= end_date:
            d0 = st_d.strftime("%Y%m%d")
            st_d += timedelta(days=2)
            d1 = st_d.strftime("%Y%m%d")
            
            windf = 'NWP_nhOSISAF_aggr_{}12-{}12.nc'.format(d0, d1)
            windf = os.path.join(wind_dir,windf)
            wind_files.append(windf)
            icef = 'ice_drift_nh_polstere-625_multi-oi_{}1200-{}1200.nc'.format(d0, d1)
            icef = os.path.join(ice_dir,icef)
            ice_files.append(icef)
            
        wind_files_size = np.size(wind_files)
        ice_files_size = np.size(ice_files)
        
    elif season == 'winter':
        start_date = date(i-1,10,1)
        end_date = date(i,5,31)
        wind_files = []
        ice_files = []
        
        st_d = start_date
        while st_d <= end_date:
            d0 = st_d.strftime("%Y%m%d")
            st_d += timedelta(days=2)
            d1 = st_d.strftime("%Y%m%d")
            
            windf = 'NWP_nhOSISAF_aggr_{}12-{}12.nc'.format(d0, d1)
            windf = os.path.join(wind_dir,windf)
            wind_files.append(windf)
            icef = 'ice_drift_nh_polstere-625_multi-oi_{}1200-{}1200.nc'.format(d0, d1)
            icef = os.path.join(ice_dir,icef)
            ice_files.append(icef)
            
        wind_files_size = np.size(wind_files)
        ice_files_size = np.size(ice_files)
        #print(ice_files)
    else:
        raise 'Error'	# elaborate !!!

    return (wind_files, ice_files, wind_files_size, ice_files_size)


def load_data(wind_files, ice_files, wind_files_size, ice_files_size):      
   
    Wu = np.zeros([177*119,wind_files_size],dtype=np.float)    # 2D array to save LOCAL (X) WIND VELOCITY values 
    Wv = np.zeros([177*119,wind_files_size],dtype=np.float)    # 2D array to save MERIDIONAL (Y) WIND VELOCITY values 
    Iu = ma.zeros([177*119,ice_files_size],dtype=np.float)    # 2D array to save LOCAL (X) ICE VELOCITY values 
    Iv = ma.zeros([177*119,ice_files_size],dtype=np.float)    # 2D array to save MERIDIONAL (X) ICE VELOCITY values

    #Wu[:,:] = np.nan
    #Wv[:,:] = np.nan
    #Iu[:,:] = np.nan
    #Iv[:,:] = np.nan

    n = -1       # initial value of counter  
    nn = -1

    for windf in wind_files:
        with Dataset(windf, mode='r') as datum1:
	    lat1 = datum1.variables['lat1'][:]
  	    lon1 = datum1.variables['lon1'][:]
	    wu = datum1.variables['u_wind_avg'][:]
	    wv = datum1.variables['v_wind_avg'][:]
    
        gsiz = lat1.size
        gshp = lat1.shape
    
        n += 1
        Wu[:,n] = wu.flat           
        Wv[:,n] = wv.flat      
    #print(Wu)
    for icef in ice_files:
        with Dataset(icef, mode='r') as datum2:
            lat1 = datum2.variables['lat1'][:]
            lon1 = datum2.variables['lon1'][:]
            iu = datum2.variables['dX'][:]
            try:
                iv = datum2.variables['dY_v1p4'][:]
            except KeyError:
                iv = datum2.variables['dY'][:]
    
        iu = iu[0,:,:]*(1000.0/(2*86400))       # convert to m/s from km/2days
        iv = iv[0,:,:]*(1000.0/(2*86400))       # before transforming to columns
        
        nn += 1
        Iu[:,n] = ma.ravel(iu)             # saves LOCAL (X) ICE VELOCITY values in columns 
        Iv[:,n] = ma.ravel(iv)             # saves MERIDIONAL (Y) ICE VELOCITY values in columns 
    #print(Iu)  
    return (Wu, Wv, Iu, Iv)


def inversion(Wu, Wv, Iu, Iv, yr_range, season, m):

    C = np.zeros([177*119,len(yr_range)],dtype=np.complex)     # 2D array to save A and C values 
    A = np.zeros([177*119,len(yr_range)],dtype=np.complex)     # as they're computed for each grid point (rows) for every year (colunmns)
    C[:,:].real = np.nan                                 	# changes zeros to nans
    C[:,:].imag = np.nan
    A[:,:].real = np.nan
    A[:,:].imag = np.nan

    RMSresR = np.zeros([177*119,len(yr_range)],dtype=np.float)    # 2D array to save LOCAL (X) residual values 
    RMSresI = np.zeros([177*119,len(yr_range)],dtype=np.float)    # 2D array to save MERIDIONAL (Y) residual values 
    RMSresR[:,:] = np.nan
    RMSresI[:,:] = np.nan
    
    for element in range(177*119):         # loop over every gridpoint
        
        #print ('element: ', element)
        
        WU = Wu[element,:]              # extracts one by one every row (which represents wind values for each grid point)
        WV = Wv[element,:]        
        IU = Iu[element,:]        
        IV = Iv[element,:]              # extracts one by one every row (which represents ice drift values for each grid point)

        k = ma.getmask(IU)                # checks for nans

        WC = np.array(WU,dtype=np.complex)    
        WC.imag = (WV)
        IC = np.array(IU,dtype=np.complex)
        IC.imag = (IV)

        WC = WC[~k]                     # removes nans
        IC = IC[~k]

	if season == 'summer':        		# if summer --> len(IC) > 5 _ minimum lengh of time series to compute
        	if len(IC) > 5:                
            		ID = np.eye(len(WC))
            		IC = IC[:,np.newaxis]                             # transforms the horizontal array into a vertical one
            		WC = [WC,np.ones_like(WC,dtype=np.complex)]       # attaches an array of ones to the list of complex velocity values
            		WC = np.mat(WC)                                   # converts a list of a 1D list and 1D array into a 2 row matrix
            		WC = WC.T                                         # transposes the matrix into a 2 column one (1 x 2)
            		WCH = WC.getH()                                   # transforms it into a Hermitian (conjugate transpose) matrix
                 
            		AC = (inv(WCH*WC)*WCH*IC)         # computes the inversion // 1x2 matrix w/ values of A and C (in complex form)

            		A[element,m] = AC[0,0]       # saves A values by columns
            		C[element,m] = AC[1,0]       # appends C values by columns

            		Res = (ID - (WC*(inv(WCH*WC))*WCH))*IC

            		RMSresR[element,m] = np.sqrt(np.mean(np.square(Res.real)))         
            		RMSresI[element,m] = np.sqrt(np.mean(np.square(Res.imag)))         

	elif season == 'winter':       		 # if winter --> len(IC) > 10 _ minimum lengh of time series to compute
        	if len(IC) > 10:                
            		ID = np.eye(len(WC))
            		IC = IC[:,np.newaxis]                             # transforms the horizontal array into a vertical one
            		WC = [WC,np.ones_like(WC,dtype=np.complex)]       # attaches an array of ones to the list of complex velocity values
            		WC = np.mat(WC)                                   # converts a list of a 1D list and 1D array into a 2 row matrix
            		WC = WC.T                                         # transposes the matrix into a 2 column one (1 x 2)
            		WCH = WC.getH()                                   # transforms it into a Hermitian (conjugate transpose) matrix
                 
            		AC = (inv(WCH*WC)*WCH*IC)         # computes the inversion // 1x2 matrix w/ values of A and C (in complex form)

            		A[element,m] = AC[0,0]       # saves A values by columns
            		C[element,m] = AC[1,0]       # appends C values by columns

            		Res = (ID - (WC*(inv(WCH*WC))*WCH))*IC

            		RMSresR[element,m] = np.sqrt(np.mean(np.square(Res.real)))         
            		RMSresI[element,m] = np.sqrt(np.mean(np.square(Res.imag)))         
    #print(A)
    return (A, C, RMSresR, RMSresI)


def average(A, C, RMSresR, RMSresI):
    
    CC = np.zeros([177*119,],dtype=np.complex)     # 2D array to save A and C values 
    AA = np.zeros([177*119,],dtype=np.complex)     # as they're computed for each grid point (rows) for every year (colunmns)
    CC[:,].real = np.nan                                 # changes zeros to nans
    CC[:,].imag = np.nan
    AA[:,].real = np.nan
    AA[:,].imag = np.nan

    RMS_ResR = np.zeros([177*119,],dtype=np.complex)    # 2D array to save LOCAL (X) RMS residual values
    RMS_ResI = np.zeros([177*119,],dtype=np.complex)    # 2D array to save MERIDIONAL (Y) RMS residual values 
    RMS_ResR[:,] = np.nan
    RMS_ResI[:,] = np.nan

    for ii in range(177*119):                  # loops over every grid point to compute the mean over the whole time series
        CC[ii] = np.nanmean(C[ii,:])
        AA[ii] = np.nanmean(A[ii,:])

    for ii in range(177*119):
        RMS_ResR[ii] = np.nanmean(RMSresR[ii,:])         # saves LOCAL (X) ICE VELOCITY values in columns           
        RMS_ResI[ii] = np.nanmean(RMSresI[ii,:])         # saves LOCAL (X) ICE VELOCITY values in columns           

    absA = np.absolute(AA)                           # |A| informs of coupling btw wind/ice (and internal ice stresses)
    ThA = np.angle(AA,deg=True)              # angle btw wind/ice motion vectors
       
    absA = np.reshape(absA,(177,119))             # reshapes the A parameter 1D array into 2D for plotting with contourf
    ThA = np.reshape(ThA,(177,119))              # angle btw wind/ice motion vectors
    
    modulus = abs(CC)
    
    c = CC/modulus                       # normalized vectors
    c = np.reshape(c,(177,119))
    
    AA = np.reshape(A,(177,119))
    CC = np.reshape(C,(177,119))
    
    RMS_ResR = np.reshape(RMS_ResR,(177,119))
    RMS_ResI = np.reshape(RMS_ResI,(177,119))
    
    real_C = CC.real
    imag_C = CC.imag
    real_A = AA.real
    imag_A = AA.imag
    real_c = c.real
    imag_c = c.imag
    
    return (real_A, imag_A, real_C, imag_C, real_c, imag_c, absA, ThA, RMS_ResR, RMS_ResI)

    
def compute_params(start_yr, stop_yr, season, wind_dir, ice_dir, out_dir='.'):
 
    yr_range = range(int(start_yr), int(stop_yr))

    m = -1
    for i in yr_range:
        m += 1
        wind_files, ice_files, wind_files_size, ice_files_size = load_files(i, season, wind_dir, ice_dir)
        Wu, Wv, Iu, Iv = load_data(wind_files, ice_files, wind_files_size, ice_files_size)
        A, C, RMSresR, RMSresI = inversion(Wu, Wv, Iu, Iv, yr_range, season, m)
    real_A, imag_A, real_C, imag_C, real_c, imag_c, absA, ThA, RMS_ResR, RMS_ResI = average(A, C, RMSresR, RMSresI)

    hemis = 'nh'
    dfile = write_params(start_yr, stop_yr, season, hemis, real_A, imag_A, real_C, imag_C, real_c, imag_c, absA, ThA, RMS_ResR, RMS_ResI, out_dir)
    print ('Done with {}'.format(dfile))
    
    
if __name__ == '__main__':
    
    import argparse
    from argparse import RawDescriptionHelpFormatter
    
    p = argparse.ArgumentParser('params_from_inv', formatter_class=RawDescriptionHelpFormatter)
    p.add_argument('start_yr', help='starting date (integer of format YYYY) for computing parameters')
    p.add_argument('stop_yr', help='ending date (integer of format YYYY) for computing parameters')
    p.add_argument('season', help='string ("summer" or "winter")')
    p.add_argument('--wind_dir', help='where to look for wind files')
    p.add_argument('--ice_dir', help='where to look for ice files')
    p.add_argument('--out_dir', help='where to write output file')
    args = p.parse_args()
        
    compute_params(args.start_yr, args.stop_yr, args.season, args.wind_dir, args.ice_dir, out_dir=args.out_dir)    
    
    sys.exit('Done with params_from_inv')
