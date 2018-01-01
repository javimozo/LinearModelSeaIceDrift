#!/usr/bin/env python
"""
Created on Thu Jun 14 16:48:36 2017

@author: javier
"""

import os
import sys
import numpy as np
from netCDF4 import Dataset
from WriteLRSID import write_icedrift
from datetime import datetime, timedelta, time

def compute_drift(dt, wind_dir, etcd='./etc/', out_dir='.'):

    dt0 = datetime.combine(dt,time(12))
    dt1 = dt0 + timedelta(days=2)
    windf = 'NWP_nhOSISAF_aggr_{:%Y%m%d%H}-{:%Y%m%d%H}.nc'.format(dt0,dt1)
    windf = os.path.join(wind_dir,windf)
    
    # load information from wind file
    with Dataset(windf, mode='r') as datum1:
        lat1 = datum1.variables['lat1'][:]
        lon1 = datum1.variables['lon1'][:]
        wu = datum1.variables['u_wind_avg'][:]
        wv = datum1.variables['v_wind_avg'][:]
        sic = datum1.variables['sic_avg'][:]

    gsiz = lat1.size
    gshp = lat1.shape

    # applying the model requires complex numbers
    wc = np.array(wu,dtype=np.complex)
    wc.imag = wv

    Wc = np.zeros([gsiz,],dtype=np.complex)
    Wc[:,].real = np.nan                        # changes zeros to nans
    Wc[:,].imag = np.nan

    Wc[:,] = wc.flat

    # load the free-drift parameters
    with Dataset(os.path.join(etcd,'AcTh__2014-2015.nc'), mode='r') as datum2:
        a_real = datum2.variables['A_real'][:]
        a_imag = datum2.variables['A_imag'][:]
        c_real = datum2.variables['C_real'][:]
        c_imag = datum2.variables['C_imag'][:]

    a = np.array(a_real,dtype=np.complex)           # recombine real and imaginary parts of
    a.imag = a_imag                                 # both model parameters
    c = np.array(c_real,dtype=np.complex)
    c.imag = c_imag

    ac = np.zeros((gsiz,2),dtype=np.complex)       # 2D array to save A nad C values 
    ac[:,:].real = np.nan                           # changes zeros to nans
    ac[:,:].imag = np.nan
    ac[:,0] = a.flat
    ac[:,1] = c.flat

    II = np.zeros([gsiz,],dtype=np.complex)    # 1D array to save iceDrift values 
    II[:,].real = np.nan                        # changes zeros to nans
    II[:,].imag = np.nan

    # loop through the grid cells and apply the model
    for i in range(gsiz,):         # loop over every gridpoint
        WC = Wc[i,]              # extracts one by one every row (which represents wind values for each grid point)

        WC = [WC,np.ones_like(WC,dtype=np.complex)]       # attaches an array of ones to the list of complex velocity values
        WC = np.mat(WC)

        AC = ac[i,:]
        AC = np.mat(AC)
        AC = AC.T

        I = WC*AC
        II[i,] = I[0,0]#

    modulus = abs(II)
    ic = II/modulus          # normalized vectors
    ic = np.reshape(ic,gshp)

    # get to the dX dY components (TODO: the 2*24 is hard coded)
    tspan_s = (dt1-dt0).total_seconds()
    dx = np.ma.masked_invalid(ic.real*(tspan_s / 1000.0))
    dy = np.ma.masked_invalid(ic.imag*(tspan_s / 1000.0))

    # write netCDF 
    meth = 'wind'
    hemis = 'nh'
    dfile = write_icedrift(out_dir,meth,hemis,dt0,dt1,dx,dy,)
    print ("Done with {}".format(dfile,))


if __name__ == '__main__':

    import argparse
    from argparse import RawDescriptionHelpFormatter

    p = argparse.ArgumentParser("icedrift_from_winds",formatter_class=RawDescriptionHelpFormatter)
    p.add_argument('DATE',help='date (YYYYMMDD) for running free-drift model')
    p.add_argument('--wind_dir',help='where to look for wind files.',default=None)
    p.add_argument('--out_dir',help='where to write obs and conf files.',default='./case')
    args = p.parse_args()

    if not os.path.exists(args.out_dir):
        os.mkdir(args.out_dir)

    dt = datetime.strptime(args.DATE,"%Y%m%d")
    compute_drift(dt, args.wind_dir, out_dir = args.out_dir,)

    sys.exit('Done with icedrift_from_winds')




