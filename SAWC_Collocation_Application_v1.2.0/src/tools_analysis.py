##########################################################################
#
# PYTHON 3 FUNCTIONS FOR tools_analysis
#
# CONTRIBUTOS:  Katherine E. Lukens             NOAA/NESDIS/STAR, CISESS at U. of Maryland
#               Kevin Garrett                   NOAA/NWS/OSTI
#               Kayo Ide                        U. of Maryland
#               David Santek                    CIMSS at U. of Wisconsin-Madison
#               Brett Hoover                    NOAA/NWS/NCEP/EMC, Lynker Technologies
#               David Huber                     NOAA/NWS/NCEP/EMC, Redline Performance Solutions, LLC
#               Ross N. Hoffman                 NOAA/NESDIS/STAR, CISESS at U. of Maryland
#               Hui Liu                         NOAA/NESDIS/STAR, CISESS at U. of Maryland
#
# Built with the following conda environment:
#
# name: bhoover-obs_match_3d
# channels:
#   - conda-forge
#   - defaults
# dependencies:
#   - python=3
#   - numpy
#   - pandas
#   - pynio
#   - matplotlib
#   - cartopy
#   - jupyter
#   - netCDF4
#   - scikit-learn
#   - dask
#   - geopy
#   - pip
#   - pip:
#
###########################################################################
#
# Import modules
#

import sys
import math
import numpy as np #....................................................... Array module
import datetime as dt #.................................................... Datetime module
import time #.............................................................. Time module
from datetime import datetime
import warnings

import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.colors import Normalize 
from scipy.interpolate import interpn
from cartopy.util import add_cyclic_point

from statistics import mean

#
###########################################################################
#
# STATISTICAL ANALYSIS functions for collocation program
#

fill = -999.0

# -------------------------------------------------------------------------
# Compute Horizontal Line-of-Sight (HLOS) Wind for using Aeolus azimuth angle
#
#	INPUTS:
#		x_azm ............................... Aeolus azimuth angle (HLOS direction)
#		x_spd ............................... Aeolus wind speed
#		y_dir ............................... Non-Aeolus wind direction
#		y_spd ............................... Non-Aeolus wind speed
#
#	OUTPUTS:
#		hlos ................................ Non-Aeolus wind projected onto Aeolus HLOS direction (x_azm)
#
def compute_hlos(x_azm,x_spd,y_dir,y_spd):

  if x_spd==np.nan or y_spd==np.nan:
    hlos = np.nan
    return hlos

	# compute HLOS wind velocity of y using x HLOS angle
  sindir = math.sin(y_dir*(math.pi/180.0))
  cosdir = math.cos(y_dir*(math.pi/180.0))
  
  sinazm = -1.0*math.sin(x_azm*(math.pi/180.0))
  cosazm = -1.0*math.cos(x_azm*(math.pi/180.0))
  
  u = -1.0*y_spd*sindir
  v = -1.0*y_spd*cosdir
  
  hlos = u*sinazm + v*cosazm
  
  return hlos
    
# -------------------------------------------------------------------------
# Reassign all missing values to NaN
#	Applies to collocation pairs using any dataset
#
#       INPUTS:
#               sdrv_spd ............................ Wind velocity array 
#
#       OUTPUTS:
#               drv_spd ............................. Wind velocity array with NaN as missing value
#
def to_nan(sdrv_spd):

    # Reassign array value to NaN if value <= fill value (missing)
    tdrv_spd  = [np.nan if sdrv_spd[i]==fill else sdrv_spd[i] for i in range(np.size(sdrv_spd))]
    del sdrv_spd

    drv_spd = np.asarray(tdrv_spd)

    return drv_spd
    
# -------------------------------------------------------------------------
# Super-ob (average) all Dependent observations that match the same Driver observation
#
def superob_matches(idxD,Dyr,Dmm,Ddy,Dhr,Dmn,Dlat,Dlon,Dprs,Dhgt,Dspd,Ddir,tlat,tlon,tprs,thgt,tspd,tdir,Dwcm,twcm):

    sqdrv_yr = []; sqdrv_mm = []; sqdrv_dy = []; sqdrv_hr = []; sqdrv_mn = []
    sqdrv_lat = []; sqdrv_lon = []; sqdrv_prs = []; sqdrv_hgt = []; sqdrv_spd = []; sqdrv_dir = []
    sqt_lat = []; sqt_lon = []; sqt_prs = []; sqt_hgt = []; sqt_spd = []; sqt_dir = []
    sqdrv_wcm = []; sqt_wcm = []
    
    idxD_uniq = list(set(idxD))
    for imax in range(np.size(idxD_uniq)):
      tmp_idx = np.where(idxD==idxD_uniq[imax])
      if np.size(tmp_idx)>0:
        ttlat = tlat[tmp_idx]
        ttlon = tlon[tmp_idx]
        ttprs = tprs[tmp_idx]
        tthgt = thgt[tmp_idx]
        ttspd = tspd[tmp_idx]
        ttdir = tdir[tmp_idx]
        ttwcm = twcm[tmp_idx]
		# average DEPENDENT obs per DRIVER ob.
        mean_lat = np.mean(ttlat)
        mean_lon = np.mean(ttlon)
        mean_prs = np.mean(ttprs)
        mean_hgt = np.mean(tthgt)
        mean_spd = np.mean(ttspd)
        mean_dir = np.mean(ttdir)
        mean_wcm = np.mean(ttwcm)
        del ttlat,ttlon,ttprs,tthgt,ttspd,ttdir,ttwcm
		# append to output arrays
			# DEPENDENT
        sqt_lat.append(mean_lat)
        sqt_lon.append(mean_lon)
        sqt_prs.append(mean_prs)
        sqt_hgt.append(mean_hgt)
        sqt_spd.append(mean_spd)
        sqt_dir.append(mean_dir)
        sqt_wcm.append(mean_wcm)
        del mean_lat,mean_lon,mean_prs,mean_hgt,mean_spd,mean_dir,mean_wcm
			# DRIVER
        sDyr  = Dyr[tmp_idx]        
        sDmm  = Dmm[tmp_idx]        
        sDdy  = Ddy[tmp_idx]        
        sDhr  = Dhr[tmp_idx]        
        sDmn  = Dmn[tmp_idx]        
        sDlat = Dlat[tmp_idx]        
        sDlon = Dlon[tmp_idx]        
        sDprs = Dprs[tmp_idx]        
        sDhgt = Dhgt[tmp_idx]        
        sDspd = Dspd[tmp_idx]        
        sDdir = Ddir[tmp_idx]        
        sDwcm = Dwcm[tmp_idx]        

        sqdrv_yr.append(sDyr[0])
        sqdrv_mm.append(sDmm[0])
        sqdrv_dy.append(sDdy[0])
        sqdrv_hr.append(sDhr[0])
        sqdrv_mn.append(sDmn[0])
        sqdrv_lat.append(sDlat[0])
        sqdrv_lon.append(sDlon[0])
        sqdrv_prs.append(sDprs[0])
        sqdrv_hgt.append(sDhgt[0])
        sqdrv_spd.append(sDspd[0])
        sqdrv_dir.append(sDdir[0])
        sqdrv_wcm.append(sDwcm[0])
        del sDyr,sDmm,sDdy,sDhr,sDmn,sDlat,sDlon,sDprs,sDhgt,sDspd,sDdir,sDwcm
      del tmp_idx

	# convert list to np.array
    qdrv_yr  = np.asarray(sqdrv_yr )
    qdrv_mm  = np.asarray(sqdrv_mm )
    qdrv_dy  = np.asarray(sqdrv_dy )
    qdrv_hr  = np.asarray(sqdrv_hr )
    qdrv_mn  = np.asarray(sqdrv_mn )
    qdrv_lat = np.asarray(sqdrv_lat)
    qdrv_lon = np.asarray(sqdrv_lon)
    qdrv_prs = np.asarray(sqdrv_prs)
    qdrv_hgt = np.asarray(sqdrv_hgt)
    qdrv_spd = np.asarray(sqdrv_spd)
    qdrv_dir = np.asarray(sqdrv_dir)
    qdrv_wcm = np.asarray(sqdrv_wcm)

    qt_lat   = np.asarray(sqt_lat)
    qt_lon   = np.asarray(sqt_lon)
    qt_prs   = np.asarray(sqt_prs)
    qt_hgt   = np.asarray(sqt_hgt)
    qt_spd   = np.asarray(sqt_spd)
    qt_dir   = np.asarray(sqt_dir)
    qt_wcm   = np.asarray(sqt_wcm)

    return qdrv_yr,qdrv_mm,qdrv_dy,qdrv_hr,qdrv_mn,qdrv_lat,qdrv_lon,qdrv_prs,qdrv_hgt,qdrv_spd,qdrv_dir,qt_lat,qt_lon,qt_prs,qt_hgt,qt_spd,qt_dir,qdrv_wcm,qt_wcm

# -------------------------------------------------------------------------
# Convert pressure to height
#	Convert pressure to pressure altitude (height) following NWS formulation (https://www.weather.gov/media/epz/wxcalc/pressureAltitude.pdf)
#
# 	Input	
#		prs .............................. Pressure in hPa
# 	Output 	
#		hgt .............................. Height in km
#
def prs_to_hgt(prs):	  

    	# check pressure units and convert to hPa
    if max(prs) > 10000.:
      prs = prs/100.

	# convert pressure to height
    hgt = np.nan * np.ones_like(prs)
    for i in range(np.size(prs)):
      hgt[i] = 145366.45 * (1.0 - (prs[i]/1013.25)**0.190284)			# convert hPa (mb) to feet
      hgt[i] = hgt[i] * 0.3048							# convert to meters
      hgt[i] = hgt[i] / 1000.0							# convert to km

    return hgt

# -------------------------------------------------------------------------
# Find number density of matched winds within grid cells as defined by (x,y), and conform resulting array to 2 dimensions
#
def var_to_2d_counts(xstr,ystr,tx,ty,txvar0,yvar0,txvar1,yvar1):

	# convert lists to arrays
    x     = np.asarray(tx)
    y     = np.asarray(ty)
    xvar0 = np.asarray(txvar0)		# DRIVER
    xvar1 = np.asarray(txvar1)		# DEPENDENT

        # sort y to be in ascending order
    y.sort()

        # loop to compute sums per grid cell (determined by x,y arrays)
		# 'nsum' contains the number of obs within each grid cell, as determined by
		# the conditions within np.where() in the loop.
    nsum0 = np.zeros([len(y),len(x)], dtype=float)
    nsum1 = np.zeros([len(y),len(x)], dtype=float)
    for j in range(len(y)-1):
      if ystr.find("Pressure")!=-1 or ystr.find("Height")!=-1:
        ydiff = abs(y[j+1] - y[j])/2.0                  # half the difference between each gridded z value
        ymin  = y[j]-ydiff
        ymax  = y[j]+ydiff
        del ydiff
      else:
        ymin  = y[j]
        ymax  = y[j+1]

		# 'xidx' is an array containing the number of obs that satisfy the conditions within np.where()
      xidx0 = np.zeros([len(x)], dtype=int)
      xidx1 = np.zeros([len(x)], dtype=int)

      if xstr.find("Time")!=-1:
        xidx0 = [np.size(np.where((xvar0>=x[k])*(yvar0>=ymin)*(yvar0<ymax))) for k in range(len(x))]
        xidx1 = [np.size(np.where((xvar1>=x[k])*(yvar1>=ymin)*(yvar1<ymax))) for k in range(len(x))]
      else:
        xidx0 = [np.size(np.where((xvar0>=x[k])*(xvar0<x[k+1])*(yvar0>=ymin)*(yvar0<ymax))) for k in range(len(x)-1)]
        xidx1 = [np.size(np.where((xvar1>=x[k])*(xvar1<x[k+1])*(yvar1>=ymin)*(yvar1<ymax))) for k in range(len(x)-1)]
			# compute density for cyclic grid point
        xidx0_last = [np.size(np.where((xvar0>=x[len(x)-1])*(yvar0>=ymin)*(yvar0<ymax)))]
        xidx1_last = [np.size(np.where((xvar1>=x[len(x)-1])*(yvar1>=ymin)*(yvar1<ymax)))]
        xidx0 = np.append(xidx0,xidx0_last,axis=0)
        xidx1 = np.append(xidx1,xidx1_last,axis=0)
        del xidx0_last,xidx1_last

      shapeidx0 = np.shape(xidx0)		# shape of xidx
      shapeidx1 = np.shape(xidx1)		# shape of xidx

      nsum0[j,0:shapeidx0[0]] = xidx0             # write 'xidx' to each j dimension of 'nsum'
      nsum1[j,0:shapeidx1[0]] = xidx1             # write 'xidx' to each j dimension of 'nsum'

      del xidx0,xidx1,shapeidx0,shapeidx1
      del ymin,ymax

        # mean per grid cell
    warnings.filterwarnings('ignore', category=RuntimeWarning)

	# fill grid cells where nsum=0 with nan
    nsum_nan0 = np.where(nsum0==0,np.nan,nsum0)
    nsum_nan1 = np.where(nsum1==0,np.nan,nsum1)

    return nsum_nan0,nsum_nan1   	# return counts per grid cell
    
# -------------------------------------------------------------------------
# Find observation number density for a single dataset within grid cells as defined by (x,y), and conform resulting array to 2 dimensions
#
def var_to_2d_counts_1dset(xstr,ystr,tx,ty,txvar0,tyvar0):

	# convert lists to arrays
    x     = np.asarray(tx)
    y     = np.asarray(ty)
    xvar0 = np.asarray(txvar0)
    yvar0 = np.asarray(tyvar0)

        # loop to compute sums per grid cell (determined by x,y arrays)
		# 'nsum' contains the number of obs within each grid cell, as determined by
		# the conditions within np.where() in the loop.
    nsum0 = np.zeros([len(y),len(x)], dtype=float)
    for j in range(len(y)-1):
      if ystr.find("Pressure")!=-1 or ystr.find("Height")!=-1:
        ydiff = abs(y[j+1] - y[j])/2.0                  # half the difference between each gridded z value
        ymin  = y[j]-ydiff
        ymax  = y[j]+ydiff
        del ydiff
      else:
        ymin  = y[j]
        ymax  = y[j+1]

		# 'xidx' is an array containing the number of obs that satisfy the conditions within np.where()
      xidx0 = np.zeros([len(x)], dtype=int)
      if xstr.find("Time")!=-1:
        xidx0 = [np.size(np.where((xvar0>=x[k])*(yvar0>=ymin)*(yvar0<ymax))) for k in range(len(x))]
      else:
        xidx0 = [np.size(np.where((xvar0>=x[k])*(xvar0<x[k+1])*(yvar0>=ymin)*(yvar0<ymax))) for k in range(len(x)-1)]
                        # compute density for cyclic grid point
        xidx0_last = [np.size(np.where((xvar0>=x[len(x)-1])*(yvar0>=ymin)*(yvar0<ymax)))]
        xidx0 = np.append(xidx0,xidx0_last,axis=0)
        del xidx0_last

      sizeidx0 = np.size(xidx0)			# size of xidx
      nsum0[j,0:sizeidx0] = xidx0		# write 'xidx' to each j dimension of 'nsum'

      del xidx0,sizeidx0
      del ymin,ymax

        # mean per grid cell
    warnings.filterwarnings('ignore', category=RuntimeWarning)

	# fill grid cells where nsum=0 with nan
    nsum_nan0 = np.where(nsum0==0,np.nan,nsum0)

    return nsum_nan0	   	# return counts per grid cell

# -------------------------------------------------------------------------
# Conform variable to 2 dimensions
#
# Input:
#	xstr .......................... name of new x dimension
#	ystr .......................... name of new y dimension
#	x ............................. new x dimension
#	y ............................. new y dimension
#	xdim .......................... 1d x dimension variable to be transformed
#	ydim ........................., 1d y dimension variable to be transformed
#	var0 .......................... variable from dataset 0 to be transformed to 2d grid (x,y)
#	var1 .......................... variable from dataset 1 to be transformed to 2d grid (x,y)
#
# Output:
#	mean0 ......................... mean of var0 on 2d grid
#	mean1 ......................... mean of var1 on 2d grid
#	meandiff ...................... mean difference (var1 - var0) on 2d grid
#	nSD_nan ....................... standard deviation of meandiff on 2d grid
#
def var_to_2d(xstr,ystr,x,y,xdim,ydim,var0,var1):

	# sort y to be in ascending order
    y.sort()

	# loop to compute sums per grid cell (determined by x,y arrays)
    SDdiff2d  = np.zeros([len(y),len(x)], dtype=float)
    sumdiff2d = np.zeros([len(y),len(x)], dtype=float)
    sumx2d    = np.zeros([len(y),len(x)], dtype=float)
    sumy2d    = np.zeros([len(y),len(x)], dtype=float)
    nsum      = np.zeros([len(y),len(x)], dtype=float)
    for j in range(len(y)-1):
      if ystr.find("Pressure")!=-1 or ystr.find("Height")!=-1:
        ydiff = abs(y[j+1] - y[j])/2.0  		# half the difference between each gridded z value
        ymin  = y[j]-ydiff
        ymax  = y[j]+ydiff
        del ydiff
      else:
        ymin  = y[j]
        ymax  = y[j+1]

      for k in range(len(x)-1):
        if xstr.find("Time")!=-1:
          xidx = np.where((xdim==x[k])*(ydim>=ymin)*(ydim<ymax))
        else:
          xidx = np.where((xdim>=x[k])*(xdim<x[k+1])*(ydim>=ymin)*(ydim<ymax))
        sumx2d[j,k]    = np.sum(var0[xidx])
        sumy2d[j,k]    = np.sum(var1[xidx])
        sumdiff2d[j,k] = np.sum(var1[xidx]-var0[xidx])
        SDdiff2d[j,k]  = np.std(var1[xidx]-var0[xidx])
        nsum[j,k]      = np.size(xidx)
        del xidx
                        # compute density for cyclic grid point
      xidx = np.where((xdim>=x[len(x)-1])*(ydim>=ymin)*(ydim<ymax))
      sumx2d[j,k]    += np.sum(var0[xidx])
      sumy2d[j,k]    += np.sum(var1[xidx])
      sumdiff2d[j,k] += np.sum(var1[xidx]-var0[xidx])
      SDdiff2d[j,k]  += np.std(var1[xidx]-var0[xidx])
      nsum[j,k]      += np.size(xidx)
      del xidx
      del ymin,ymax
	    
	# mean per grid cell 
    warnings.filterwarnings('ignore', category=RuntimeWarning)
    mean0    = sumx2d / nsum
    mean1    = sumy2d / nsum
    meandiff = sumdiff2d / nsum

    nSD_nan = SDdiff2d
    for j in range(len(y)):
      for k in range(np.size(x)):
        if nSD_nan[j,k]==0: nSD_nan[j,k]=np.nan

    return mean0,mean1,meandiff,nSD_nan 	# return mean x, mean y, mean diff, and SD of diff on 2d grid
    
# -------------------------------------------------------------------------
# Conform variable to 2 dimensions
#	One dataset only
#
# Input:
#       xstr .......................... name of new x dimension
#       ystr .......................... name of new y dimension
#       x ............................. new x dimension
#       y ............................. new y dimension
#       xdim .......................... 1d x dimension variable to be transformed
#       ydim ........................., 1d y dimension variable to be transformed
#       var0 .......................... variable from dataset to be transformed to 2d grid (x,y)
#
# Output:
#       mean0 ......................... mean of var0 on 2d grid
#
def var_to_2d_1dset(xstr,ystr,x,y,xdim,ydim,var0):

	# sort y to be in ascending order
    y.sort()

	# loop to compute sums per grid cell (determined by x,y arrays)
    sumx2d    = np.zeros([len(y),len(x)], dtype=float)
    nsum      = np.zeros([len(y),len(x)], dtype=float)
    for j in range(len(y)-1):
      if ystr.find("Pressure")!=-1 or ystr.find("Height")!=-1:
        ydiff = abs(y[j+1] - y[j])/2.0  		# half the difference between each gridded z value
        ymin  = y[j]-ydiff
        ymax  = y[j]+ydiff
        del ydiff
      else:
        ymin  = y[j]
        ymax  = y[j+1]

      for k in range(len(x)-1):
        if xstr.find("Time")!=-1:
          xidx = np.where((xdim==x[k])*(ydim>=ymin)*(ydim<ymax))
        else:
          xidx = np.where((xdim>=x[k])*(xdim<x[k+1])*(ydim>=ymin)*(ydim<ymax))
        sumx2d[j,k] = np.sum(var0[xidx])
        nsum[j,k]   = np.size(xidx)
        del xidx
                        # compute density for cyclic grid point
      xidx = np.where((xdim>=x[len(x)-1])*(ydim>=ymin)*(ydim<ymax))
      sumx2d[j,k] += np.sum(var0[xidx])
      nsum[j,k]   += np.size(xidx)
      del xidx
      del ymin,ymax
	    
	# mean per grid cell 
    warnings.filterwarnings('ignore', category=RuntimeWarning)
    mean0    = sumx2d / nsum

    return mean0 			# return mean x on 2D grid
      
# -------------------------------------------------------------------------
# Interpolate to given pressure level
#	Assumes a log-linear relationship
#
#	Input
#   		vcoord_data .......................... 1D array of vertical level values (e.g., pressure from a radiosonde)
#    		interp_var ........................... 1D array of the variable to be interpolated to all pressure levels
#    		interp_levels ........................ 1D array containing veritcal levels to interpolate to
#
#    	Output
#    		interp_data .......................... 1D array that contains the interpolated variable on the interp_levels
#
# Source: https://unidata.github.io/python-training/gallery/observational_data_cross_section/
#
def vert_interp(vcoord_data, interp_var, interp_levels):

    # Make veritcal coordinate data and grid level log variables
    lnp 	  = np.log(vcoord_data)
    lnp_intervals = np.log(interp_levels)

    # Use numpy to interpolate from observed levels to grid levels
    interp_data   = np.interp(lnp_intervals[::-1], lnp[::-1], interp_var[::-1])[::-1]

    # Mask for missing data (generally only near the surface)
    mask_low 		   = interp_levels > vcoord_data[0]
    mask_high 		   = interp_levels < vcoord_data[-1]
    interp_data[mask_low]  = interp_var[0]
    interp_data[mask_high] = interp_var[-1]

    return interp_data


#############################################################################
