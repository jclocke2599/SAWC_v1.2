
###########################################################################
#
# PYTHON 3 FUNCTIONS FOR match_dependencies
#
# CONTRIBUTORS: Katherine E. Lukens             NOAA/NESDIS/STAR, CISESS at U. of Maryland
#               Kevin Garrett                   NOAA/NWS/OSTI
#               Kayo Ide                        U. of Maryland
#               David Santek                    CIMSS at U. of Wisconsin-Madison
#               Brett Hoover                    NOAA/NWS/NCEP/EMC, Lynker Technologies
#               David Huber                     NOAA/NWS/NCEP/EMC, Redline Performance Solutions, LLC
#               Ross N. Hoffman                 NOAA/NESDIS/STAR, CISESS at U. of Maryland
#               Hui Liu                         NOAA/NESDIS/STAR, CISESS at U. of Maryland
#
# Developed by: Brett Hoover
#               David Santek
#
# Modified by:  Katherine E. Lukens
#
# History:
#       2021            B. Hoover, D. Santek    Created/developed program.
#       2023-04-14      K.E. Lukens             Finalized for upload to SAWC archive.
#	2023-07-15	D. Huber, K.E. Lukens	Optimized program to improve performance.
#
# Built with the following conda environment:
#
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
import numpy as np #....................................................... Array module
import dask.array as da #.................................................. Dask module
from dask.distributed import Client, LocalCluster #........................ Dask client modules (LocalCluster == runs on a local machine)
import datetime as dt #.................................................... Datetime module
from geopy.distance import geodesic #...................................... Geodesic distance module
import time #.............................................................. Time module
import netCDF4 #........................................................... netCDF4 module
from netCDF4 import Dataset #.............................................. netCDF4 Dataset submodule
#
###########################################################################
#
# Block mappable functions: Functions that can be applied to Dask arrays
# using da.map_blocks() to distribute the function across multiple
# processors. Used in obs_match_3d().
#
###########################################################################
#
# QUALITY CONTROL FUNCTIONS
#
# Fix longitude: Longitudes need to be in 0-to-360 format. This function
#                will search for any longitude < 0. (implying a -180-to-
#                180 format) and perform a correction to transform to a
#                0-to-360 format.
#
#    INPUTS:
#            lon ............................. longitudes [n]
#
#    OUTPUTS:
#            lon ............................. longitudes (corrected) [n]
#
#
def fix_lon(lon):
    idx = np.where(lon < 0.)
    lon[idx] = lon[idx] + 360.
    
    return lon
# Fix year: Some dates may include an hour-catetory with values > 23.
#           This can happen when, for example, observations in a dataset
#           cross over to the next day and an hour = 24 or 25 is used as
#           as a convention to signifity that the observation appears on
#           hour = 0 or 1 of the day following the time-stamp of the obs
#           file. These need to be corrected to the appropriate date-
#           stamp before being converted into a datetime variable.
#    INPUTS:
#            yr ...,.......................... year [n]
#            mm ...,.......................... month [n]
#            dy ...,.......................... day [n]
#            hr ...,.......................... hour [n]
#            mn ...,.......................... minute [n]
#
#    OUTPUTS:
#            yr .............................. year (corrected) [n]
#

def fix_yr(yr,mm,dy,hr,mn):

 idx = []

 for j in hr:
  if float(j) >= 24.: 
   idx.append(j)
    
#    idx = np.where(hr>=24) #.......................................... Indices where hr=24
  if (np.size(idx)>0):
       for i in idx[0]:
           d = dt.datetime( #......................................... Correct datetime (hr=23 + 1 hr)
                             int(yr[i]) , 
                             int(mm[i]) , 
                             int(dy[i]) , 
                             23         ,
                             int(mn[i])
                           ) + dt.timedelta(hours=int(hr[i])-23)
           yr[i] = float(d.year)
           mm[i] = float(d.month)
           dy[i] = float(d.day)
           hr[i] = float(d.hour)
           mn[i] = float(d.minute)
  return yr
# Fix month Some dates may include an hour-catetory with values > 23.
#           This can happen when, for example, observations in a dataset
#           cross over to the next day and an hour = 24 or 25 is used as
#           as a convention to signifity that the observation appears on
#           hour = 0 or 1 of the day following the time-stamp of the obs
#           file. These need to be corrected to the appropriate date-
#           stamp before being converted into a datetime variable.
#    INPUTS:
#            yr ...,.......................... year [n]
#            mm ...,.......................... month [n]
#            dy ...,.......................... day [n]
#            hr ...,.......................... hour [n]
#            mn ...,.......................... minute [n]
#
#    OUTPUTS:
#            mm .............................. month (corrected) [n]
#
def fix_mm(yr,mm,dy,hr,mn):
 for l in hr:
  if float(l) == 24.:
   idx.append(l)

 #   idx = np.where(hr==24.) #.......................................... Indices where hr=24
   if (np.size(idx)>0):
        for i in idx[0]:
            d = dt.datetime( #......................................... Correct datetime (hr=23 + 1 hr)
                             int(yr[i]) , 
                             int(mm[i]) , 
                             int(dy[i]) , 
                             23         ,
                             int(mn[i])
                           ) + dt.timedelta(hours=int(hr[i])-23)
            yr[i] = float(d.year)
            mm[i] = float(d.month)
            dy[i] = float(d.day)
            hr[i] = float(d.hour)
            mn[i] = float(d.minute)
 return mm
# Fix day:  Some dates may include an hour-catetory with values > 23.
#           This can happen when, for example, observations in a dataset
#           cross over to the next day and an hour = 24 or 25 is used as
#           as a convention to signifity that the observation appears on
#           hour = 0 or 1 of the day following the time-stamp of the obs
#           file. These need to be corrected to the appropriate date-
#           stamp before being converted into a datetime variable.
#    INPUTS:
#            yr ...,.......................... year [n]
#            mm ...,.......................... month [n]
#            dy ...,.......................... day [n]
#            hr ...,.......................... hour [n]
#            mn ...,.......................... minute [n]
#
#    OUTPUTS:
#            dy .............................. day (corrected) [n]
#
def fix_dy(yr,mm,dy,hr,mn):
  for k in hr:
   if float(k) == 24.:
    idx.append(k)

#    idx = np.where(hr==24.) #.......................................... Indices where hr=24
    if (np.size(idx)>0):
        for i in idx[0]:
            d = dt.datetime( #......................................... Correct datetime (hr=23 + 1 hr)
                             int(yr[i]) , 
                             int(mm[i]) , 
                             int(dy[i]) , 
                             23         ,
                             int(mn[i])
                           ) + dt.timedelta(hours=int(hr[i])-23)
            yr[i] = float(d.year)
            mm[i] = float(d.month)
            dy[i] = float(d.day)
            hr[i] = float(d.hour)
            mn[i] = float(d.minute)
  return dy
# Fix hour: Some dates may include an hour-catetory with values > 23.
#           This can happen when, for example, observations in a dataset
#           cross over to the next day and an hour = 24 or 25 is used as
#           as a convention to signifity that the observation appears on
#           hour = 0 or 1 of the day following the time-stamp of the obs
#           file. These need to be corrected to the appropriate date-
#           stamp before being converted into a datetime variable.
#    INPUTS:
#            yr ...,.......................... year [n]
#            mm ...,.......................... month [n]
#            dy ...,.......................... day [n]
#            hr ...,.......................... hour [n]
#            mn ...,.......................... minute [n]
#
#    OUTPUTS:
#            hr .............................. hour (corrected) [n]
#
def fix_hr(yr,mm,dy,hr,mn):
  for m in hr:
   if float(m) == 24.:
    idx.append(m)

#    idx = np.where(hr==24.) #.......................................... Indices where hr=24
    if (np.size(idx)>0):
        for i in idx[0]:
            d = dt.datetime( #......................................... Correct datetime (hr=23 + 1 hr)
                             int(yr[i]) , 
                             int(mm[i]) , 
                             int(dy[i]) , 
                             23         ,
                             int(mn[i])
                           ) + dt.timedelta(hours=int(hr[i])-23)
            yr[i] = float(d.year)
            mm[i] = float(d.month)
            dy[i] = float(d.day)
            hr[i] = float(d.hour)
            mn[i] = float(d.minute)
  return hr
# Fix minute: Some dates may include an hour-catetory with values > 23.
#           This can happen when, for example, observations in a dataset
#           cross over to the next day and an hour = 24 or 25 is used as
#           as a convention to signifity that the observation appears on
#           hour = 0 or 1 of the day following the time-stamp of the obs
#           file. These need to be corrected to the appropriate date-
#           stamp before being converted into a datetime variable.
#    INPUTS:
#            yr ...,.......................... year [n]
#            mm ...,.......................... month [n]
#            dy ...,.......................... day [n]
#            hr ...,.......................... hour [n]
#            mn ...,.......................... minute [n]
#
#    OUTPUTS:
#            mn .............................. minute (corrected) [n]
#
def fix_mn(yr,mm,dy,hr,mn):
  for n in hr:
   if float(n) == 24.:
    idx.append(n)

  #  idx = np.where(hr==24.) #.......................................... Indices where hr=24
    if (np.size(idx)>0):
        for i in idx[0]:
            d = dt.datetime( #......................................... Correct datetime (hr=23 + 1 hr)
                             int(yr[i]) , 
                             int(mm[i]) , 
                             int(dy[i]) , 
                             23         ,
                             int(mn[i])
                           ) + dt.timedelta(hours=int(hr[i])-23)
            yr[i] = float(d.year)
            mm[i] = float(d.month)
            dy[i] = float(d.day)
            hr[i] = float(d.hour)
            mn[i] = float(d.minute)
  return mn
#
###########################################################################
#
# Match observations: The algorithm nests matching criteria within one
# another to increase efficiency. The order is time --> pressure -->
# distance, such that the pressure check isn't performed unless the time
# check is passed, and the distance check isn't performed unless the
# pressure check is passed. The distance check is the most computationally
# expensive of the three checks, which is why it is placed in the inner-
# most nest.
#
#    INPUTS:
#            latd ............................ latitudes (Dask array)  [n_d]
#            lond ............................ longitudes (Dask array) [n_d]
#            prsd ............................ pressures (Dask array)  [n_d]
#            yrd ............................. years (Dask array)      [n_d]
#            mmd ............................. months (Dask array)     [n_d]
#            dyd ............................. days (Dask array)       [n_d]
#            hrd ............................. hours (Dask array)      [n_d]
#            mnd ............................. minutes (Dask array)    [n_d]
#            latn ............................ latitudes (np array)    [n_n]
#            lonn ............................ longitudes (np array)   [n_n]
#            prsn ............................ pressures (np array)    [n_n]
#            yrn ............................. years (np array)        [n_n]
#            mmn ............................. months (np array)       [n_n]
#            dyn ............................. days (np array)         [n_n]
#            hrn ............................. hours (np array)        [n_n]
#            mnn ............................. minutes (np array)      [n_n]
#            dst_max ......................... maximum distance (km)
#            prs_max ......................... maximum log10(P)-diff (hPa)
#            tim_max ......................... maximum time-diff (km)
#            n_max ........................... maximum number of matches
#    OUTPUTS:
#            matches ......................... matching indices of np-array to each index of Dask-array [n_d,n_max]
#
def match_obs( 
               latd ,    		# latitudes ....... dask-array
               lond ,    		# lontitudes ...... dask-array
               prsd ,    		# pressures ....... dask-array
               yrd  ,    		# years ........... dask-array
               mmd  ,    		# months .......... dask-array
               dyd  ,    		# days ............ dask-array
               hrd  ,    		# hours ........... dask-array
               mnd  ,    		# minutes ......... dask-array
               latn ,    		# latitudes drv ... numpy-array
               lonn ,    		# lontitudes drv .. numpy-array
               prsn ,    		# pressures drv ... numpy-array
               yrn  ,    		# years drv ....... numpy-array
               mmn  ,    		# months drv ...... numpy-array
               dyn  ,    		# days drv ........ numpy-array
               hrn  ,    		# hours drv ....... numpy-array
               mnn  ,    		# minutes drv ..... numpy-array
               tdst_max , 		# maximum distance (units=km)
               tprs_max , 		# maximum pressure-difference (units=log(hPa))
               ttim_max , 		# maximum time-difference (units=minutes)
               n_max   , 		# maximum number of matches reported per dask array ob (fewer allowed matches == more efficient) 
	       hgtd  ,   		# height .......... dask-array
	       hgtn  ,   		# height drv ...... numpy-array
	       thgt_max ,   		# maximum height-difference (units=km) 
					#	'matches'=matched indices, 
					#	'dependents'=integer values pertaining to each dependent dataset per matched index 
               parallelize_set_1 ,	# flag indicating which dataset is the dask-array (parallelized) (False=DRIVER, or True=DEPENDENTS)  
               usePH			# integer indicating which dependent datasets each array space belongs to
             ):
   
	# Create output arrays
    results = np.nan * np.ones((np.size(latd),n_max,6),dtype=latd.dtype)  #... List of matching indices, dim size = 6 (for the following 6 vars) (initialized to nan)
    matches = np.nan * np.ones((np.size(latd),n_max),dtype=latd.dtype)    #... List of numpy indices matching dask ob (initialized to nan)
    GCDs    = np.nan * np.ones((np.size(latd),n_max),dtype=latd.dtype)    #... List of great circle distances for matching dask ob (initialized to nan)
    DPs     = np.nan * np.ones((np.size(latd),n_max),dtype=latd.dtype)    #... List of pressure differences for matching dask ob (initialized to nan)
    DPslog  = np.nan * np.ones((np.size(latd),n_max),dtype=latd.dtype)    #... List of log10(pressure) differences for matching dask ob (initialized to nan)
    HTs     = np.nan * np.ones((np.size(latd),n_max),dtype=latd.dtype)    #... List of height differences for matching dask ob (initialized to nan)
    DTs     = np.nan * np.ones((np.size(latd),n_max),dtype=latd.dtype)    #... List of time differences for matching dask ob (initialized to nan)

    	# Loop through dask array (dependehap datasets)
    for i in range(np.size(latd)):
        m = [] #........................................................... List of numpy indices matching i-th dask ob (initialized to empty)
        mGCD = []
        mDP = []
        mDPlog = []
        mHT = []
        mDT = []
        # Loop through numpy array (driver dataset)
        for j in range(np.size(latn)):
	          #Check if lats are within 1 degree of each other. If not, do NOT compute GCdist.
          lat_check = abs(latn[j] - latd[i])

          if lat_check <= 1.0:
            # Time constraint
            t1 = dt.datetime( #........................ Datetime of i-th dask ob
                              int(yrd[i]) , 
                              int(mmd[i]) , 
                              int(dyd[i]) , 
                              int(hrd[i]) , 
                              int(mnd[i])
                            )
            t2 = dt.datetime( #........................ Datetime of j-th numpy ob
                              int(yrn[j]) , 
                              int(mmn[j]) , 
                              int(dyn[j]) , 
                              int(hrn[j]) , 
                              int(mnn[j])
                            )
            tim_diff = abs(t2-t1) #.................... Time-difference in seconds: 
				  #			abs() will allow tim_diff to be expressed the same regardless of whether t1 is before or after t2
            tim_diff = tim_diff.seconds / 60. #........ Time-difference in minutes

            if np.size(ttim_max)==np.size(latd):
              tim_max = ttim_max[i]
            elif np.size(ttim_max)==np.size(latn):
              tim_max = ttim_max[j]

		# get actual (not absolute) time difference (negative value indicates time before DRIVER time)
            tt1 = (t1 - dt.datetime(1970,1,1)).total_seconds()
            tt2 = (t2 - dt.datetime(1970,1,1)).total_seconds()
            if parallelize_set_1==True:
		# da are DEPENDENTS
              tmp = tt1 - tt2
            elif parallelize_set_1==False:
		# np are DEPENDENTS
              tmp = tt2 - tt1
            del tt1,tt2
            if tmp < 0:
              tim_diff_out = -1.0*tim_diff
            else:
              tim_diff_out = tim_diff
            del tmp

            if tim_diff <= tim_max:
              if (np.size(usePH)==np.size(prsd) and usePH[i]==0) or (np.size(usePH)==np.size(prsn) and usePH[j]==0):
		#Apply pressure difference constraint if:
		#	DRIVER and DEPENDENT have pressure variable (even if one/both also have height)

                	# assign temporary pressure variables
                if not np.isnan(prsn[j]):
                  t_prsn = prsn[j]
                else:
                  t_prsn = np.nan

                if not da.isnan(prsd[i]):
                  t_prsd = prsd[i]
                else:
                  t_prsd = np.nan

                # Pressure constraint
                if not np.isnan(t_prsn) and not da.isnan(t_prsd):
                  prs_diff = abs(np.log10(t_prsn)-np.log10(t_prsd)) #...... log10(Pressure) difference in hPa

                  if np.size(ttim_max)==np.size(latd):
                    prs_max = tprs_max[i]
                  elif np.size(ttim_max)==np.size(latn):
                    prs_max = tprs_max[j]

                  if parallelize_set_1==True:
                    prs_diff_out_log = np.log10(t_prsd) - np.log10(t_prsn)
                    prs_diff_out     = t_prsd - t_prsn
                  elif parallelize_set_1==False:
                    prs_diff_out_log = np.log10(t_prsn) - np.log10(t_prsd)
                    prs_diff_out     = t_prsn - t_prsd

                  if prs_diff <= prs_max:
                    # Distance constraint
                    dst_diff = geodesic( #................................. Distance between obs in km
                                         (latd[i] , lond[i]) , 
                                         (latn[j] , lonn[j])
                                       ).km

                    if np.size(ttim_max)==np.size(latd):
                      dst_max = tdst_max[i]
                    elif np.size(ttim_max)==np.size(latn):
                      dst_max = tdst_max[j]

                    if dst_diff <= dst_max:
                        # MATCH
                        m.append(j)
                        mGCD.append(dst_diff)
                        mDP.append(prs_diff_out)
                        mDPlog.append(prs_diff_out_log)
                        mHT.append(np.nan)
                        mDT.append(tim_diff_out)
	
              elif (np.size(usePH)==np.size(prsd) and usePH[i]==1) or (np.size(usePH)==np.size(prsn) and usePH[j]==1):
		#Apply height difference constraint if:
		#	DRIVER and DEPENDENT only have height variable, or
		#	DRIVER has height only, or
		#	DEPENDENT has height only

                # Check if height array is missing, and set temporary height variables ('t_hgtn' and 't_hgtd') for difference calc.
		#	This applies if dataset only has pressure variable
                if not np.isnan(hgtn[j]):
                  t_hgtn = hgtn[j]
                elif np.isnan(hgtn).all() and not np.isnan(prsn[j]):
			#if height is missing and pressure is not, convert pressure to pressure altitude (height)
			#  following NWS formulation (https://www.weather.gov/media/epz/wxcalc/pressureAltitude.pdf)
                  t_hgtn = 145366.45 * (1.0 - (prsn[j]/1013.25)**0.190284)	# units = feet
                  t_hgtn = t_hgtn * 0.3048					# convert to meters
                  t_hgtn = t_hgtn / 1000.0					# convert to km
                else:
                  t_hgtn = np.nan

                if not da.isnan(hgtd[i]):
                  t_hgtd = hgtd[i]
                elif da.isnan(hgtd).all() and not da.isnan(prsd[i]):
                        #if height is missing and pressure is not, convert pressure to pressure altitude (height)
                        #  following NWS formulation (https://www.weather.gov/media/epz/wxcalc/pressureAltitude.pdf)
                  t_hgtd = 145366.45 * (1.0 - (prsd[i]/1013.25)**0.190284)      # units = feet
                  t_hgtd = t_hgtd * 0.3048                                      # convert to meters
                  t_hgtd = t_hgtd / 1000.0					# convert to km
                else:
                   t_hgtd = np.nan

	      	# Height constraint
                if not np.isnan(t_hgtn) and not da.isnan(t_hgtd):
                  hgt_diff = abs(t_hgtn-t_hgtd) #...... height difference in km 

                  if np.size(ttim_max)==np.size(latd):
                    hgt_max = thgt_max[i]
                  elif np.size(ttim_max)==np.size(latn):
                    hgt_max = thgt_max[j]

                  if parallelize_set_1==True:
                    hgt_diff_out = t_hgtd - t_hgtn
                  elif parallelize_set_1==False:
                    hgt_diff_out = t_hgtn - t_hgtd

                  if hgt_diff <= hgt_max:
                    # Distance constraint
                    dst_diff = geodesic( #................................. Distance between obs in km
                                         (latd[i] , lond[i]) , 
                                         (latn[j] , lonn[j])
                                       ).km

                    if np.size(ttim_max)==np.size(latd):
                      dst_max = tdst_max[i]
                    elif np.size(ttim_max)==np.size(latn):
                      dst_max = tdst_max[j]

                    if dst_diff <= dst_max:
                        # MATCH
                        m.append(j)
                        mGCD.append(dst_diff)
                        mDP.append(np.nan)
                        mDPlog.append(np.nan)
                        mHT.append(hgt_diff_out)
                        mDT.append(tim_diff_out)

        # Each element of matches is a list of numpy array indices matching
        # the i-th observation of dask array (or an empty list for no
        # matches). For each match up to n_max, fill in i-th row of matches
        # array with numpy array match indices
        for k in range(min([n_max,len(m)])):
            matches[i,k] = m[k]
            GCDs[i,k] = mGCD[k]
            DPs[i,k] = mDP[k]
            DPslog[i,k] = mDPlog[k]
            HTs[i,k] = mHT[k]
            DTs[i,k] = mDT[k]
		# 'results' contains the above 6 arrays for each i
            results[i,k,0:6] = [m[k], mGCD[k], mDP[k], mDPlog[k], mHT[k], mDT[k]]

    #print("Did we make it to return?", results.shape)
    return results
#
###########################################################################
#
# PYTHON 3 FUNCTION: obs_match_3d
#
# DESCRIPTION
#
# This function will perform matching between 2 observational datasets,
# producing a set of two lists containing index numbers matching obs
# from set-1 with obs from set-2. Multiple matches from one ob are
# allowed. The matching is based on 3 criteria (hence, 3d match):
#
#    1) Obs must sample within tim_max minutes of each other
#    2) Obs must sample within prs_max log10(pressure) of each other,
#       this constraint is tighter at lower pressures
#    3) Obs must sample within dst_max km of each other
#
# This function assumes set-1 is the driver, meaning that every ob in
# set-1 is independently compared to all obs in set-2 to define matches.
# Based on this architecture, we can represent set-1 as chunked arrays
# that are evaluated with a distributed computing method (Dask). Since
# all obs in set-2 must be used for comparison, set-2 is not represented
# as a chunked array and instead is treated as a single numPy array (or
# an array with one chunk, as Dask interprets it). The chunking of set-1
# allows for scalability of the matching function across processors.
#
# INPUTS
#
# This function is intended to be used to perform matches on datasets that
# may be distributed across multiple input files. So the inputs to this
# function are the full arrays of data for set-1 and set-2, rather than
# file names. A function can be called prior to obs_match_3d to collect
# data from multiple input files and pass the total arrays to obs_match_3d.
#
# Input variables
#
# lat_1_np = Latitudes of set-1 obs (single-dimension numPy array.....deg)
# lon_1_np = Longitudes of set-1 obs (single-dimension numPy array....deg)
# prs_1_np = Pressure of set-1 obs (single-dimension numPy array......hPa)
# yr_1_np  = Year of set-1 obs (single-dimension numPy array..........YYYY)
# mm_1_np  = Month of set-1 obs (single-dimension numPy array.........MM)
# dy_1_np  = Day of set-1 obs (single-dimension numPy array...........DD)
# hr_1_np  = Hour of set-1 obs (single-dimension numPy array..........HH)
# mn_1_np  = Minute of set-1 obs (single-dimension numPy array........MN)
# lat_2_np = Latitudes of set-2 obs (single-dimension numPy array.....deg)
# lon_2_np = Longitudes of set-2 obs (single-dimension numPy array....deg)
# prs_2_np = Pressure of set-2 obs (single-dimension numPy array......hPa)
# yr_2_np  = Year of set-2 obs (single-dimension numPy array..........YYYY)
# mm_2_np  = Month of set-2 obs (single-dimension numPy array.........MM)
# dy_2_np  = Day of set-2 obs (single-dimension numPy array...........DD)
# hr_2_np  = Hour of set-2 obs (single-dimension numPy array..........HH)
# mn_2_np  = Minute of set-2 obs (single-dimension numPy array........MN)
# dst_max  = Maximum allowable distance between matches...............km
# prs_max  = Maximum allowable log10(P) distance between matches......hPa
# tim_max  = Maximum allowable time between matches...................mins
# n_max    = Maximum allowable number of matches to search set........ (int)
# nproc    = Number of processors to use..............................(int)
#
# OUTPUTS
#
# The function outputs two lists corresponding to index numbers of set-1
# and set-2 obs that pass the 3d match. The i-th element of one list is
# matched to the i-th element of the other list, establishing a match.
# Duplicates can appear in either list, representing multiple matches of
# an ob in a list with obs in the other list, allowing for the possibility
# of superobbing.
#
# Output Variables
#
# match_1_list = List of set-1 indices that match set-2
# match_2_list = List of set-1 indices that match set-1
#
###########################################################################
def obs_match_3d(
                  lat_1_np ,
                  lon_1_np ,
                  prs_1_np ,
                  yr_1_np  ,
                  mm_1_np  ,
                  dy_1_np  ,
                  hr_1_np  ,
                  mn_1_np  ,
                  lat_2_np ,
                  lon_2_np ,
                  prs_2_np ,
                  yr_2_np  ,
                  mm_2_np  ,
                  dy_2_np  ,
                  hr_2_np  ,
                  mn_2_np  ,
                  dst_max  ,
                  prs_max  ,
                  tim_max  ,
                  n_max    ,
                  nproc    ,
		  hgt_1_np ,
		  hgt_2_np , 
		  hgt_max  ,
                  usePH
                ):
		# NOTE: *_1_np = dependent datasets | *_2_np = driver dataset 
    
    #######################################################################
    #
    master_clock_begin = time.time() #..................................... Beginning-time of function-call
    #
    #######################################################################
    #
    # 1. OPEN A DASK CLIENT to schedule tasks
    #
    nthreads = 1
    cluster = LocalCluster(n_workers=nproc,threads_per_worker=nthreads) #......... Local cluster definition with nproc single-thread processors
    client = Client(cluster) #............................................. Dask client
    #
    #######################################################################
    #
    # 2. DATA INGEST AND QUALITY CONTROL
    #
    clock1 = time.time() #................................................. Beginning-time of task
    #
    # Convert larger of set-1 or set-2 to Dask arrays: This will optimize
    # parallelization
    #
    # Define chunk-size of parallelized arrays by dividing (as equally as possible)
    # the data in the parallelized set among nproc processors
    #
    print("Beginning chunk size") 
    if (np.size(lat_1_np) >= np.size(lat_2_np)):
        parallelize_set_1 = True #......................................... Flag: set-1 is dask, set-2 is numpy
        chksize=int(np.ceil(np.size(lat_1_np)/nproc)) #.................... Chunk-size
        print('DEPENDENT dataset(s) parallelized')
        # Chunk set-1 arrays into Dask arrays
        lat_da = da.from_array(lat_1_np,chunks=chksize) #.................. Dask array version of lat_1_np
        lon_da = da.from_array(lon_1_np,chunks=chksize) #.................. Dask array version of lon_1_np
        prs_da = da.from_array(prs_1_np,chunks=chksize) #.................. Dask array version of prs_1_np
        yr_da  = da.from_array(yr_1_np,chunks=chksize) #................... Dask array version of yr_1_np
        mm_da  = da.from_array(mm_1_np,chunks=chksize) #................... Dask array version of mm_1_np
        dy_da  = da.from_array(dy_1_np,chunks=chksize) #................... Dask array version of dy_1_np
        hr_da  = da.from_array(hr_1_np,chunks=chksize) #................... Dask array version of hr_1_np
        mn_da  = da.from_array(mn_1_np,chunks=chksize) #................... Dask array version of mn_1_np
        hgt_da = da.from_array(hgt_1_np,chunks=chksize) #.................. Dask array version of hgt_1_np
        usePH_chunk = da.from_array(usePH,chunks=chksize) #................ Dask array version of usePH

        dst_max_npda = da.from_array(dst_max,chunks=chksize)	#.......... Dask array version of dst_max
        prs_max_npda = da.from_array(prs_max,chunks=chksize)	#.......... Dask array version of prs_max
        tim_max_npda = da.from_array(tim_max,chunks=chksize)	#.......... Dask array version of tim_max
        hgt_max_npda = da.from_array(hgt_max,chunks=chksize)	#.......... Dask array version of hgt_max
        # Set-2 arrays remain numPy arrays
        lat_np = lat_2_np #................................................ lat_2_np (renamed)
        lon_np = lon_2_np #................................................ lon_2_np (renamed)
        prs_np = prs_2_np #................................................ prs_2_np (renamed)
        yr_np = yr_2_np #.................................................. yr_2_np (renamed)
        mm_np = mm_2_np #.................................................. mm_2_np (renamed)
        dy_np = dy_2_np #.................................................. dy_2_np (renamed)
        hr_np = hr_2_np #.................................................. hr_2_np (renamed)
        mn_np = mn_2_np #.................................................. mn_2_np (renamed)
        hgt_np = hgt_2_np #................................................ hgt_2_np (renamed)
        # Delete old variables from memory
        del lat_1_np, lon_1_np, prs_1_np, yr_1_np, mm_1_np, dy_1_np, hr_1_np, mn_1_np
        del lat_2_np, lon_2_np, prs_2_np, yr_2_np, mm_2_np, dy_2_np, hr_2_np, mn_2_np
        del hgt_1_np, hgt_2_np 
        del usePH 
        del dst_max, prs_max, tim_max, hgt_max	
    else:
        parallelize_set_1 = False #........................................ Flag: set-1 is numpy, set-2 is dask
        chksize=int(np.ceil(np.size(lat_2_np)/nproc)) #.................... Chunk-size
        print('DRIVER dataset parallelized')
        # Chunk set-2 arrays into Dask arrays
        lat_da = da.from_array(lat_2_np,chunks=chksize) #.................. Dask array version of lat_2_np
        lon_da = da.from_array(lon_2_np,chunks=chksize) #.................. Dask array version of lon_2_np
        prs_da = da.from_array(prs_2_np,chunks=chksize) #.................. Dask array version of prs_2_np
        yr_da  = da.from_array(yr_2_np,chunks=chksize) #................... Dask array version of yr_2_np
        mm_da  = da.from_array(mm_2_np,chunks=chksize) #................... Dask array version of mm_2_np
        dy_da  = da.from_array(dy_2_np,chunks=chksize) #................... Dask array version of dy_2_np
        hr_da  = da.from_array(hr_2_np,chunks=chksize) #................... Dask array version of hr_2_np
        mn_da  = da.from_array(mn_2_np,chunks=chksize) #................... Dask array version of mn_2_np
        hgt_da = da.from_array(hgt_2_np,chunks=chksize) #.................. Dask array version of hgt_2_np
        # Set-2 arrays remain numPy arrays
        lat_np = lat_1_np #................................................ lat_1_np (renamed)
        lon_np = lon_1_np #................................................ lon_1_np (renamed)
        prs_np = prs_1_np #................................................ prs_1_np (renamed)
        yr_np = yr_1_np #.................................................. yr_1_np (renamed)
        mm_np = mm_1_np #.................................................. mm_1_np (renamed)
        dy_np = dy_1_np #.................................................. dy_1_np (renamed)
        hr_np = hr_1_np #.................................................. hr_1_np (renamed)
        mn_np = mn_1_np #.................................................. mn_1_np (renamed)
        hgt_np = hgt_1_np #................................................ hgt_1_np (renamed)
        usePH_chunk = usePH #.............................................. usePH (renamed)

        dst_max_npda = dst_max	
        prs_max_npda = prs_max	
        tim_max_npda = tim_max	
        hgt_max_npda = hgt_max	
        # Delete old variables from memory
        del lat_1_np, lon_1_np, prs_1_np, yr_1_np, mm_1_np, dy_1_np, hr_1_np, mn_1_np
        del lat_2_np, lon_2_np, prs_2_np, yr_2_np, mm_2_np, dy_2_np, hr_2_np, mn_2_np
        del hgt_1_np, hgt_2_np	
        del usePH 
        del dst_max, prs_max, tim_max, hgt_max	
    #
    print('Chunk size:',chksize)
    print('Maximum allowable matches:',n_max)
    #
    # Quality Control: Convert longitude to 0-to-360 degree format, if 
    # necessary.
    #
    lon_da = da.map_blocks(fix_lon,lon_da)
    lon_np = fix_lon(lon_np)
    #
    # Quality Control: Convert hour>=24 dates to corresponding dates
    # with day incremented
    #
    # Dask doesn't play nice with functions that return multiple variables
    # (these return as either tuples or lists, and Dask chokes on them 
    # because they contain no dtype). So functions are defined to return
    # the corrected year, month, day, hour, and minute individually.
    yr_fix = da.map_blocks( #.............................................. Corrected yr_da
                             fix_yr , yr_da , mm_da , 
                             dy_da  , hr_da , mn_da
                            )
    mm_fix = da.map_blocks( #.............................................. Corrected mm_da
                             fix_mm , yr_da , mm_da , 
                             dy_da  , hr_da , mn_da
                            )
    dy_fix = da.map_blocks( #.............................................. Corrected dy_da
                             fix_dy , yr_da , mm_da , 
                             dy_da  , hr_da , mn_da
                            )
    hr_fix = da.map_blocks( #.............................................. Corrected hr_da
                             fix_hr , yr_da , mm_da , 
                             dy_da  , hr_da , mn_da
                            )
    mn_fix = da.map_blocks( #.............................................. Corrected mn_da
                             fix_mn , yr_da , mm_da , 
                             dy_da  , hr_da , mn_da
                            )
    yr_da = yr_fix
    mm_da = mm_fix
    dy_da = dy_fix
    hr_da = hr_fix
    mn_da = mn_fix
    #
    del yr_fix, mm_fix, dy_fix, hr_fix, mn_fix
    #
    # Apply date-fix functions to numpy arrays directly. Use same method
    # to write to new variables and move over at the end.
    yr_fix = fix_yr(yr_np,mm_np,dy_np,hr_np,mn_np) #....................... Corrected yr_np
    mm_fix = fix_mm(yr_np,mm_np,dy_np,hr_np,mn_np) #....................... Corrected mm_np
    dy_fix = fix_dy(yr_np,mm_np,dy_np,hr_np,mn_np) #....................... Corrected dy_np
    hr_fix = fix_hr(yr_np,mm_np,dy_np,hr_np,mn_np) #....................... Corrected hr_np
    mn_fix = fix_mn(yr_np,mm_np,dy_np,hr_np,mn_np) #....................... Corrected mn_np

    yr_np = yr_fix
    mm_np = mm_fix
    dy_np = dy_fix
    hr_np = hr_fix
    mn_np = mn_fix
    #
    del yr_fix, mm_fix, dy_fix, hr_fix, mn_fix
    #
    clock2 = time.time() #................................................. Ending-time of task
    clock = clock2 - clock1 #.............................................. Clock time (sec)
    print('Data ingest completed in ',clock,' seconds')
    #
    #######################################################################
    #
    # 3. APPLY MATCHING ALGORITHM
    #
    clock1 = time.time() #................................................. Beginning-time of task
    #
    #
    # Apply matching function via da.map_blocks
    #	'result' contains 6 arrays returned by 'match_obs' function
    #
    print("CALL match_obs")
    result = da.map_blocks( #............................................. Matched numpy indices (uncomputed)
                             match_obs , 
                             lat_da    , 
                             lon_da    , 
                             prs_da    , 
                             yr_da     , 
                             mm_da     , 
                             dy_da     , 
                             hr_da     , 
                             mn_da     , 
                             lat_np    , 
                             lon_np    , 
                             prs_np    , 
                             yr_np     , 
                             mm_np     , 
                             dy_np     , 
                             hr_np     , 
                             mn_np     , 
                             dst_max_npda   , 
                             prs_max_npda   , 
                             tim_max_npda   ,
                             n_max     ,
			     hgt_da    ,
			     hgt_np    ,
			     hgt_max_npda   ,
                             parallelize_set_1,
                             usePH_chunk
                           )

    #
    clock2 = time.time() #................................................. Ending-time of task
    clock = clock2 - clock1 #.............................................. Clock time (sec)
    print('Matching algorithm defined in ',clock,'seconds')
    #
    #######################################################################
    #
    # 4. COMPUTE MATCHING
    #
    #    Dask carries out task on individual chunks of dask arrays treating
    #    numpy arrays as a one-chunk array.
    #
    clock1 = time.time() #................................................. Beginning-time of task
    results_np=result.compute(num_workers=1) #......................... Matched indices of numpy array for each ob in dask array (numPy array)
    clock2 = time.time() #................................................. Ending-time of task
    clock = clock2-clock1
    print('Match computing completed in ',clock,'seconds')

	# separate 'results_np' (from 'result') into each of its 6 arrays
    ob_idx_np = results_np[:,:,0]
    GCD_np = results_np[:,:,1]
    DP_np = results_np[:,:,2]
    DPlog_np = results_np[:,:,3]
    HT_np = results_np[:,:,4]
    DT_np = results_np[:,:,5]

    #######################################################################
    #
    # 5. COMPOSE MATCHING INDICES LISTS
    #
    clock1 = time.time() #................................................. Beginning-time of task
    da_idx = [] #.......................................................... List of dask array matching indices (initialized to empty)
    np_idx = [] #.......................................................... List of numpy array matching indices (initialized to empty)
    GCD_match = []
    DP_match = []
    DPlog_match = []
    HT_match = []
    DT_match = []
    for i in range(np.ma.size(ob_idx_np,axis=0)):
        # Extract row of matches to dask array ob-i
        x = ob_idx_np[i,:]
        yGCD = GCD_np[i,:]
        yDP = DP_np[i,:]
        yDPlog = DPlog_np[i,:]
        yHT = HT_np[i,:]
        yDT = DT_np[i,:]

        # Loop though row, appending i to da_idx and row-value to np_idx
        # for all row-values that are not np.nan
        #
        # Alert user if n_max matches appeared (may be in a situation where
        # increasing n_max is appropriate)
        if np.size(np.where(np.isnan(x))) == 0:
            print('WARNING: Maximum searchable matches found for observation ',i,' consider increasing n_max.')
        for match in x:
            if not np.isnan(match):
                da_idx.append(i)
                np_idx.append(int(match))
        for match_dim in range(len(yGCD)):
            if not np.isnan(yGCD[match_dim]):
                GCD_match.append(yGCD[match_dim])
                DP_match.append(yDP[match_dim])
                DPlog_match.append(yDPlog[match_dim])
                HT_match.append(yHT[match_dim])
                DT_match.append(yDT[match_dim])

    clock2 = time.time() #................................................. Ending-time of task
    clock = clock2-clock1
    print('Match lists completed in ',clock,'seconds')
    #
    #######################################################################
    #
    master_clock_end = time.time() #....................................... Ending-time of function-call
    master_clock = master_clock_end-master_clock_begin
    print('TOTAL FUNCTION TIME: ',master_clock,'seconds');
    #
    #######################################################################
    #
    # RETURN: 
    #
    # set_1_idx = list of indices of set-1 matching set-2 obs
    # set_2_idx = list of indices of set-2 matching set-1 obs
    #
    # Whether set 1 or 2 corresponds to da_idx or np_idx depends on which
    # set was chunked
    #
    if parallelize_set_1:
        set_1_idx = da_idx
        set_2_idx = np_idx
    else:
        set_1_idx = np_idx
        set_2_idx = da_idx
    #
    return set_1_idx, set_2_idx, GCD_match, DP_match, DPlog_match, HT_match, DT_match
#
###########################################################################
#
# Output match lists to netCDF file: indices from each dataset in vals list
# are written to a netCDF file with a name corresponding to the names list.
#
#    INPUTS:
#            names ........................... list of variable names
#            vals ............................ list of index lists
#            nc_out_filenm ................... netCDF output file name
#	     src_files_list .................... list of source files: full paths + filenames of Driver file(s) and dependent Dataset files used
#	     qc_list ......................... list of QC used, if applicable
#    OUTPUTS:
#            None
#
def match_lists_to_netCDF(names,vals,nc_out_filenm,longnames,units,missing,src_path,qcflag,qclist,qcnotes,shortnames,aeol_wind_type,aeol_dset_type):
    nc_out = Dataset( #.................................................... Dataset object for output
                      nc_out_filenm    , # Dataset input: Output file name
                      "w"              , # Dataset input: Make file write-able
                      format="NETCDF4" , # Dataset input: Set output format to netCDF4
                    )

    fill = -999

    # Add File Attributes
    nc_out.title = "Indices of Collocated Winds"

	# assign date created att
    nc_out.creation_date = dt.datetime.now().strftime('%Y-%m-%d %H:%M:%S')

    # Add Data
    n_vars = len(names) #.................................................. Number of variables to output
    for i in range(n_vars):
	
      # Set dimension name
      idx_name 		= "idx_"+str(i)

      # Set main attributes
      var_name 		= names[i]
      var_longname 	= longnames[i]
      var_units 	= units[i]
      var_miss 		= missing[i]
      
      if var_name == 'drv' or var_name.startswith('dset'):
        vartype = 'S'+str(len(vals[i]))	# character string of length = len(var_vals)

        var_vals	= np.array([vals[i]],dtype=vartype)

      	# Dataset name variables (strings)
		# Set additional attributes
        var_shortname   = shortnames[i]
        if (var_shortname.find('Aeolus') != -1):
          var_aeol_wind_type	= aeol_wind_type[i]
          var_aeol_dset_type	= aeol_dset_type[i]

        var_path 	= src_path[i]
        var_qcflag 	= qcflag[i]
        var_qclist 	= qclist[i]
        var_qcnotes	= qcnotes[i]

	# Create dimensions
        nc_out.createDimension('nstrings_'+str(i),None)
        nc_out.createDimension('nchars_'+str(i),len(vals[i]))

        # Create output variable
        var_out_str = nc_out.createVariable(var_name,"S1",('nstrings_'+str(i),'nchars_'+str(i)))

	# Fill output variable with data
        var_out_str[:] = netCDF4.stringtochar(var_vals)

        # Add attributes to each variable
        var_out_str.short_name    = var_shortname
        if (var_shortname.find('Aeolus') != -1):
          var_out_str.wind_type	  	= var_aeol_wind_type
          var_out_str.dataset_type	= var_aeol_dset_type
        var_out_str.long_name     = var_longname
        var_out_str.missing_value = var_miss
        var_out_str.source_file   = var_path
        var_out_str.QC_flag       = var_qcflag
        var_out_str.QC_list       = var_qclist
        var_out_str.QC_notes      = var_qcnotes

      elif var_name.startswith('idx_') or var_name.startswith('ntypes_') or var_name.startswith('nprofiles_') or var_name.startswith('nlevels_'):
        # Index or profile-specific variables (integers)
        vartype = "i8"		# 64-bit signed integer
      else:
	# All other variables (float)
        vartype = "f8"		# 64-bit floating-point number

      if vartype == "i8" or vartype == "f8":
        # Add dimension
        idx_out  = nc_out.createDimension( #............................... Output dimension
                                           idx_name , # nc_out.createDimension input: Dimension name
                                           None    # nc_out.createDimension input: Dimension size limit ("None" == unlimited)
                                         )

        # Output Variables
        var_out = nc_out.createVariable( #................................. Output variable
                                      var_name  , # nc_out.createVariable input: Variable name 
                                      vartype   , # nc_out.createVariable input: Variable format 
                                      ( 
                                        idx_name  # nc_out.createVariable input: Variable dimension
                                      )
                                    )

        var_vals = vals[i]

        # Fill netCDF file variables
        var_out[:] = var_vals
    
      		# Add attributes to each variable
        var_out.long_name       = var_longname
        var_out.missing_value   = var_miss

        if var_units != fill and var_units != str(fill):#not np.isnan(str(var_units)):
          var_out.units	      = var_units

    # Close netCDF file
    nc_out.close()
    return

###########################################################################

