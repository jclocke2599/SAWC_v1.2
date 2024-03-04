
###########################################################################
#
# PYTHON 3 FUNCTIONS FOR read_data
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
# Import python modules
#

from os.path import exists
import sys
import numpy as np #....................................................... Array module
import datetime as dt #.................................................... Datetime module
import time #.............................................................. Time module
from netCDF4 import Dataset #.............................................. netCDF module
import math

from quality_controls import qc_aeolus
#from quality_controls import qc_aircraft
from quality_controls import qc_amv

#
###########################################################################

# FUNCTIONS for reading input datasets

#===============================================================================================
# Check if pressure or height variable exists in NetCDF index file
#
#	INPUTS:
#		dset_file ........................... Index filename
#		Pvar ................................ Pressure differences variable name
#		Hvar ................................ Height differences variable name
#
#	OUTPUTS:
#		PHmatch ............................. Array of either pressure differences or height differences
#		PHstr ............................... String indicating the variable that exists (i.e., Pressure or Height)
#
def check_PHvar_exists(dset_file,Pvar,Hvar):

	# load dataset
  data_hdl = Dataset(dset_file)

  try:
  	# check if pressure diff array exists
    PHmatch	= np.asarray( data_hdl.variables[Pvar] )
    PHstr	= "Pressure"
  except:	
  	# if pressure diff does NOT exist
    try:	
    	# check if height diff array exists
      PHmatch	= np.asarray( data_hdl.variables[Hvar] )
      PHstr	= "Height"
    except: 
      print("ERROR in check_PHvar_exists: DP_match and HT_match do not exist! At least one must exist to continue.")
      sys.exit()

  data_hdl.close()

	# Return variable
  return PHmatch,PHstr

#===============================================================================================
# Read Collocation Index File
#
#	INPUTS:
#		yyyymmddhh .......................... Current date in yyyymmddhh format
#		pct ................................. Minimum AMV quality indicator (QI) in % for QC
#		bool_qc ............................. Choice to apply Aeolus QC: True=apply QC, False=don't apply QC
#
#	OUTPUTS:
#		d_lat ............................... Latitude in degrees [-90,90]
#		d_lon ............................... Longitude in degrees [0,360]
#		d_prs ............................... Pressure in hPa
#		d_hgt ............................... Height in km
#		d_yr ................................ Year
#		d_mm ................................ Month
#		d_dy ................................ Day
#		d_hr ................................ Hour
#		d_mn ................................ Minute
#		d_pccf .............................. AMV QI (percent confidence) in %
#               indexes2 ............................ Indices of input obs
#		qc_list ............................. List of QC applied (if applicable)
#
def read_index_file(path_prefix,yyyymmddhh,idx_file_str,dir_in):

  qc_list = ""			#initialize

  yyyy = yyyymmddhh[0:4]
  mm   = yyyymmddhh[4:6]
  dd   = yyyymmddhh[6:8]

	#-------------------------------------------------
    	# Define dataset

  if dir_in == str(path_prefix)+'collocation/index_files':	# index files are located within 'path_prefix'
    dset_path = path_prefix+'/collocation/index_files/'+yyyy+'/'+mm+'/'+dd+'/'
  else:								# index files are somewhere else
    dset_path = dir_in 

  dset_filename 	= 'index.'+yyyymmddhh+idx_file_str

  dset_file   		= dset_path+dset_filename
  print("read_index_file: index file = "+dset_file)
  dset_exists = exists(dset_file)
  if dset_exists==False:
    print("ERROR in read_index_file: index file "+dset_file+" does not exist!")
    dset_shortname = "NoMatches"
    idx_drv_dset   = -1
    idx_dset       = -1
    DT_match       = 0
    GCD_match      = 0
    Vert_match     = 0
    PHstr          = "NO_MATCHES"

    return dset_shortname,idx_drv_dset,idx_dset,DT_match,GCD_match,Vert_match,PHstr

    	# Path on FTP/web archive server (for output NetCDF only)
  str_dset_path = dset_path

  dset_src = str_dset_path+dset_filename

	# Variables names
		# DEPENDENT dataset name
  dset_name_var = 'dset1'

    		# indices of matched obs
			# DRIVER indices
  idx_drv_dset_var = 'idx_drv_'
  			# DEPENDENT indices
  idx_dset_var     = 'idx_'

		# collocation differences per matched pair
  DT_match_var  = 'DT_match_drv_'
  DP_match_var  = 'DP_match_drv_'
  DPlog_match_var  = 'DPlog_match_drv_'
  HT_match_var  = 'HT_match_drv_'
  GCD_match_var = 'GCD_match_drv_'
  
    	#-------------------------------------------------
  	# Load file and extract data
  print("read_index_file: load index file")

	# load dataset
  data_hdl = Dataset(dset_file)

  	# Extract data 
		# DEPENDENT dataset names
  		# ... get dset names
  dset_name      = np.asarray( data_hdl.variables[dset_name_var] )
  		# ... LOOP through attribute names
  for attname in data_hdl.variables[dset_name_var].ncattrs():
    if attname == 'short_name':
      dset_shortname = getattr(data_hdl.variables[dset_name_var],attname)

		# indices of matched obs
  		# ... get indices of matched obs
  		# ... ... DRIVER
  idx_drv_dset  = np.asarray( data_hdl.variables[idx_drv_dset_var+dset_name_var] )
  		# ... ... DEPENDENT(S)
  idx_dset      = np.asarray( data_hdl.variables[idx_dset_var+dset_name_var] )

		# collocation differences per matched pair
  		# ... get colloc criteria arrays
  if np.size(idx_drv_dset) > 0:
		# ... ... time difference
    DT_match		= np.asarray( data_hdl.variables[DT_match_var+dset_name_var] )
		# ... ... great circle distance
    GCD_match		= np.asarray( data_hdl.variables[GCD_match_var+dset_name_var] )  
		# ... ... pressure/height difference
    Vert_match,PHstr	= check_PHvar_exists(dset_file,DP_match_var+dset_name_var,HT_match_var+dset_name_var)
  else:
    DT_match 	= [-1]
    GCD_match 	= [-1]
    Vert_match 	= [-1]
    PHstr 	= "NO_MATCHES"

  data_hdl.close()

  return dset_shortname,idx_drv_dset,idx_dset,DT_match,GCD_match,Vert_match,PHstr

#===============================================================================================
# Read AEOLUS (for collocation)
#
#	INPUTS:
#		yyyymmddhh .......................... Current date in yyyymmddhh format
#		dateB4 .............................. Day before "yyyymmddhh" in yyyymmdd format
#		driver_dset_type .................... Aeolus dataset type: orig=original (not reprocessed), for reprocessed files use baseline number (example: B10)
#		driver_wind_type .................... Aeolus wind data type: RayClear=Rayleigh-clear, MieCloud=Mie-cloudy
#		bool_drv_qc ......................... Choice to apply Aeolus QC: True=apply QC, False=don't apply QC
#
#	OUTPUTS:
#		d_lat ............................... Latitude in degrees [-90,90]
#		d_lon ............................... Longitude in degrees [0,360]
#		d_prs ............................... Pressure in hPa
#		d_hgt ............................... Height in km
#		d_yr ................................ Year
#		d_mm ................................ Month
#		d_dy ................................ Day
#		d_hr ................................ Hour
#		d_mn ................................ Minute
#		indexesD ............................ Indices of input obs
#		qc_list ............................. List of QC applied (if applicable)
#
def read_aeolus(path_prefix,yyyymmddhh,dateB4,dateA,driver_dset_type,driver_wind_type,bool_drv_qc,dsetflag,runtype,time_diff_max,idxs):

  qc_list = ""			#initialize

  yyyy = yyyymmddhh[0:4]
  mm   = yyyymmddhh[4:6]
  dd   = yyyymmddhh[6:8]
  hour = yyyymmddhh[8:10]
  
  yyyymmdd = yyyy+mm+dd
  
  yyB4 = dateB4[0:4]
  mmB4 = dateB4[4:6]
  ddB4 = dateB4[6:8]
  
  yyA  = dateA[0:4]
  mmA  = dateA[4:6]
  ddA  = dateA[6:8]
  
	#-------------------------------------
        # Find hour limits ... for Aeolus datafiles

		# convert 'time_max_diff' to rounded integer
  tqc_time = float(time_diff_max)
  if tqc_time%60 != 0.0:
    qc_time = int(math.ceil(tqc_time/60.0))             # round up to nearest or equal integer hour
  else:
    qc_time = int(tqc_time)

		# convert hour to integer
  ihr = int(hour)

		# find integer hours for before/after files for collocation
  if hour == "00":
    	# 00
    ihrB = 24 - qc_time
    ihrA = ihr + qc_time
  else:
	# 06,12,18
    ihrB = ihr - qc_time
    ihrA = ihr + qc_time

  del qc_time,tqc_time,ihr

	#-------------------------------------------------
    	# Define dataset

  if driver_dset_type=='orig':						# Aeolus original dataset (not reprocessed)
    driver_dset_type_str = 'original'
  elif driver_dset_type!='orig':					# Aeolus dataset reprocessed by ESA
    driver_dset_type_str = 'reprocessed/2'+str(driver_dset_type)	# 	Reprocessed with Baseline B10 processor

  tmp_driver_path    = '/scratch/aeolus-dataset/netcdf/'+driver_dset_type_str+'/'+yyyy+'/'+mm+'/'
  tmp_driver_path_B4 = '/scratch/aeolus-dataset/netcdf/'+driver_dset_type_str+'/'+yyB4+'/'+mmB4+'/'
  tmp_driver_path_A  = '/scratch/aeolus-dataset/netcdf/'+driver_dset_type_str+'/'+yyA+'/'+mmA+'/'


		# Paths/files
		# ...Current date
  driver_path     	= path_prefix+tmp_driver_path
  driver_filename 	= 'Aeolus.L2B.'+driver_wind_type+'.1day.'+yyyymmdd+'.nc'
  driver_file     	= driver_path+driver_filename
		# ...Before date
  driver_path_B4  	= path_prefix+tmp_driver_path_B4
  driver_filename_B4 	= 'Aeolus.L2B.'+driver_wind_type+'.1day.'+dateB4+'.nc'
  driver_file_B4  	= driver_path_B4+driver_filename_B4
		# ...After date
  driver_path_A  	= path_prefix+tmp_driver_path_A
  driver_filename_A 	= 'Aeolus.L2B.'+driver_wind_type+'.1day.'+dateA+'.nc'
  driver_file_A  	= driver_path_A+driver_filename_A

  			# initialize flag indicating if dataset exists. 0=yes, 1=no
  existB = 0		# ... B = date Before current
  existA = 0		# ... A = date After current
  driver_exists = exists(driver_file)
  if driver_exists==False:
    print("ERROR: file "+driver_file+" does not exist!")
    sys.exit()
  driverB4_exists = exists(driver_file_B4)
  if driverB4_exists==False:
    print("WARNING: file "+driver_file_B4+" does not exist!")
    existB = 1
  driverA_exists = exists(driver_file_A)
  if driverA_exists==False:
    print("WARNING: file "+driver_file_A+" does not exist!")
    existA = 1

  print("read_aeolus: AEOLUS file = "+str(driver_file))

		# Paths on FTP/web archive server (for output NetCDF only)
  str_driver_path 	= driver_path
  str_driver_path_B4 	= driver_path_B4
  str_driver_path_A	= driver_path_A

  drv_src  = str_driver_path_B4+driver_filename_B4
  drv_src += ", "+str_driver_path+driver_filename
  drv_src += ", "+str_driver_path_A+driver_filename_A

		# Variable names
  drv_lat_var 		= 'latitude'
  drv_lon_var 		= 'longitude'
  drv_prs_var 		= 'pressure'
  drv_hgt_var 		= 'height_mid'
  drv_yr_var  		= 'year'
  drv_mm_var  		= 'month'
  drv_dy_var 		= 'day'
  drv_hr_var  		= 'hour'
  drv_mn_var  		= 'minute'

  drv_err_var 		= 'HLOS_error'
  drv_len_var 		= 'length'
  drv_hgtTop_var 	= 'height_top'
  drv_hgtBot_var 	= 'height_bot'
  drv_spd_var		= 'HLOS_wind_velocity'
  drv_dir_var		= 'HLOS_azimuth_angle'

    	#-------------------------------------------------
  	# Load datasets
  
  	#`````````````````````````````````````````````````
	# CURRENT DATE
  data_hdl = Dataset(driver_file)

  tdrvC_lat = np.asarray( data_hdl.variables[drv_lat_var] )
  tdrvC_lon = np.asarray( data_hdl.variables[drv_lon_var] )
  tdrvC_prs = np.asarray( data_hdl.variables[drv_prs_var] )
  tdrvC_hgt = np.asarray( data_hdl.variables[drv_hgt_var] )
  tdrvC_yr  = np.asarray( data_hdl.variables[drv_yr_var]  )
  tdrvC_mm  = np.asarray( data_hdl.variables[drv_mm_var]  )
  tdrvC_dy  = np.asarray( data_hdl.variables[drv_dy_var]  )
  tdrvC_hr  = np.asarray( data_hdl.variables[drv_hr_var]  )
  tdrvC_mn  = np.asarray( data_hdl.variables[drv_mn_var]  )

  tdrvC_err     = np.asarray( data_hdl.variables[drv_err_var]  )
  tdrvC_len     = np.asarray( data_hdl.variables[drv_len_var]  )
  tdrvC_hgtTop  = np.asarray( data_hdl.variables[drv_hgtTop_var]  )
  tdrvC_hgtBot  = np.asarray( data_hdl.variables[drv_hgtBot_var]  )
  tdrvC_spd     = np.asarray( data_hdl.variables[drv_spd_var]  )
  tdrvC_dir     = np.asarray( data_hdl.variables[drv_dir_var]  )

  data_hdl.close() 

	# check pressure units and convert to hPa
  if max(tdrvC_prs) > 10000.:
    tdrvC_prs = tdrvC_prs/100.

	# check height units and convert to km
  if max(tdrvC_hgt) > 1000.:
    tdrvC_hgt = tdrvC_hgt/1000.

	# check height Top units and convert to km
  if max(tdrvC_hgtTop) > 1000.:
    tdrvC_hgtTop = tdrvC_hgtTop/1000.

	# check height Bottom units and convert to km
  if max(tdrvC_hgtBot) > 1000.:
    tdrvC_hgtBot = tdrvC_hgtBot/1000.

  if existB == 0:
  	#`````````````````````````````````````````````````
	# DATE BEFORE
    data_hdl = Dataset(driver_file_B4)

    tdrvB4_lat = np.asarray( data_hdl.variables[drv_lat_var] )
    tdrvB4_lon = np.asarray( data_hdl.variables[drv_lon_var] )
    tdrvB4_prs = np.asarray( data_hdl.variables[drv_prs_var] )
    tdrvB4_hgt = np.asarray( data_hdl.variables[drv_hgt_var] )
    tdrvB4_yr  = np.asarray( data_hdl.variables[drv_yr_var]  )
    tdrvB4_mm  = np.asarray( data_hdl.variables[drv_mm_var]  )
    tdrvB4_dy  = np.asarray( data_hdl.variables[drv_dy_var]  )
    tdrvB4_hr  = np.asarray( data_hdl.variables[drv_hr_var]  )
    tdrvB4_mn  = np.asarray( data_hdl.variables[drv_mn_var]  )

    tdrvB4_err  = np.asarray( data_hdl.variables[drv_err_var]  )
    tdrvB4_len  = np.asarray( data_hdl.variables[drv_len_var]  )
    tdrvB4_hgtTop = np.asarray( data_hdl.variables[drv_hgtTop_var]  )
    tdrvB4_hgtBot = np.asarray( data_hdl.variables[drv_hgtBot_var]  )
    tdrvB4_spd  = np.asarray( data_hdl.variables[drv_spd_var]  )
    tdrvB4_dir  = np.asarray( data_hdl.variables[drv_dir_var]  )

    data_hdl.close()

      # check pressure units and convert to hPa
    if max(tdrvB4_prs) > 10000.:
      tdrvB4_prs = tdrvB4_prs/100.

      # check height units and convert to km
    if max(tdrvB4_hgt) > 1000.:
      tdrvB4_hgt = tdrvB4_hgt/1000.

      # check height Top units and convert to km
    if max(tdrvB4_hgtTop) > 1000.:
      tdrvB4_hgtTop = tdrvB4_hgtTop/1000.
  
      # check height Bottom units and convert to km
    if max(tdrvB4_hgtBot) > 1000.:
      tdrvB4_hgtBot = tdrvB4_hgtBot/1000.
    
  if existA == 0:
    	#`````````````````````````````````````````````````
	# DATE AFTER
    data_hdl = Dataset(driver_file_A)

    tdrvA_lat = np.asarray( data_hdl.variables[drv_lat_var] )
    tdrvA_lon = np.asarray( data_hdl.variables[drv_lon_var] )
    tdrvA_prs = np.asarray( data_hdl.variables[drv_prs_var] )
    tdrvA_hgt = np.asarray( data_hdl.variables[drv_hgt_var] )
    tdrvA_yr  = np.asarray( data_hdl.variables[drv_yr_var]  )
    tdrvA_mm  = np.asarray( data_hdl.variables[drv_mm_var]  )
    tdrvA_dy  = np.asarray( data_hdl.variables[drv_dy_var]  )
    tdrvA_hr  = np.asarray( data_hdl.variables[drv_hr_var]  )
    tdrvA_mn  = np.asarray( data_hdl.variables[drv_mn_var]  )

    tdrvA_err  = np.asarray( data_hdl.variables[drv_err_var]  )
    tdrvA_len  = np.asarray( data_hdl.variables[drv_len_var]  )
    tdrvA_hgtTop = np.asarray( data_hdl.variables[drv_hgtTop_var]  )
    tdrvA_hgtBot = np.asarray( data_hdl.variables[drv_hgtBot_var]  )
    tdrvA_spd  = np.asarray( data_hdl.variables[drv_spd_var]  )
    tdrvA_dir  = np.asarray( data_hdl.variables[drv_dir_var]  )

    data_hdl.close()

      # check pressure units and convert to hPa
    if max(tdrvA_prs) > 10000.:
      tdrvA_prs = tdrvA_prs/100.

      # check height units and convert to km
    if max(tdrvA_hgt) > 1000.:
      tdrvA_hgt = tdrvA_hgt/1000.

      # check height Top units and convert to km
    if max(tdrvA_hgtTop) > 1000.:
      tdrvA_hgtTop = tdrvA_hgtTop/1000.
  
      # check height Bottom units and convert to km
    if max(tdrvA_hgtBot) > 1000.:
      tdrvA_hgtBot = tdrvA_hgtBot/1000.

	#-------------------------------------------------
      	# Extract 6-hr range from Aeolus day arrays: +/- 3 hours around center analysis hour of current date (00, 06, 12, or 18 UTC)

  i_tdrvC_hr  = tdrvC_hr.astype("int")
  if existB == 0:
    i_tdrvB4_hr = tdrvB4_hr.astype("int")
  if existA == 0:
    i_tdrvA_hr  = tdrvA_hr.astype("int")
  
  if hour == "00":
    tiHHC1 = np.where(((i_tdrvC_hr >= 0) * (i_tdrvC_hr < 3)))
    siHHC1 = np.asarray(tiHHC1)
    iHHC1  = siHHC1.flatten()

    if existB == 0:
      tiHHC2 = np.where(((i_tdrvB4_hr >= 21) * (i_tdrvB4_hr < 24)))  
      siHHC2 = np.asarray(tiHHC2)
      iHHC2  = siHHC2.flatten()

    hdrvC1_lat     = tdrvC_lat[iHHC1]
    hdrvC1_lon     = tdrvC_lon[iHHC1]
    hdrvC1_prs     = tdrvC_prs[iHHC1]
    hdrvC1_hgt     = tdrvC_hgt[iHHC1]
    hdrvC1_yr      = tdrvC_yr [iHHC1]
    hdrvC1_mm      = tdrvC_mm [iHHC1]
    hdrvC1_dy      = tdrvC_dy [iHHC1]
    hdrvC1_hr      = tdrvC_hr [iHHC1]
    hdrvC1_mn      = tdrvC_mn [iHHC1]
    hdrvC1_err	   = tdrvC_err[iHHC1]
    hdrvC1_len	   = tdrvC_len[iHHC1]
    hdrvC1_hgtTop  = tdrvC_hgtTop[iHHC1]
    hdrvC1_hgtBot  = tdrvC_hgtBot[iHHC1]
    hdrvC1_spd     = tdrvC_spd[iHHC1]
    hdrvC1_dir     = tdrvC_dir[iHHC1]

    if existB == 0:
      hdrvC2_lat	= tdrvB4_lat[iHHC2]
      hdrvC2_lon	= tdrvB4_lon[iHHC2]
      hdrvC2_prs	= tdrvB4_prs[iHHC2]
      hdrvC2_hgt	= tdrvB4_hgt[iHHC2]
      hdrvC2_yr		= tdrvB4_yr [iHHC2]
      hdrvC2_mm		= tdrvB4_mm [iHHC2]
      hdrvC2_dy		= tdrvB4_dy [iHHC2]
      hdrvC2_hr		= tdrvB4_hr [iHHC2]
      hdrvC2_mn		= tdrvB4_mn [iHHC2]
      hdrvC2_err	= tdrvB4_err[iHHC2]
      hdrvC2_len	= tdrvB4_len[iHHC2]
      hdrvC2_hgtTop 	= tdrvB4_hgtTop[iHHC2]
      hdrvC2_hgtBot 	= tdrvB4_hgtBot[iHHC2]
      hdrvC2_spd    	= tdrvB4_spd[iHHC2]
      hdrvC2_dir    	= tdrvB4_dir[iHHC2]
    
    	# Append current arrays to before date arrays
    if existB == 0:
      hdrvC_lat    = np.append(hdrvC2_lat,hdrvC1_lat,axis=0)
      hdrvC_lon    = np.append(hdrvC2_lon,hdrvC1_lon,axis=0)
      hdrvC_prs    = np.append(hdrvC2_prs,hdrvC1_prs,axis=0)
      hdrvC_hgt    = np.append(hdrvC2_hgt,hdrvC1_hgt,axis=0)
      hdrvC_yr     = np.append(hdrvC2_yr ,hdrvC1_yr ,axis=0)
      hdrvC_mm     = np.append(hdrvC2_mm ,hdrvC1_mm ,axis=0)
      hdrvC_dy     = np.append(hdrvC2_dy ,hdrvC1_dy ,axis=0)
      hdrvC_hr     = np.append(hdrvC2_hr ,hdrvC1_hr ,axis=0)
      hdrvC_mn     = np.append(hdrvC2_mn ,hdrvC1_mn ,axis=0)
      hdrvC_err    = np.append(hdrvC2_err,hdrvC1_err,axis=0)
      hdrvC_len    = np.append(hdrvC2_len,hdrvC1_len,axis=0)
      hdrvC_hgtTop = np.append(hdrvC2_hgtTop,hdrvC1_hgtTop,axis=0)
      hdrvC_hgtBot = np.append(hdrvC2_hgtBot,hdrvC1_hgtBot,axis=0)
      hdrvC_spd    = np.append(hdrvC2_spd,hdrvC1_spd,axis=0)
      hdrvC_dir    = np.append(hdrvC2_dir,hdrvC1_dir,axis=0)
    else:
      hdrvC_lat    = hdrvC1_lat
      hdrvC_lon    = hdrvC1_lon
      hdrvC_prs    = hdrvC1_prs
      hdrvC_hgt    = hdrvC1_hgt
      hdrvC_yr     = hdrvC1_yr 
      hdrvC_mm     = hdrvC1_mm 
      hdrvC_dy     = hdrvC1_dy 
      hdrvC_hr     = hdrvC1_hr 
      hdrvC_mn     = hdrvC1_mn 
      hdrvC_err    = hdrvC1_err
      hdrvC_len    = hdrvC1_len
      hdrvC_hgtTop = hdrvC1_hgtTop
      hdrvC_hgtBot = hdrvC1_hgtBot
      hdrvC_spd    = hdrvC1_spd
      hdrvC_dir    = hdrvC1_dir
    
    if dsetflag == "dep":
    	# keep ihrB hrs before CURRENT and ihrA hrs after CURRENT 6hrs
      tiHHA  = np.where(((i_tdrvC_hr >= 3) * (i_tdrvC_hr < (3+ihrA))))
      siHHA  = np.asarray(tiHHA)
      iHHA   = siHHA.flatten()
      
      if existB == 0:
        tiHHB4 = np.where(((i_tdrvB4_hr >= (21-ihrB)) * (i_tdrvB4_hr < 21)))
        siHHB4 = np.asarray(tiHHB4)
        iHHB4  = siHHB4.flatten()

      hdrvA_lat     = tdrvC_lat[iHHA]
      hdrvA_lon     = tdrvC_lon[iHHA]
      hdrvA_prs     = tdrvC_prs[iHHA]
      hdrvA_hgt     = tdrvC_hgt[iHHA]
      hdrvA_yr	    = tdrvC_yr [iHHA]
      hdrvA_mm	    = tdrvC_mm [iHHA]
      hdrvA_dy	    = tdrvC_dy [iHHA]
      hdrvA_hr	    = tdrvC_hr [iHHA]
      hdrvA_mn	    = tdrvC_mn [iHHA]
      hdrvA_err     = tdrvC_err[iHHA]
      hdrvA_len     = tdrvC_len[iHHA]
      hdrvA_hgtTop  = tdrvC_hgtTop[iHHA]
      hdrvA_hgtBot  = tdrvC_hgtBot[iHHA]
      hdrvA_spd     = tdrvC_spd[iHHA]
      hdrvA_dir     = tdrvC_dir[iHHA]

      if existB == 0:
        hdrvB4_lat     = tdrvB4_lat[iHHB4]
        hdrvB4_lon     = tdrvB4_lon[iHHB4]
        hdrvB4_prs     = tdrvB4_prs[iHHB4]
        hdrvB4_hgt     = tdrvB4_hgt[iHHB4]
        hdrvB4_yr      = tdrvB4_yr [iHHB4]
        hdrvB4_mm      = tdrvB4_mm [iHHB4]
        hdrvB4_dy      = tdrvB4_dy [iHHB4]
        hdrvB4_hr      = tdrvB4_hr [iHHB4]
        hdrvB4_mn      = tdrvB4_mn [iHHB4]
        hdrvB4_err     = tdrvB4_err[iHHB4]
        hdrvB4_len     = tdrvB4_len[iHHB4]
        hdrvB4_hgtTop  = tdrvB4_hgtTop[iHHB4]
        hdrvB4_hgtBot  = tdrvB4_hgtBot[iHHB4]
        hdrvB4_spd     = tdrvB4_spd[iHHB4]
        hdrvB4_dir     = tdrvB4_dir[iHHB4]

		# Append current arrays to before date arrays
      if existB == 0:
        drv_lat    = np.append(hdrvB4_lat ,hdrvC_lat,axis=0)
        drv_lon    = np.append(hdrvB4_lon ,hdrvC_lon,axis=0)
        drv_prs    = np.append(hdrvB4_prs ,hdrvC_prs,axis=0)
        drv_hgt    = np.append(hdrvB4_hgt ,hdrvC_hgt,axis=0)
        drv_yr     = np.append(hdrvB4_yr  ,hdrvC_yr ,axis=0)
        drv_mm     = np.append(hdrvB4_mm  ,hdrvC_mm ,axis=0)
        drv_dy     = np.append(hdrvB4_dy  ,hdrvC_dy ,axis=0)
        drv_hr     = np.append(hdrvB4_hr  ,hdrvC_hr ,axis=0)
        drv_mn     = np.append(hdrvB4_mn  ,hdrvC_mn ,axis=0)
        drv_err	   = np.append(hdrvB4_err ,hdrvC_err,axis=0)
        drv_len	   = np.append(hdrvB4_len ,hdrvC_len,axis=0)
        drv_hgtTop = np.append(hdrvB4_hgtTop,hdrvC_hgtTop,axis=0)
        drv_hgtBot = np.append(hdrvB4_hgtBot,hdrvC_hgtBot,axis=0)
        drv_spd    = np.append(hdrvB4_spd ,hdrvC_spd,axis=0)
        drv_dir    = np.append(hdrvB4_dir ,hdrvC_dir,axis=0)
      else:
        drv_lat    = hdrvC_lat
        drv_lon    = hdrvC_lon
        drv_prs    = hdrvC_prs
        drv_hgt    = hdrvC_hgt
        drv_yr     = hdrvC_yr 
        drv_mm     = hdrvC_mm 
        drv_dy     = hdrvC_dy 
        drv_hr     = hdrvC_hr 
        drv_mn     = hdrvC_mn 
        drv_err	   = hdrvC_err
        drv_len	   = hdrvC_len
        drv_hgtTop = hdrvC_hgtTop
        drv_hgtBot = hdrvC_hgtBot
        drv_spd    = hdrvC_spd
        drv_dir    = hdrvC_dir
      
      drv_lat    = np.append(drv_lat ,hdrvA_lat,axis=0)
      drv_lon    = np.append(drv_lon ,hdrvA_lon,axis=0)
      drv_prs    = np.append(drv_prs ,hdrvA_prs,axis=0)
      drv_hgt    = np.append(drv_hgt ,hdrvA_hgt,axis=0)
      drv_yr     = np.append(drv_yr  ,hdrvA_yr ,axis=0)
      drv_mm     = np.append(drv_mm  ,hdrvA_mm ,axis=0)
      drv_dy     = np.append(drv_dy  ,hdrvA_dy ,axis=0)
      drv_hr     = np.append(drv_hr  ,hdrvA_hr ,axis=0)
      drv_mn     = np.append(drv_mn  ,hdrvA_mn ,axis=0)
      drv_err	 = np.append(drv_err ,hdrvA_err   ,axis=0)
      drv_len	 = np.append(drv_len ,hdrvA_len   ,axis=0)
      drv_hgtTop = np.append(drv_hgtTop,hdrvA_hgtTop,axis=0)
      drv_hgtBot = np.append(drv_hgtBot,hdrvA_hgtBot,axis=0)
      drv_spd    = np.append(drv_spd ,hdrvA_spd,axis=0)
      drv_dir    = np.append(drv_dir ,hdrvA_dir,axis=0)
      
    else:
    	# dsetflag = "drv"
      drv_lat    = hdrvC_lat
      drv_lon    = hdrvC_lon
      drv_prs    = hdrvC_prs
      drv_hgt    = hdrvC_hgt
      drv_yr     = hdrvC_yr 
      drv_mm     = hdrvC_mm 
      drv_dy     = hdrvC_dy 
      drv_hr     = hdrvC_hr 
      drv_mn     = hdrvC_mn 
      drv_err	 = hdrvC_err   
      drv_len	 = hdrvC_len   
      drv_hgtTop = hdrvC_hgtTop
      drv_hgtBot = hdrvC_hgtBot
      drv_spd    = hdrvC_spd
      drv_dir    = hdrvC_dir
      
  elif hour == "06":
    tiHHC1 = np.where(((i_tdrvC_hr >= 3) * (i_tdrvC_hr < 9)))
    siHHC1 = np.asarray(tiHHC1)
    iHHC1  = siHHC1.flatten()

    hdrvC_lat	  = tdrvC_lat[iHHC1]
    hdrvC_lon	  = tdrvC_lon[iHHC1]
    hdrvC_prs	  = tdrvC_prs[iHHC1]
    hdrvC_hgt	  = tdrvC_hgt[iHHC1]
    hdrvC_yr	  = tdrvC_yr [iHHC1]
    hdrvC_mm	  = tdrvC_mm [iHHC1]
    hdrvC_dy	  = tdrvC_dy [iHHC1]
    hdrvC_hr	  = tdrvC_hr [iHHC1]
    hdrvC_mn	  = tdrvC_mn [iHHC1]
    hdrvC_err	  = tdrvC_err[iHHC1]
    hdrvC_len	  = tdrvC_len[iHHC1]
    hdrvC_hgtTop  = tdrvC_hgtTop[iHHC1]
    hdrvC_hgtBot  = tdrvC_hgtBot[iHHC1]
    hdrvC_spd     = tdrvC_spd[iHHC1]
    hdrvC_dir     = tdrvC_dir[iHHC1]

    if dsetflag == "dep":
    	# keep ihrB before CURRENT and ihrA after CURRENT 6hrs
      tiHHA   = np.where(((i_tdrvC_hr >= 9) * (i_tdrvC_hr < (9+ihrA))))
      siHHA   = np.asarray(tiHHA)
      iHHA    = siHHA.flatten()
      
      if (3-ihrB) < 0:
        dhrB = abs(3 - ihrB)
      
        tiHHB42 = np.where(((i_tdrvC_hr >= 0) * (i_tdrvC_hr < 3)))
        siHHB42 = np.asarray(tiHHB42)
        iHHB42  = siHHB42.flatten()
      
        if existB == 0:
          tiHHB41 = np.where(((i_tdrvB4_hr >= (24-dhrB)) * (i_tdrvB4_hr < 24)))  
          siHHB41 = np.asarray(tiHHB41)
          iHHB41  = siHHB41.flatten()

        del dhrB
      else:
        tiHHB42 = np.where(((i_tdrvC_hr >= (3-ihrB)) * (i_tdrvC_hr < 3)))
        siHHB42 = np.asarray(tiHHB42)
        iHHB42  = siHHB42.flatten()

      hdrvA_lat     = tdrvC_lat[iHHA]
      hdrvA_lon     = tdrvC_lon[iHHA]
      hdrvA_prs     = tdrvC_prs[iHHA]
      hdrvA_hgt     = tdrvC_hgt[iHHA]
      hdrvA_yr	    = tdrvC_yr [iHHA]
      hdrvA_mm	    = tdrvC_mm [iHHA]
      hdrvA_dy	    = tdrvC_dy [iHHA]
      hdrvA_hr	    = tdrvC_hr [iHHA]
      hdrvA_mn	    = tdrvC_mn [iHHA]
      hdrvA_err     = tdrvC_err[iHHA]
      hdrvA_len     = tdrvC_len[iHHA]
      hdrvA_hgtTop  = tdrvC_hgtTop[iHHA]
      hdrvA_hgtBot  = tdrvC_hgtBot[iHHA]
      hdrvA_spd     = tdrvC_spd[iHHA]
      hdrvA_dir     = tdrvC_dir[iHHA]

      if existB == 0:
        hdrvB41_lat     = tdrvB4_lat[iHHB41]
        hdrvB41_lon     = tdrvB4_lon[iHHB41]
        hdrvB41_prs     = tdrvB4_prs[iHHB41]
        hdrvB41_hgt     = tdrvB4_hgt[iHHB41]
        hdrvB41_yr      = tdrvB4_yr [iHHB41]
        hdrvB41_mm      = tdrvB4_mm [iHHB41]
        hdrvB41_dy      = tdrvB4_dy [iHHB41]
        hdrvB41_hr      = tdrvB4_hr [iHHB41]
        hdrvB41_mn      = tdrvB4_mn [iHHB41]
        hdrvB41_err     = tdrvB4_err[iHHB41]
        hdrvB41_len     = tdrvB4_len[iHHB41]
        hdrvB41_hgtTop  = tdrvB4_hgtTop[iHHB41]
        hdrvB41_hgtBot  = tdrvB4_hgtBot[iHHB41]
        hdrvB41_spd     = tdrvB4_spd[iHHB41]
        hdrvB41_dir     = tdrvB4_dir[iHHB41]
      
      hdrvB42_lat     = tdrvC_lat[iHHB42]
      hdrvB42_lon     = tdrvC_lon[iHHB42]
      hdrvB42_prs     = tdrvC_prs[iHHB42]
      hdrvB42_hgt     = tdrvC_hgt[iHHB42]
      hdrvB42_yr      = tdrvC_yr [iHHB42]
      hdrvB42_mm      = tdrvC_mm [iHHB42]
      hdrvB42_dy      = tdrvC_dy [iHHB42]
      hdrvB42_hr      = tdrvC_hr [iHHB42]
      hdrvB42_mn      = tdrvC_mn [iHHB42]
      hdrvB42_err     = tdrvC_err[iHHB42]
      hdrvB42_len     = tdrvC_len[iHHB42]
      hdrvB42_hgtTop  = tdrvC_hgtTop[iHHB42]
      hdrvB42_hgtBot  = tdrvC_hgtBot[iHHB42]
      hdrvB42_spd     = tdrvC_spd[iHHB42]
      hdrvB42_dir     = tdrvC_dir[iHHB42]

		# Append current arrays to before date arrays
      if existB == 0:
        drv_lat    = np.append(hdrvB41_lat ,hdrvB42_lat,axis=0)
        drv_lon    = np.append(hdrvB41_lon ,hdrvB42_lon,axis=0)
        drv_prs    = np.append(hdrvB41_prs ,hdrvB42_prs,axis=0)
        drv_hgt    = np.append(hdrvB41_hgt ,hdrvB42_hgt,axis=0)
        drv_yr     = np.append(hdrvB41_yr  ,hdrvB42_yr ,axis=0)
        drv_mm     = np.append(hdrvB41_mm  ,hdrvB42_mm ,axis=0)
        drv_dy     = np.append(hdrvB41_dy  ,hdrvB42_dy ,axis=0)
        drv_hr     = np.append(hdrvB41_hr  ,hdrvB42_hr ,axis=0)
        drv_mn     = np.append(hdrvB41_mn  ,hdrvB42_mn ,axis=0)
        drv_err	   = np.append(hdrvB41_err ,hdrvB42_err,axis=0)
        drv_len	   = np.append(hdrvB41_len ,hdrvB42_len,axis=0)
        drv_hgtTop = np.append(hdrvB41_hgtTop,hdrvB42_hgtTop,axis=0)
        drv_hgtBot = np.append(hdrvB41_hgtBot,hdrvB42_hgtBot,axis=0)
        drv_spd    = np.append(hdrvB41_spd ,hdrvB42_spd,axis=0)
        drv_dir    = np.append(hdrvB41_dir ,hdrvB42_dir,axis=0)
      else:
        drv_lat    = hdrvB42_lat
        drv_lon    = hdrvB42_lon
        drv_prs    = hdrvB42_prs
        drv_hgt    = hdrvB42_hgt
        drv_yr     = hdrvB42_yr 
        drv_mm     = hdrvB42_mm 
        drv_dy     = hdrvB42_dy 
        drv_hr     = hdrvB42_hr 
        drv_mn     = hdrvB42_mn 
        drv_err	   = hdrvB42_err
        drv_len	   = hdrvB42_len
        drv_hgtTop = hdrvB42_hgtTop
        drv_hgtBot = hdrvB42_hgtBot
        drv_spd    = hdrvB42_spd
        drv_dir    = hdrvB42_dir
      
      drv_lat    = np.append(drv_lat ,hdrvC_lat,axis=0)
      drv_lon    = np.append(drv_lon ,hdrvC_lon,axis=0)
      drv_prs    = np.append(drv_prs ,hdrvC_prs,axis=0)
      drv_hgt    = np.append(drv_hgt ,hdrvC_hgt,axis=0)
      drv_yr     = np.append(drv_yr  ,hdrvC_yr ,axis=0)
      drv_mm     = np.append(drv_mm  ,hdrvC_mm ,axis=0)
      drv_dy     = np.append(drv_dy  ,hdrvC_dy ,axis=0)
      drv_hr     = np.append(drv_hr  ,hdrvC_hr ,axis=0)
      drv_mn     = np.append(drv_mn  ,hdrvC_mn ,axis=0)
      drv_err	 = np.append(drv_err ,hdrvC_err   ,axis=0)
      drv_len	 = np.append(drv_len ,hdrvC_len   ,axis=0)
      drv_hgtTop = np.append(drv_hgtTop,hdrvC_hgtTop,axis=0)
      drv_hgtBot = np.append(drv_hgtBot,hdrvC_hgtBot,axis=0)
      drv_spd    = np.append(drv_spd ,hdrvC_spd,axis=0)
      drv_dir    = np.append(drv_dir ,hdrvC_dir,axis=0)
      
      drv_lat    = np.append(drv_lat ,hdrvA_lat,axis=0)
      drv_lon    = np.append(drv_lon ,hdrvA_lon,axis=0)
      drv_prs    = np.append(drv_prs ,hdrvA_prs,axis=0)
      drv_hgt    = np.append(drv_hgt ,hdrvA_hgt,axis=0)
      drv_yr     = np.append(drv_yr  ,hdrvA_yr ,axis=0)
      drv_mm     = np.append(drv_mm  ,hdrvA_mm ,axis=0)
      drv_dy     = np.append(drv_dy  ,hdrvA_dy ,axis=0)
      drv_hr     = np.append(drv_hr  ,hdrvA_hr ,axis=0)
      drv_mn     = np.append(drv_mn  ,hdrvA_mn ,axis=0)
      drv_err	 = np.append(drv_err ,hdrvA_err   ,axis=0)
      drv_len	 = np.append(drv_len ,hdrvA_len   ,axis=0)
      drv_hgtTop = np.append(drv_hgtTop,hdrvA_hgtTop,axis=0)
      drv_hgtBot = np.append(drv_hgtBot,hdrvA_hgtBot,axis=0)
      drv_spd    = np.append(drv_spd ,hdrvA_spd,axis=0)
      drv_dir    = np.append(drv_dir ,hdrvA_dir,axis=0)
      
    else:
    	# dsetflag = "drv"
      drv_lat    = hdrvC_lat
      drv_lon    = hdrvC_lon
      drv_prs    = hdrvC_prs
      drv_hgt    = hdrvC_hgt
      drv_yr     = hdrvC_yr 
      drv_mm     = hdrvC_mm 
      drv_dy     = hdrvC_dy 
      drv_hr     = hdrvC_hr 
      drv_mn     = hdrvC_mn 
      drv_err	 = hdrvC_err   
      drv_len	 = hdrvC_len   
      drv_hgtTop = hdrvC_hgtTop
      drv_hgtBot = hdrvC_hgtBot  
      drv_spd    = hdrvC_spd
      drv_dir    = hdrvC_dir   
      
  elif hour == "12":
    tiHHC1 = np.where(((i_tdrvC_hr >= 9) * (i_tdrvC_hr < 15)))
    siHHC1 = np.asarray(tiHHC1)
    iHHC1  = siHHC1.flatten()

    hdrvC_lat	  = tdrvC_lat[iHHC1]
    hdrvC_lon	  = tdrvC_lon[iHHC1]
    hdrvC_prs	  = tdrvC_prs[iHHC1]
    hdrvC_hgt	  = tdrvC_hgt[iHHC1]
    hdrvC_yr	  = tdrvC_yr [iHHC1]
    hdrvC_mm	  = tdrvC_mm [iHHC1]
    hdrvC_dy	  = tdrvC_dy [iHHC1]
    hdrvC_hr	  = tdrvC_hr [iHHC1]
    hdrvC_mn	  = tdrvC_mn [iHHC1]
    hdrvC_err	  = tdrvC_err[iHHC1]
    hdrvC_len	  = tdrvC_len[iHHC1]
    hdrvC_hgtTop  = tdrvC_hgtTop[iHHC1]
    hdrvC_hgtBot  = tdrvC_hgtBot[iHHC1]
    hdrvC_spd     = tdrvC_spd[iHHC1]
    hdrvC_dir     = tdrvC_dir[iHHC1]

    if dsetflag == "dep":
    	# keep ihrB before CURRENT and ihrA after CURRENT 6hrs
      tiHHA  = np.where(((i_tdrvC_hr >= 15) * (i_tdrvC_hr < (15+ihrA)))) 
      tiHHB4 = np.where(((i_tdrvC_hr >= (9-ihrB)) * (i_tdrvC_hr < 9)))
    
      siHHA  = np.asarray(tiHHA)
      iHHA   = siHHA.flatten()
      siHHB4 = np.asarray(tiHHB4)
      iHHB4  = siHHB4.flatten()

      hdrvA_lat     = tdrvC_lat[iHHA]
      hdrvA_lon     = tdrvC_lon[iHHA]
      hdrvA_prs     = tdrvC_prs[iHHA]
      hdrvA_hgt     = tdrvC_hgt[iHHA]
      hdrvA_yr	    = tdrvC_yr [iHHA]
      hdrvA_mm	    = tdrvC_mm [iHHA]
      hdrvA_dy	    = tdrvC_dy [iHHA]
      hdrvA_hr	    = tdrvC_hr [iHHA]
      hdrvA_mn	    = tdrvC_mn [iHHA]
      hdrvA_err     = tdrvC_err[iHHA]
      hdrvA_len     = tdrvC_len[iHHA]
      hdrvA_hgtTop  = tdrvC_hgtTop[iHHA]
      hdrvA_hgtBot  = tdrvC_hgtBot[iHHA]
      hdrvA_spd     = tdrvC_spd[iHHA]
      hdrvA_dir     = tdrvC_dir[iHHA]

      hdrvB4_lat     = tdrvC_lat[iHHB4]
      hdrvB4_lon     = tdrvC_lon[iHHB4]
      hdrvB4_prs     = tdrvC_prs[iHHB4]
      hdrvB4_hgt     = tdrvC_hgt[iHHB4]
      hdrvB4_yr      = tdrvC_yr [iHHB4]
      hdrvB4_mm      = tdrvC_mm [iHHB4]
      hdrvB4_dy      = tdrvC_dy [iHHB4]
      hdrvB4_hr      = tdrvC_hr [iHHB4]
      hdrvB4_mn      = tdrvC_mn [iHHB4]
      hdrvB4_err     = tdrvC_err[iHHB4]
      hdrvB4_len     = tdrvC_len[iHHB4]
      hdrvB4_hgtTop  = tdrvC_hgtTop[iHHB4]
      hdrvB4_hgtBot  = tdrvC_hgtBot[iHHB4]
      hdrvB4_spd     = tdrvC_spd[iHHB4]
      hdrvB4_dir     = tdrvC_dir[iHHB4]

		# Append current arrays to before date arrays
      drv_lat    = np.append(hdrvB4_lat ,hdrvC_lat,axis=0)
      drv_lon    = np.append(hdrvB4_lon ,hdrvC_lon,axis=0)
      drv_prs    = np.append(hdrvB4_prs ,hdrvC_prs,axis=0)
      drv_hgt    = np.append(hdrvB4_hgt ,hdrvC_hgt,axis=0)
      drv_yr     = np.append(hdrvB4_yr  ,hdrvC_yr ,axis=0)
      drv_mm     = np.append(hdrvB4_mm  ,hdrvC_mm ,axis=0)
      drv_dy     = np.append(hdrvB4_dy  ,hdrvC_dy ,axis=0)
      drv_hr     = np.append(hdrvB4_hr  ,hdrvC_hr ,axis=0)
      drv_mn     = np.append(hdrvB4_mn  ,hdrvC_mn ,axis=0)
      drv_err	 = np.append(hdrvB4_err ,hdrvC_err   ,axis=0)
      drv_len	 = np.append(hdrvB4_len ,hdrvC_len   ,axis=0)
      drv_hgtTop = np.append(hdrvB4_hgtTop,hdrvC_hgtTop,axis=0)
      drv_hgtBot = np.append(hdrvB4_hgtBot,hdrvC_hgtBot,axis=0)
      drv_spd    = np.append(hdrvB4_spd ,hdrvC_spd,axis=0)
      drv_dir    = np.append(hdrvB4_dir ,hdrvC_dir,axis=0)
      
      drv_lat    = np.append(drv_lat ,hdrvA_lat,axis=0)
      drv_lon    = np.append(drv_lon ,hdrvA_lon,axis=0)
      drv_prs    = np.append(drv_prs ,hdrvA_prs,axis=0)
      drv_hgt    = np.append(drv_hgt ,hdrvA_hgt,axis=0)
      drv_yr     = np.append(drv_yr  ,hdrvA_yr ,axis=0)
      drv_mm     = np.append(drv_mm  ,hdrvA_mm ,axis=0)
      drv_dy     = np.append(drv_dy  ,hdrvA_dy ,axis=0)
      drv_hr     = np.append(drv_hr  ,hdrvA_hr ,axis=0)
      drv_mn     = np.append(drv_mn  ,hdrvA_mn ,axis=0)
      drv_err	 = np.append(drv_err ,hdrvA_err   ,axis=0)
      drv_len	 = np.append(drv_len ,hdrvA_len   ,axis=0)
      drv_hgtTop = np.append(drv_hgtTop,hdrvA_hgtTop,axis=0)
      drv_hgtBot = np.append(drv_hgtBot,hdrvA_hgtBot,axis=0)
      drv_spd    = np.append(drv_spd ,hdrvA_spd,axis=0)
      drv_dir    = np.append(drv_dir ,hdrvA_dir,axis=0)
      
    else:
    	# dsetflag = "drv"
      drv_lat    = hdrvC_lat
      drv_lon    = hdrvC_lon
      drv_prs    = hdrvC_prs
      drv_hgt    = hdrvC_hgt
      drv_yr     = hdrvC_yr 
      drv_mm     = hdrvC_mm 
      drv_dy     = hdrvC_dy 
      drv_hr     = hdrvC_hr 
      drv_mn     = hdrvC_mn 
      drv_err	 = hdrvC_err   
      drv_len	 = hdrvC_len   
      drv_hgtTop = hdrvC_hgtTop
      drv_hgtBot = hdrvC_hgtBot 
      drv_spd    = hdrvC_spd
      drv_dir    = hdrvC_dir    
 
  elif hour == "18":
    tiHHC1 = np.where(((i_tdrvC_hr >= 15) * (i_tdrvC_hr < 21)))
    siHHC1 = np.asarray(tiHHC1)
    iHHC1  = siHHC1.flatten()

    hdrvC_lat	  = tdrvC_lat[iHHC1]
    hdrvC_lon	  = tdrvC_lon[iHHC1]
    hdrvC_prs	  = tdrvC_prs[iHHC1]
    hdrvC_hgt	  = tdrvC_hgt[iHHC1]
    hdrvC_yr	  = tdrvC_yr [iHHC1]
    hdrvC_mm	  = tdrvC_mm [iHHC1]
    hdrvC_dy	  = tdrvC_dy [iHHC1]
    hdrvC_hr	  = tdrvC_hr [iHHC1]
    hdrvC_mn	  = tdrvC_mn [iHHC1]
    hdrvC_err	  = tdrvC_err[iHHC1]
    hdrvC_len	  = tdrvC_len[iHHC1]
    hdrvC_hgtTop  = tdrvC_hgtTop[iHHC1]
    hdrvC_hgtBot  = tdrvC_hgtBot[iHHC1]
    hdrvC_spd     = tdrvC_spd[iHHC1]
    hdrvC_dir     = tdrvC_dir[iHHC1]

    if dsetflag == "dep":
    	# keep ihrB before CURRENT and ihrA after CURRENT 6hrs
      if (21+ihrA) >= 24:
        dhrA = abs(24 - (21 + ihrA))

        tiHHA1 = np.where(((i_tdrvC_hr >= 21) * (i_tdrvC_hr < 24)))
        siHHA1 = np.asarray(tiHHA1)
        iHHA1  = siHHA1.flatten()
      
        if existA == 0:
          tiHHA2 = np.where(((i_tdrvA_hr >= 0) * (i_tdrvA_hr < dhrA)))  
          siHHA2 = np.asarray(tiHHA2)
          iHHA2  = siHHA2.flatten()
      
      tiHHB4 = np.where(((i_tdrvC_hr >= (15-ihrB)) * (i_tdrvC_hr < 15)))
      siHHB4 = np.asarray(tiHHB4)
      iHHB4  = siHHB4.flatten()

      hdrvA1_lat     = tdrvC_lat[iHHA1]
      hdrvA1_lon     = tdrvC_lon[iHHA1]
      hdrvA1_prs     = tdrvC_prs[iHHA1]
      hdrvA1_hgt     = tdrvC_hgt[iHHA1]
      hdrvA1_yr	     = tdrvC_yr [iHHA1]
      hdrvA1_mm	     = tdrvC_mm [iHHA1]
      hdrvA1_dy	     = tdrvC_dy [iHHA1]
      hdrvA1_hr	     = tdrvC_hr [iHHA1]
      hdrvA1_mn	     = tdrvC_mn [iHHA1]
      hdrvA1_err     = tdrvC_err[iHHA1]
      hdrvA1_len     = tdrvC_len[iHHA1]
      hdrvA1_hgtTop  = tdrvC_hgtTop[iHHA1]
      hdrvA1_hgtBot  = tdrvC_hgtBot[iHHA1]
      hdrvA1_spd     = tdrvC_spd[iHHA1]
      hdrvA1_dir     = tdrvC_dir[iHHA1]

      if existA == 0:
        hdrvA2_lat     = tdrvA_lat[iHHA2]
        hdrvA2_lon     = tdrvA_lon[iHHA2]
        hdrvA2_prs     = tdrvA_prs[iHHA2]
        hdrvA2_hgt     = tdrvA_hgt[iHHA2]
        hdrvA2_yr      = tdrvA_yr [iHHA2]
        hdrvA2_mm      = tdrvA_mm [iHHA2]
        hdrvA2_dy      = tdrvA_dy [iHHA2]
        hdrvA2_hr      = tdrvA_hr [iHHA2]
        hdrvA2_mn      = tdrvA_mn [iHHA2]
        hdrvA2_err     = tdrvA_err[iHHA2]
        hdrvA2_len     = tdrvA_len[iHHA2]
        hdrvA2_hgtTop  = tdrvA_hgtTop[iHHA2]
        hdrvA2_hgtBot  = tdrvA_hgtBot[iHHA2]
        hdrvA2_spd     = tdrvA_spd[iHHA2]
        hdrvA2_dir     = tdrvA_dir[iHHA2]
      
      hdrvB4_lat     = tdrvC_lat[iHHB4]
      hdrvB4_lon     = tdrvC_lon[iHHB4]
      hdrvB4_prs     = tdrvC_prs[iHHB4]
      hdrvB4_hgt     = tdrvC_hgt[iHHB4]
      hdrvB4_yr      = tdrvC_yr [iHHB4]
      hdrvB4_mm      = tdrvC_mm [iHHB4]
      hdrvB4_dy      = tdrvC_dy [iHHB4]
      hdrvB4_hr      = tdrvC_hr [iHHB4]
      hdrvB4_mn      = tdrvC_mn [iHHB4]
      hdrvB4_err     = tdrvC_err[iHHB4]
      hdrvB4_len     = tdrvC_len[iHHB4]
      hdrvB4_hgtTop  = tdrvC_hgtTop[iHHB4]
      hdrvB4_hgtBot  = tdrvC_hgtBot[iHHB4]
      hdrvB4_spd     = tdrvC_spd[iHHB4]
      hdrvB4_dir     = tdrvC_dir[iHHB4]

		# Append current arrays to before date arrays
      drv_lat    = np.append(hdrvB4_lat ,hdrvC_lat,axis=0)
      drv_lon    = np.append(hdrvB4_lon ,hdrvC_lon,axis=0)
      drv_prs    = np.append(hdrvB4_prs ,hdrvC_prs,axis=0)
      drv_hgt    = np.append(hdrvB4_hgt ,hdrvC_hgt,axis=0)
      drv_yr     = np.append(hdrvB4_yr  ,hdrvC_yr ,axis=0)
      drv_mm     = np.append(hdrvB4_mm  ,hdrvC_mm ,axis=0)
      drv_dy     = np.append(hdrvB4_dy  ,hdrvC_dy ,axis=0)
      drv_hr     = np.append(hdrvB4_hr  ,hdrvC_hr ,axis=0)
      drv_mn     = np.append(hdrvB4_mn  ,hdrvC_mn ,axis=0)
      drv_err	 = np.append(hdrvB4_err        ,hdrvC_err   ,axis=0)
      drv_len	 = np.append(hdrvB4_len        ,hdrvC_len   ,axis=0)
      drv_hgtTop = np.append(hdrvB4_hgtTop,hdrvC_hgtTop,axis=0)
      drv_hgtBot = np.append(hdrvB4_hgtBot,hdrvC_hgtBot,axis=0)
      drv_spd    = np.append(hdrvB4_spd ,hdrvC_spd,axis=0)
      drv_dir    = np.append(hdrvB4_dir ,hdrvC_dir,axis=0)
      
      drv_lat    = np.append(drv_lat ,hdrvA1_lat,axis=0)
      drv_lon    = np.append(drv_lon ,hdrvA1_lon,axis=0)
      drv_prs    = np.append(drv_prs ,hdrvA1_prs,axis=0)
      drv_hgt    = np.append(drv_hgt ,hdrvA1_hgt,axis=0)
      drv_yr     = np.append(drv_yr  ,hdrvA1_yr ,axis=0)
      drv_mm     = np.append(drv_mm  ,hdrvA1_mm ,axis=0)
      drv_dy     = np.append(drv_dy  ,hdrvA1_dy ,axis=0)
      drv_hr     = np.append(drv_hr  ,hdrvA1_hr ,axis=0)
      drv_mn     = np.append(drv_mn  ,hdrvA1_mn ,axis=0)
      drv_err	 = np.append(drv_err ,hdrvA1_err   ,axis=0)
      drv_len	 = np.append(drv_len ,hdrvA1_len   ,axis=0)
      drv_hgtTop = np.append(drv_hgtTop,hdrvA1_hgtTop,axis=0)
      drv_hgtBot = np.append(drv_hgtBot,hdrvA1_hgtBot,axis=0)
      drv_spd    = np.append(drv_spd ,hdrvA1_spd,axis=0)
      drv_dir    = np.append(drv_dir ,hdrvA1_dir,axis=0)
      
      if existA == 0:
        drv_lat    = np.append(drv_lat ,hdrvA2_lat,axis=0)
        drv_lon    = np.append(drv_lon ,hdrvA2_lon,axis=0)
        drv_prs    = np.append(drv_prs ,hdrvA2_prs,axis=0)
        drv_hgt    = np.append(drv_hgt ,hdrvA2_hgt,axis=0)
        drv_yr     = np.append(drv_yr  ,hdrvA2_yr ,axis=0)
        drv_mm     = np.append(drv_mm  ,hdrvA2_mm ,axis=0)
        drv_dy     = np.append(drv_dy  ,hdrvA2_dy ,axis=0)
        drv_hr     = np.append(drv_hr  ,hdrvA2_hr ,axis=0)
        drv_mn     = np.append(drv_mn  ,hdrvA2_mn ,axis=0)
        drv_err	 = np.append(drv_err ,hdrvA2_err   ,axis=0)
        drv_len	 = np.append(drv_len ,hdrvA2_len   ,axis=0)
        drv_hgtTop = np.append(drv_hgtTop,hdrvA2_hgtTop,axis=0)
        drv_hgtBot = np.append(drv_hgtBot,hdrvA2_hgtBot,axis=0)
        drv_spd    = np.append(drv_spd ,hdrvA2_spd,axis=0)
        drv_dir    = np.append(drv_dir ,hdrvA2_dir,axis=0)
      
    else:
    	# dsetflag = "drv"
      drv_lat    = hdrvC_lat
      drv_lon    = hdrvC_lon
      drv_prs    = hdrvC_prs
      drv_hgt    = hdrvC_hgt
      drv_yr     = hdrvC_yr 
      drv_mm     = hdrvC_mm 
      drv_dy     = hdrvC_dy 
      drv_hr     = hdrvC_hr 
      drv_mn     = hdrvC_mn 
      drv_err	 = hdrvC_err   
      drv_len	 = hdrvC_len   
      drv_hgtTop = hdrvC_hgtTop
      drv_hgtBot = hdrvC_hgtBot     
      drv_spd    = hdrvC_spd
      drv_dir    = hdrvC_dir

	#----------------------------------------
	# if 'runtype' = 'match', get indices of matches and apply QC if bool_drv_qc=True.
	# if 'runtype' = 'plot', SKIP bool if-block
	
  if runtype == "match":
    if bool_drv_qc:
      tindexesDC,qc_list = qc_aeolus(driver_wind_type,drv_prs,drv_err,drv_len,drv_hgtTop,drv_hgtBot)
      sindexesDC = np.asarray(tindexesDC)

      	#'indexesDC' = indices of appended array (not indices from Aeolus source files)
	#      It includes both current date (C) and before date (B) indices. This combo is saved in 
      	#      the collocation index files (for 00 UTC files only). The plotting code knows this and applies
      	#      the same approach as here (find 6-hr window, then append B and C indices) to find indices of 
      	#      collocated pairs.
      indexesDC  = sindexesDC.flatten()

      d_lat = drv_lat[indexesDC]
      d_lon = drv_lon[indexesDC]
      d_prs = drv_prs[indexesDC]
      d_hgt = drv_hgt[indexesDC]
      d_yr  = drv_yr[indexesDC]
      d_mm  = drv_mm[indexesDC]
      d_dy  = drv_dy[indexesDC]
      d_hr  = drv_hr[indexesDC]
      d_mn  = drv_mn[indexesDC]
      d_err = drv_err[indexesDC]
      d_len = drv_len[indexesDC]
      d_spd = drv_spd[indexesDC]
      d_dir = drv_dir[indexesDC]

    elif not bool_drv_qc:
      # Do not apply AEOLUS QC

      sindexesDC = np.asarray(np.where(drv_lat==drv_lat))	      #get all indices
      indexesDC  = sindexesDC.flatten()

      qc_list = "No QC applied"

      d_lat = drv_lat 
      d_lon = drv_lon 
      d_prs = drv_prs 
      d_hgt = drv_hgt 
      d_yr  = drv_yr  
      d_mm  = drv_mm  
      d_dy  = drv_dy  
      d_hr  = drv_hr  
      d_mn  = drv_mn 
      d_err = drv_err
      d_len = drv_len
      d_spd = drv_spd
      d_dir = drv_dir
      
  elif runtype == "plot":

    indexesDC = [0]

    qc_list = "No QC applied"

    d_lat = drv_lat [idxs]
    d_lon = drv_lon [idxs]
    d_prs = drv_prs [idxs]
    d_hgt = drv_hgt [idxs]
    d_yr  = drv_yr  [idxs]
    d_mm  = drv_mm  [idxs]
    d_dy  = drv_dy  [idxs]
    d_hr  = drv_hr  [idxs]
    d_mn  = drv_mn  [idxs]
    d_err = drv_err [idxs]
    d_len = drv_len [idxs]
    d_spd = drv_spd [idxs]
    d_dir = drv_dir [idxs]

	#----------------------------------------

  del drv_lat
  del drv_lon
  del drv_prs
  del drv_hgt
  del drv_yr 
  del drv_mm 
  del drv_dy 
  del drv_hr 
  del drv_mn
  del drv_err
  del drv_len
  del drv_hgtTop
  del drv_hgtBot
  del drv_spd
  del drv_dir

  indexesD     = indexesDC

  if np.size(indexesD)<=0:
    print("ERROR: No obs available to match on "+str(yyyymmddhh)+" in "+driver_file)
    sys.exit()
 
  	# Return variables to MAIN
  return d_lat,d_lon,d_prs,d_hgt,d_yr,d_mm,d_dy,d_hr,d_mn,indexesD,qc_list,drv_src,d_err,d_len,d_spd,d_dir

#===============================================================================================
# Read AIRCRAFT (from NCEP)
#	Don't QC
#
#	INPUTS:
#		yyyymmddhh .......................... Current date in yyyymmddhh format
#		bool_qc ......................... Choice to apply Aeolus QC: True=apply QC, False=don't apply QC
#
#	OUTPUTS:
#		d_lat ............................... Latitude in degrees [-90,90]
#		d_lon ............................... Longitude in degrees [0,360]
#		d_prs ............................... Pressure in hPa
#		d_hgt ............................... Height in km
#		d_yr ................................ Year
#		d_mm ................................ Month
#		d_dy ................................ Day
#		d_hr ................................ Hour
#		d_mn ................................ Minute
#		qc_list ............................. List of QC applied (if applicable)
#
def read_aircraft(path_prefix,yyyymmddhh,dateB4,dateA,bool_qc,dsetflag,runtype,time_diff_max,idxs):

  qc_list = ""			#initialize

  yyyy = yyyymmddhh[0:4]
  mm   = yyyymmddhh[4:6]
  dd   = yyyymmddhh[6:8]
  hour = yyyymmddhh[8:10]

  bufrtype = 'aircft'

	#-------------------------------------
        # Find hour limits ... for Aeolus datafiles

  if hour == "00":
    hB4 = "18"
    hA  = "06"
  elif hour == "06":
    hB4 = "00"
    hA  = "12"
  elif hour == "12":
    hB4 = "06"
    hA  = "18"
  elif hour == "18":
    hB4 = "12"
    hA  = "00"

  if hour == "00":
    yyyymmddhhB4 = dateB4+hB4
    yyB4 = dateB4[0:4]
    mmB4 = dateB4[4:6]
    ddB4 = dateB4[6:8]
  else:
    yyyymmddhhB4 = yyyy+mm+dd+hB4
    yyB4 = yyyy
    mmB4 = mm
    ddB4 = dd
    
  if hour == "18":
    yyyymmddhhA = dateA+hA
    yyA  = dateA[0:4]
    mmA  = dateA[4:6]
    ddA  = dateA[6:8]
  else:
    yyyymmddhhA = yyyy+mm+dd+hA
    yyA = yyyy
    mmA = mm
    ddA = dd

		# convert 'time_max_diff' to rounded integer
  tqc_time = float(time_diff_max)
  if tqc_time%60 != 0.0:
    qc_time = int(math.ceil(tqc_time/60.0))             # round up to nearest or equal integer hour
  else:
    qc_time = int(tqc_time)

		# convert hour to integer
  ihr = int(hour)

		# find integer hours for before/after files for collocation
  if hour == "00":
    	# 00
    ihrB = 24 - qc_time
    ihrA = ihr + qc_time
  else:
	# 06,12,18
    ihrB = ihr - qc_time
    ihrA = ihr + qc_time

  del qc_time,tqc_time,ihr

	#-------------------------------------------------
    	# Define dataset

  tmp_dset1_path    = '/scratch/atmos-nc-dataset/aircraft/'+yyyy+'/'+mm+'/'+dd+'/'
  tmp_dset1_path_B4 = '/scratch/atmos-nc-dataset/aircraft/'+yyB4+'/'+mmB4+'/'+ddB4+'/'
  tmp_dset1_path_A  = '/scratch/atmos-nc-dataset/aircraft/'+yyA+'/'+mmA+'/'+ddA+'/'

		# Path/file
		# ...Current date
  dset1_path   		= path_prefix+tmp_dset1_path
  dset1_filename 	= 'gdas.'+yyyymmddhh+'.'+str(bufrtype)+'.tm00.bufr_d.nc4'
  dset1_file   		= dset1_path+dset1_filename
  		# ...Before date
  dset1_path_B4		= path_prefix+tmp_dset1_path_B4
  dset1_filename_B4 	= 'gdas.'+yyyymmddhhB4+'.'+str(bufrtype)+'.tm00.bufr_d.nc4'
  dset1_file_B4		= dset1_path_B4+dset1_filename_B4
  		# ...After date
  dset1_path_A 		= path_prefix+tmp_dset1_path_A
  dset1_filename_A 	= 'gdas.'+yyyymmddhhA+'.'+str(bufrtype)+'.tm00.bufr_d.nc4'
  dset1_file_A 		= dset1_path_A+dset1_filename_A
  
  			# initialize flag indicating if dataset exists. 0=yes, 1=no
  existB = 0		# ... B = date Before current
  existA = 0		# ... A = date After current
  dset1_exists = exists(dset1_file)
  if dset1_exists==False:
    print("ERROR: file "+dset1_file+" does not exist!")
    sys.exit()							#exit script immediately
  dset1_exists_B4 = exists(dset1_file_B4)
  if dset1_exists_B4==False:
    print("WARNING: 'before' file "+dset1_file_B4+" does not exist!")
    existB = 1		# set flag to 1=file doesn't exist
  dset1_exists_A = exists(dset1_file_A)
  if dset1_exists_A==False:
    print("WARNING: 'after' file "+dset1_file_A+" does not exist!")
    existA = 1		# set flag to 1=file doesn't exist

		# Path on FTP/web archive server (for output NetCDF only)
  str_dset1_path    = dset1_path
  str_dset1_path_B4 = dset1_path_B4
  str_dset1_path_A  = dset1_path_A

  dset1_src  = str_dset1_path_B4+dset1_filename_B4
  dset1_src += ", "+str_dset1_path+dset1_filename
  dset1_src += ", "+str_dset1_path_A+dset1_filename_A

		# Variable names
  dset1_lat_var = 'latitude'
  dset1_lon_var = 'longitude'
  dset1_prs_var = 'pressure'
  dset1_yr_var  = 'year'
  dset1_mm_var  = 'month'
  dset1_dy_var  = 'day'
  dset1_hr_var  = 'hour'
  dset1_mn_var  = 'minutes'
  
  dset1_ht1_var = 'flight_level'
  dset1_ht2_var = 'height'
  dset1_ht3_var = 'aircraft_altitude'
  
  dset1_spd_var = 'wind_speed'
  dset1_dir_var = 'wind_direction'
		
    	#-------------------------------------------------
  	# Load dataset
	#	AIRCRAFT data is divided into Groups
  
  		#```````````````````````````````````````````````
  		# CURRENT date
  data_hdl = Dataset(dset1_file)

  grps = list(data_hdl.groups)

	# populate full arrays with Group1 (grps(0))
  tdC_lat = np.asarray( data_hdl.groups[grps[0]].variables[dset1_lat_var] )
  tdC_lon = np.asarray( data_hdl.groups[grps[0]].variables[dset1_lon_var] )
  tdC_yr  = np.asarray( data_hdl.groups[grps[0]].variables[dset1_yr_var] )
  tdC_mm  = np.asarray( data_hdl.groups[grps[0]].variables[dset1_mm_var] )
  tdC_dy  = np.asarray( data_hdl.groups[grps[0]].variables[dset1_dy_var] )
  tdC_hr  = np.asarray( data_hdl.groups[grps[0]].variables[dset1_hr_var] )
  tdC_mn  = np.asarray( data_hdl.groups[grps[0]].variables[dset1_mn_var] )

  if grps[0]=='NC004006' or grps[0]=='NC004009':
    tdC_hgt = np.asarray( data_hdl.groups[grps[0]].variables[dset1_ht2_var]  )
  elif grps[0]=='NC004004': 
    tdC_hgt = np.asarray( data_hdl.groups[grps[0]].variables[dset1_ht3_var]  )
  else:
    tdC_hgt = np.asarray( data_hdl.groups[grps[0]].variables[dset1_ht1_var]  )
  
  tdC_spd = np.asarray( data_hdl.groups[grps[0]].variables[dset1_spd_var]  )
  tdC_dir = np.asarray( data_hdl.groups[grps[0]].variables[dset1_dir_var]  )

	# append data from remaining Groups
  for x in range(len(grps)-1):
      tdset1_lat = np.asarray( data_hdl.groups[grps[x+1]].variables[dset1_lat_var] )
      tdset1_lon = np.asarray( data_hdl.groups[grps[x+1]].variables[dset1_lon_var] )
      tdset1_yr  = np.asarray( data_hdl.groups[grps[x+1]].variables[dset1_yr_var]  )
      tdset1_mm  = np.asarray( data_hdl.groups[grps[x+1]].variables[dset1_mm_var]  )
      tdset1_dy  = np.asarray( data_hdl.groups[grps[x+1]].variables[dset1_dy_var]  )
      tdset1_hr  = np.asarray( data_hdl.groups[grps[x+1]].variables[dset1_hr_var]  )
      tdset1_mn  = np.asarray( data_hdl.groups[grps[x+1]].variables[dset1_mn_var]  )
      tdset1_spd = np.asarray( data_hdl.groups[grps[x+1]].variables[dset1_spd_var]  )
      tdset1_dir = np.asarray( data_hdl.groups[grps[x+1]].variables[dset1_dir_var]  )

      tdC_lat = np.append(tdC_lat,tdset1_lat,axis=0)
      tdC_lon = np.append(tdC_lon,tdset1_lon,axis=0)
      tdC_yr  = np.append(tdC_yr, tdset1_yr, axis=0)
      tdC_mm  = np.append(tdC_mm, tdset1_mm, axis=0)
      tdC_dy  = np.append(tdC_dy, tdset1_dy, axis=0)
      tdC_hr  = np.append(tdC_hr, tdset1_hr, axis=0)
      tdC_mn  = np.append(tdC_mn, tdset1_mn, axis=0)
      tdC_spd = np.append(tdC_spd,tdset1_spd,axis=0)
      tdC_dir = np.append(tdC_dir,tdset1_dir,axis=0)

      if grps[x+1]=='NC004006' or grps[x+1]=='NC004009':
        z = np.asarray( data_hdl.groups[grps[x+1]].variables[dset1_ht2_var]  )
      elif grps[x+1]=='NC004004':
        z = np.asarray( data_hdl.groups[grps[x+1]].variables[dset1_ht3_var]  )
      else:
        z = np.asarray( data_hdl.groups[grps[x+1]].variables[dset1_ht1_var]  )

      tdC_hgt = np.append(tdC_hgt,z,axis=0)

      del tdset1_lat,tdset1_lon,tdset1_yr,tdset1_mm,tdset1_dy,tdset1_hr,tdset1_mn,tdset1_spd,tdset1_dir,z

  data_hdl.close()
  del grps

  if dsetflag == "dep":
    if existB == 0:
		#```````````````````````````````````````````````
  		# BEFORE date
      data_hdl = Dataset(dset1_file_B4)

      grps = list(data_hdl.groups)

		# populate full arrays with Group1 (grps(0))
      tdB_lat = np.asarray( data_hdl.groups[grps[0]].variables[dset1_lat_var] )
      tdB_lon = np.asarray( data_hdl.groups[grps[0]].variables[dset1_lon_var] )
      tdB_yr  = np.asarray( data_hdl.groups[grps[0]].variables[dset1_yr_var] )
      tdB_mm  = np.asarray( data_hdl.groups[grps[0]].variables[dset1_mm_var] )
      tdB_dy  = np.asarray( data_hdl.groups[grps[0]].variables[dset1_dy_var] )
      tdB_hr  = np.asarray( data_hdl.groups[grps[0]].variables[dset1_hr_var] )
      tdB_mn  = np.asarray( data_hdl.groups[grps[0]].variables[dset1_mn_var] )
      
      if grps[0]=='NC004006' or grps[0]=='NC004009':
        tdB_hgt = np.asarray( data_hdl.groups[grps[0]].variables[dset1_ht2_var]  )
      else:
        tdB_hgt = np.asarray( data_hdl.groups[grps[0]].variables[dset1_ht1_var]  )
      
      tdB_spd = np.asarray( data_hdl.groups[grps[0]].variables[dset1_spd_var]  )
      tdB_dir = np.asarray( data_hdl.groups[grps[0]].variables[dset1_dir_var]  )
      
          	# append data from remaining Groups
      for x in range(len(grps)-1):
        tdset1_lat = np.asarray( data_hdl.groups[grps[x+1]].variables[dset1_lat_var] )
        tdset1_lon = np.asarray( data_hdl.groups[grps[x+1]].variables[dset1_lon_var] )
        tdset1_yr  = np.asarray( data_hdl.groups[grps[x+1]].variables[dset1_yr_var]  )
        tdset1_mm  = np.asarray( data_hdl.groups[grps[x+1]].variables[dset1_mm_var]  )
        tdset1_dy  = np.asarray( data_hdl.groups[grps[x+1]].variables[dset1_dy_var]  )
        tdset1_hr  = np.asarray( data_hdl.groups[grps[x+1]].variables[dset1_hr_var]  )
        tdset1_mn  = np.asarray( data_hdl.groups[grps[x+1]].variables[dset1_mn_var]  )
        tdset1_spd = np.asarray( data_hdl.groups[grps[x+1]].variables[dset1_spd_var]  )
        tdset1_dir = np.asarray( data_hdl.groups[grps[x+1]].variables[dset1_dir_var]  )

        tdB_lat = np.append(tdB_lat,tdset1_lat,axis=0)
        tdB_lon = np.append(tdB_lon,tdset1_lon,axis=0)
        tdB_yr  = np.append(tdB_yr, tdset1_yr, axis=0)
        tdB_mm  = np.append(tdB_mm, tdset1_mm, axis=0)
        tdB_dy  = np.append(tdB_dy, tdset1_dy, axis=0)
        tdB_hr  = np.append(tdB_hr, tdset1_hr, axis=0)
        tdB_mn  = np.append(tdB_mn, tdset1_mn, axis=0)
        tdB_spd = np.append(tdB_spd, tdset1_spd, axis=0)
        tdB_dir = np.append(tdB_dir, tdset1_dir, axis=0)
      
        if grps[x+1]=='NC004006' or grps[x+1]=='NC004009':
          z = np.asarray( data_hdl.groups[grps[x+1]].variables[dset1_ht2_var]  )
        elif grps[x+1]=='NC004004':
          z = np.asarray( data_hdl.groups[grps[x+1]].variables[dset1_ht3_var]  )
        else:
          z = np.asarray( data_hdl.groups[grps[x+1]].variables[dset1_ht1_var]  )
      
        tdB_hgt = np.append(tdB_hgt,z,axis=0)
        
        del tdset1_lat,tdset1_lon,tdset1_yr,tdset1_mm,tdset1_dy,tdset1_hr,tdset1_mn,tdset1_spd,tdset1_dir,z
      
      data_hdl.close()
      del grps

    if existA == 0:
        	#```````````````````````````````````````````````
        	# AFTER date
      data_hdl = Dataset(dset1_file_A)

      grps = list(data_hdl.groups)

		# populate full arrays with Group1 (grps(0))
      tdA_lat = np.asarray( data_hdl.groups[grps[0]].variables[dset1_lat_var] )
      tdA_lon = np.asarray( data_hdl.groups[grps[0]].variables[dset1_lon_var] )
      tdA_yr  = np.asarray( data_hdl.groups[grps[0]].variables[dset1_yr_var] )
      tdA_mm  = np.asarray( data_hdl.groups[grps[0]].variables[dset1_mm_var] )
      tdA_dy  = np.asarray( data_hdl.groups[grps[0]].variables[dset1_dy_var] )
      tdA_hr  = np.asarray( data_hdl.groups[grps[0]].variables[dset1_hr_var] )
      tdA_mn  = np.asarray( data_hdl.groups[grps[0]].variables[dset1_mn_var] )
      
      if grps[0]=='NC004006' or grps[0]=='NC004009':
        tdA_hgt = np.asarray( data_hdl.groups[grps[0]].variables[dset1_ht2_var]  )
      else:
        tdA_hgt = np.asarray( data_hdl.groups[grps[0]].variables[dset1_ht1_var]  )
      
      tdA_spd = np.asarray( data_hdl.groups[grps[0]].variables[dset1_spd_var]  )
      tdA_dir = np.asarray( data_hdl.groups[grps[0]].variables[dset1_dir_var]  )
      
          	# append data from remaining Groups
      for x in range(len(grps)-1):
        tdset1_lat = np.asarray( data_hdl.groups[grps[x+1]].variables[dset1_lat_var] )
        tdset1_lon = np.asarray( data_hdl.groups[grps[x+1]].variables[dset1_lon_var] )
        tdset1_yr  = np.asarray( data_hdl.groups[grps[x+1]].variables[dset1_yr_var]  )
        tdset1_mm  = np.asarray( data_hdl.groups[grps[x+1]].variables[dset1_mm_var]  )
        tdset1_dy  = np.asarray( data_hdl.groups[grps[x+1]].variables[dset1_dy_var]  )
        tdset1_hr  = np.asarray( data_hdl.groups[grps[x+1]].variables[dset1_hr_var]  )
        tdset1_mn  = np.asarray( data_hdl.groups[grps[x+1]].variables[dset1_mn_var]  )
        tdset1_spd = np.asarray( data_hdl.groups[grps[x+1]].variables[dset1_spd_var]  )
        tdset1_dir = np.asarray( data_hdl.groups[grps[x+1]].variables[dset1_dir_var]  )

        tdA_lat = np.append(tdA_lat,tdset1_lat,axis=0)
        tdA_lon = np.append(tdA_lon,tdset1_lon,axis=0)
        tdA_yr  = np.append(tdA_yr, tdset1_yr, axis=0)
        tdA_mm  = np.append(tdA_mm, tdset1_mm, axis=0)
        tdA_dy  = np.append(tdA_dy, tdset1_dy, axis=0)
        tdA_hr  = np.append(tdA_hr, tdset1_hr, axis=0)
        tdA_mn  = np.append(tdA_mn, tdset1_mn, axis=0)
        tdA_spd = np.append(tdA_spd, tdset1_spd, axis=0)
        tdA_dir = np.append(tdA_dir, tdset1_dir, axis=0)

        if grps[x+1]=='NC004006' or grps[x+1]=='NC004009':
          z = np.asarray( data_hdl.groups[grps[x+1]].variables[dset1_ht2_var]  )
        elif grps[x+1]=='NC004004':
          z = np.asarray( data_hdl.groups[grps[x+1]].variables[dset1_ht3_var]  )
        else:
          z = np.asarray( data_hdl.groups[grps[x+1]].variables[dset1_ht1_var]  )

        tdA_hgt = np.append(tdA_hgt,z,axis=0)
        
        del tdset1_lat,tdset1_lon,tdset1_yr,tdset1_mm,tdset1_dy,tdset1_hr,tdset1_mn,tdset1_spd,tdset1_dir,z

      data_hdl.close()
      del grps

	#-------------------------------------------------
      	# Append current arrays to before date arrays

		# get BEFORE and AFTER obs
    if hour=="00":
    	# keep ihrB hrs before CURRENT and ihrA hrs after CURRENT 6hrs
      tiHHA  = np.where(((tdC_hr >= 3) * (tdC_hr < (3+ihrA))))
      siHHA  = np.asarray(tiHHA)
      iHHA   = siHHA.flatten()

      stdA_lat = tdC_lat[iHHA]
      stdA_lon = tdC_lon[iHHA]
      stdA_yr  = tdC_yr [iHHA]
      stdA_mm  = tdC_mm [iHHA]
      stdA_dy  = tdC_dy [iHHA]
      stdA_hr  = tdC_hr [iHHA]
      stdA_mn  = tdC_mn [iHHA]
      stdA_spd = tdC_spd[iHHA]
      stdA_dir = tdC_dir[iHHA]
      stdA_hgt = tdC_hgt[iHHA]
	
      if existB == 0:
        tiHHB4 = np.where(((tdB_hr >= (21-ihrB)) * (tdB_hr < 21)))
        siHHB4 = np.asarray(tiHHB4)
        iHHB4  = siHHB4.flatten()
	
        stdB_lat = tdB_lat[iHHB4]
        stdB_lon = tdB_lon[iHHB4]
        stdB_yr  = tdB_yr [iHHB4]
        stdB_mm  = tdB_mm [iHHB4]
        stdB_dy  = tdB_dy [iHHB4]
        stdB_hr  = tdB_hr [iHHB4]
        stdB_mn  = tdB_mn [iHHB4]
        stdB_spd = tdB_spd[iHHB4]
        stdB_dir = tdB_dir[iHHB4]
        stdB_hgt = tdB_hgt[iHHB4]
	
    elif hour=="06":
    	# keep ihrB before CURRENT and ihrA after CURRENT 6hrs
      tiHHA   = np.where(((tdC_hr >= 9) * (tdC_hr < (9+ihrA))))
      siHHA   = np.asarray(tiHHA)
      iHHA    = siHHA.flatten()
      
      stdA_lat = tdC_lat[iHHA]
      stdA_lon = tdC_lon[iHHA]
      stdA_yr  = tdC_yr [iHHA]
      stdA_mm  = tdC_mm [iHHA]
      stdA_dy  = tdC_dy [iHHA]
      stdA_hr  = tdC_hr [iHHA]
      stdA_mn  = tdC_mn [iHHA]
      stdA_spd = tdC_spd[iHHA]
      stdA_dir = tdC_dir[iHHA]
      stdA_hgt = tdC_hgt[iHHA]
      
      if (3-ihrB) < 0:
        dhrB = abs(3 - ihrB)
      
        tiHHB4 = np.where(((tdC_hr >= 0) * (tdC_hr < 3)))
        siHHB4 = np.asarray(tiHHB4)
        iHHB4  = siHHB4.flatten()
      
        if existB == 0:
          tiHHB41 = np.where(((tdB_hr >= (24-dhrB)) * (tdB_hr < 24)))  
          siHHB41 = np.asarray(tiHHB41)
          iHHB41  = siHHB41.flatten()
	  
          stdB_lat = np.append(tdB_lat[iHHB41], tdC_lat[iHHB4], axis=0)
          stdB_lon = np.append(tdB_lon[iHHB41], tdC_lon[iHHB4], axis=0)
          stdB_yr  = np.append(tdB_yr [iHHB41], tdC_yr [iHHB4], axis=0)
          stdB_mm  = np.append(tdB_mm [iHHB41], tdC_mm [iHHB4], axis=0)
          stdB_dy  = np.append(tdB_dy [iHHB41], tdC_dy [iHHB4], axis=0)
          stdB_hr  = np.append(tdB_hr [iHHB41], tdC_hr [iHHB4], axis=0)
          stdB_mn  = np.append(tdB_mn [iHHB41], tdC_mn [iHHB4], axis=0)
          stdB_spd = np.append(tdB_spd[iHHB41], tdC_spd[iHHB4], axis=0)
          stdB_dir = np.append(tdB_dir[iHHB41], tdC_dir[iHHB4], axis=0)
          stdB_hgt = np.append(tdB_hgt[iHHB41], tdC_hgt[iHHB4], axis=0)
        else:
          stdB_lat = tdC_lat[iHHB4]
          stdB_lon = tdC_lon[iHHB4]
          stdB_yr  = tdC_yr [iHHB4]
          stdB_mm  = tdC_mm [iHHB4]
          stdB_dy  = tdC_dy [iHHB4]
          stdB_hr  = tdC_hr [iHHB4]
          stdB_mn  = tdC_mn [iHHB4]
          stdB_spd = tdC_spd[iHHB4]
          stdB_dir = tdC_dir[iHHB4]
          stdB_hgt = tdC_hgt[iHHB4]

        del dhrB
	
      else:
        tiHHB4 = np.where(((tdC_hr >= (3-ihrB)) * (tdC_hr < 3)))
        siHHB4 = np.asarray(tiHHB4)
        iHHB4  = siHHB4.flatten()
	
        stdB_lat = tdC_lat[iHHB4]
        stdB_lon = tdC_lon[iHHB4]
        stdB_yr  = tdC_yr [iHHB4]
        stdB_mm  = tdC_mm [iHHB4]
        stdB_dy  = tdC_dy [iHHB4]
        stdB_hr  = tdC_hr [iHHB4]
        stdB_mn  = tdC_mn [iHHB4]
        stdB_spd = tdC_spd[iHHB4]
        stdB_dir = tdC_dir[iHHB4]
        stdB_hgt = tdC_hgt[iHHB4]
	
    elif hour=="12":
    	# keep ihrB before CURRENT and ihrA after CURRENT 6hrs
      tiHHA  = np.where(((tdC_hr >= 15) * (tdC_hr < (15+ihrA)))) 
      tiHHB4 = np.where(((tdC_hr >= (9-ihrB)) * (tdC_hr < 9)))
    
      siHHA  = np.asarray(tiHHA)
      iHHA   = siHHA.flatten()
      siHHB4 = np.asarray(tiHHB4)
      iHHB4  = siHHB4.flatten()
      
      stdA_lat = tdC_lat[iHHA]
      stdA_lon = tdC_lon[iHHA]
      stdA_yr  = tdC_yr [iHHA]
      stdA_mm  = tdC_mm [iHHA]
      stdA_dy  = tdC_dy [iHHA]
      stdA_hr  = tdC_hr [iHHA]
      stdA_mn  = tdC_mn [iHHA]
      stdA_spd = tdC_spd[iHHA]
      stdA_dir = tdC_dir[iHHA]
      stdA_hgt = tdC_hgt[iHHA]
      
      stdB_lat = tdC_lat[iHHB4]
      stdB_lon = tdC_lon[iHHB4]
      stdB_yr  = tdC_yr [iHHB4]
      stdB_mm  = tdC_mm [iHHB4]
      stdB_dy  = tdC_dy [iHHB4]
      stdB_hr  = tdC_hr [iHHB4]
      stdB_mn  = tdC_mn [iHHB4]
      stdB_spd = tdC_spd[iHHB4]
      stdB_dir = tdC_dir[iHHB4]
      stdB_hgt = tdC_hgt[iHHB4]
    
    elif hour=="18":
    	# keep ihrB before CURRENT and ihrA after CURRENT 6hrs
      if (21+ihrA) >= 24:
        dhrA = abs(24 - (21 + ihrA))

        tiHHA = np.where(((tdC_hr >= 21) * (tdC_hr < 24)))
        siHHA = np.asarray(tiHHA)
        iHHA  = siHHA.flatten()
      
        if existA == 0:
          tiHHA2 = np.where(((tdA_hr >= 0) * (tdA_hr < dhrA)))  
          siHHA2 = np.asarray(tiHHA2)
          iHHA2  = siHHA2.flatten()
	  
          stdA_lat = np.append(tdC_lat[iHHA], tdA_lat[iHHA2], axis=0)
          stdA_lon = np.append(tdC_lon[iHHA], tdA_lon[iHHA2], axis=0)
          stdA_yr  = np.append(tdC_yr [iHHA], tdA_yr [iHHA2], axis=0)
          stdA_mm  = np.append(tdC_mm [iHHA], tdA_mm [iHHA2], axis=0)
          stdA_dy  = np.append(tdC_dy [iHHA], tdA_dy [iHHA2], axis=0)
          stdA_hr  = np.append(tdC_hr [iHHA], tdA_hr [iHHA2], axis=0)
          stdA_mn  = np.append(tdC_mn [iHHA], tdA_mn [iHHA2], axis=0)
          stdA_spd = np.append(tdC_spd[iHHA], tdA_spd[iHHA2], axis=0)
          stdA_dir = np.append(tdC_dir[iHHA], tdA_dir[iHHA2], axis=0)
          stdA_hgt = np.append(tdC_hgt[iHHA], tdA_hgt[iHHA2], axis=0)
        else:
          stdA_lat = tdC_lat[iHHA]
          stdA_lon = tdC_lon[iHHA]
          stdA_yr  = tdC_yr [iHHA]
          stdA_mm  = tdC_mm [iHHA]
          stdA_dy  = tdC_dy [iHHA]
          stdA_hr  = tdC_hr [iHHA]
          stdA_mn  = tdC_mn [iHHA]
          stdA_spd = tdC_spd[iHHA]
          stdA_dir = tdC_dir[iHHA]
          stdA_hgt = tdC_hgt[iHHA]
      
      tiHHB4 = np.where(((tdC_hr >= (15-ihrB)) * (tdC_hr < 15)))
      siHHB4 = np.asarray(tiHHB4)
      iHHB4  = siHHB4.flatten()
    
      stdB_lat = tdC_lat[iHHB4]
      stdB_lon = tdC_lon[iHHB4]
      stdB_yr  = tdC_yr [iHHB4]
      stdB_mm  = tdC_mm [iHHB4]
      stdB_dy  = tdC_dy [iHHB4]
      stdB_hr  = tdC_hr [iHHB4]
      stdB_mn  = tdC_mn [iHHB4]
      stdB_spd = tdC_spd[iHHB4]
      stdB_dir = tdC_dir[iHHB4]
      stdB_hgt = tdC_hgt[iHHB4]
    
		# append CURRENT to BEFORE
    if existB == 0:
      td_lat = np.append(stdB_lat,tdC_lat,axis=0)
      td_lon = np.append(stdB_lon,tdC_lon,axis=0)
      td_yr  = np.append(stdB_yr ,tdC_yr ,axis=0)
      td_mm  = np.append(stdB_mm ,tdC_mm ,axis=0)
      td_dy  = np.append(stdB_dy ,tdC_dy ,axis=0)
      td_hr  = np.append(stdB_hr ,tdC_hr ,axis=0)
      td_mn  = np.append(stdB_mn ,tdC_mn ,axis=0)
      td_spd = np.append(stdB_spd,tdC_spd,axis=0)
      td_dir = np.append(stdB_dir,tdC_dir,axis=0)
      td_hgt = np.append(stdB_hgt,tdC_hgt,axis=0)
    else:
      td_lat = tdC_lat
      td_lon = tdC_lon
      td_yr  = tdC_yr 
      td_mm  = tdC_mm 
      td_dy  = tdC_dy 
      td_hr  = tdC_hr 
      td_mn  = tdC_mn 
      td_spd = tdC_spd
      td_dir = tdC_dir
      td_hgt = tdC_hgt

  		# append AFTER to CURRENT
    if existA == 0:
      td_lat = np.append(td_lat,stdA_lat,axis=0)
      td_lon = np.append(td_lon,stdA_lon,axis=0)
      td_yr  = np.append(td_yr ,stdA_yr ,axis=0)
      td_mm  = np.append(td_mm ,stdA_mm ,axis=0)
      td_dy  = np.append(td_dy ,stdA_dy ,axis=0)
      td_hr  = np.append(td_hr ,stdA_hr ,axis=0)
      td_mn  = np.append(td_mn ,stdA_mn ,axis=0)
      td_spd = np.append(td_spd,stdA_spd,axis=0)
      td_dir = np.append(td_dir,stdA_dir,axis=0)
      td_hgt = np.append(td_hgt,stdA_hgt,axis=0)

  else:
  	# dsetflag = "drv"
    td_lat = tdC_lat
    td_lon = tdC_lon
    td_yr  = tdC_yr 
    td_mm  = tdC_mm 
    td_dy  = tdC_dy 
    td_hr  = tdC_hr 
    td_mn  = tdC_mn 
    td_spd = tdC_spd
    td_dir = tdC_dir
    td_hgt = tdC_hgt  
  
   	#----------------------------------------
	# if 'runtype' = 'match', get indices of matches and apply QC if bool_drv_qc=True.
	# if 'runtype' = 'plot', SKIP bool if-block
	
  if runtype == "match":
    if bool_qc:
      print("Apply AIRCRAFT QC: TBA")
    elif not bool_qc:
      # Do not apply AIRCRAFT QC

      sindexesDC = np.asarray(np.where(td_lat==td_lat))           #get all indices
      indexes1   = sindexesDC.flatten()

      qc_list = "No QC applied"

      d_lat = td_lat
      d_lon = td_lon
      d_yr  = td_yr
      d_mm  = td_mm
      d_dy  = td_dy
      d_hr  = td_hr
      d_mn  = td_mn
      d_hgt = td_hgt
      d_spd = td_spd
      d_dir = td_dir

  elif runtype == "plot":

    indexes1 = [0]

    qc_list = "No QC applied"

    d_lat = td_lat[idxs]
    d_lon = td_lon[idxs]
    d_yr  = td_yr [idxs]
    d_mm  = td_mm [idxs]
    d_dy  = td_dy [idxs]
    d_hr  = td_hr [idxs]
    d_mn  = td_mn [idxs]
    d_hgt = td_hgt[idxs]
    d_spd = td_spd[idxs]
    d_dir = td_dir[idxs]

	#----------------------------------------
	
	# create pressure array but fill with missing -999.
	#	AIRCRAFT pressures not available.
  d_prs = np.nan * np.ones_like(d_lat)

	# check height units and convert to km
  if max(d_hgt) > 1000.:
    d_hgt = d_hgt/1000.

	# Return variables to MAIN
  return d_lat,d_lon,d_yr,d_mm,d_dy,d_hr,d_mn,d_hgt,d_prs,indexes1,qc_list,dset1_src,d_spd,d_dir 

#===============================================================================================
# Read AMV (from NCEP)
#
#	INPUTS:
#		yyyymmddhh .......................... Current date in yyyymmddhh format
#		pct ................................. Minimum AMV quality indicator (QI) in % for QC
#		bool_qc ............................. Choice to apply Aeolus QC: True=apply QC, False=don't apply QC
#
#	OUTPUTS:
#		d_lat ............................... Latitude in degrees [-90,90]
#		d_lon ............................... Longitude in degrees [0,360]
#		d_prs ............................... Pressure in hPa
#		d_hgt ............................... Height in km
#		d_yr ................................ Year
#		d_mm ................................ Month
#		d_dy ................................ Day
#		d_hr ................................ Hour
#		d_mn ................................ Minute
#		d_pccf .............................. AMV QI (percent confidence) in %
#               indexes2 ............................ Indices of input obs
#		qc_list ............................. List of QC applied (if applicable)
#
def read_amv_ncep(path_prefix,yyyymmddhh,dateB4,dateA,bool_qc,pct,qi_choice,dsetflag,runtype,time_diff_max,idxs):

  qc_list = ""			#initialize

  yyyy = yyyymmddhh[0:4]
  mm   = yyyymmddhh[4:6]
  dd   = yyyymmddhh[6:8]
  hour = yyyymmddhh[8:10]

	#-------------------------------------
        # Find hour limits ... for Aeolus datafiles

  if hour == "00":
    hB4 = "18"
    hA  = "06"
  elif hour == "06":
    hB4 = "00"
    hA  = "12"
  elif hour == "12":
    hB4 = "06"
    hA  = "18"
  elif hour == "18":
    hB4 = "12"
    hA  = "00"

  if hour == "00":
    yyyymmddhhB4 = dateB4+hB4
    yyB4 = dateB4[0:4]
    mmB4 = dateB4[4:6]
    ddB4 = dateB4[6:8]
  else:
    yyyymmddhhB4 = yyyy+mm+dd+hB4
    yyB4 = yyyy
    mmB4 = mm
    ddB4 = dd
    
  if hour == "18":
    yyyymmddhhA = dateA+hA
    yyA  = dateA[0:4]
    mmA  = dateA[4:6]
    ddA  = dateA[6:8]
  else:
    yyyymmddhhA = yyyy+mm+dd+hA
    yyA = yyyy
    mmA = mm
    ddA = dd

		# convert 'time_max_diff' to rounded integer
  tqc_time = float(time_diff_max)
  if tqc_time%60 != 0.0:
    qc_time = int(math.ceil(tqc_time/60.0))             # round up to nearest or equal integer hour
  else:
    qc_time = int(tqc_time)

		# convert hour to integer
  ihr = int(hour)

		# find integer hours for before/after files for collocation
  if hour == "00":
    	# 00
    ihrB = 24 - qc_time
    ihrA = ihr + qc_time
  else:
	# 06,12,18
    ihrB = ihr - qc_time
    ihrA = ihr + qc_time

  del qc_time,tqc_time,ihr

	#-------------------------------------------------
    	# Define dataset

  tmp_dset2_path    = '/scratch/atmos-nc-dataset/AMV/'+yyyy+'/'+mm+'/'+dd+'/'
  tmp_dset2_path_B4 = '/scratch/atmos-nc-dataset/AMV/'+yyB4+'/'+mmB4+'/'+ddB4+'/'
  tmp_dset2_path_A  = '/scratch/atmos-nc-dataset/AMV/'+yyA+'/'+mmA+'/'+ddA+'/'

		# Path/file
		# ...Current date
  dset2_path   		= path_prefix+tmp_dset2_path
  dset2_filename 	= 'gdas.'+yyyymmddhh+'.satwnd.tm00.bufr_d.nc4'
  dset2_file   		= dset2_path+dset2_filename
  		# ...Before date
  dset2_path_B4		= path_prefix+tmp_dset2_path_B4
  dset2_filename_B4 	= 'gdas.'+yyyymmddhhB4+'.satwnd.tm00.bufr_d.nc4'
  dset2_file_B4		= dset2_path_B4+dset2_filename_B4
  		# ...After date
  dset2_path_A		= path_prefix+tmp_dset2_path_A
  dset2_filename_A 	= 'gdas.'+yyyymmddhhA+'.satwnd.tm00.bufr_d.nc4'
  dset2_file_A		= dset2_path_A+dset2_filename_A

  			# initialize flag indicating if dataset exists. 0=yes, 1=no
  existB = 0		# ... B = date Before current
  existA = 0		# ... A = date After current
  dset2_exists = exists(dset2_file)
  if dset2_exists==False:
    print("ERROR: file "+dset2_file+" does not exist!")
    sys.exit()
  dset2_existsB = exists(dset2_file_B4)
  if dset2_existsB==False:
    print("WARNING: 'before' file "+dset2_file_B4+" does not exist!")
    existB = 1
  dset2_existsA = exists(dset2_file_A)
  if dset2_existsA==False:
    print("WARNING: 'after' file "+dset2_file_A+" does not exist!")
    existA = 1

		# Path on FTP/web archive server (for output NetCDF only)
  str_dset2_path    = dset2_path
  str_dset2_path_B4 = dset2_path_B4
  str_dset2_path_A  = dset2_path_A

  dset2_src  = str_dset2_path_B4+dset2_filename_B4
  dset2_src += ", "+str_dset2_path+dset2_filename
  dset2_src += ", "+str_dset2_path_A+dset2_filename_A

		# Variable names
  dset2_lat_var  = 'latitude'
  dset2_lon_var  = 'longitude'
  dset2_prs_var  = 'pressure'
  dset2_yr_var   = 'year'
  dset2_mm_var   = 'month'
  dset2_dy_var   = 'day'
  dset2_hr_var   = 'hour'
  dset2_mn_var   = 'minutes'
  if qi_choice == 'YES_FC':  
    dset2_pccf_var  = 'percent_confidence_yes_forecast'
  else:
    dset2_pccf_var  = 'percent_confidence_no_forecast'

  dset2_satname_var = 'satellite_name'
  dset2_wcm_var     = 'wind_calculation_method'
  dset2_ham_var     = 'height_assignment_method'
  dset2_spd_var     = 'wind_speed'
  dset2_dir_var     = 'wind_direction'  
		
    	#-------------------------------------------------
  	# Load dataset

		#```````````````````````````````````````````````
  		# CURRENT date
  data_hdl = Dataset(dset2_file)

  tdC_lat = np.asarray( data_hdl.variables[dset2_lat_var] )
  tdC_lon = np.asarray( data_hdl.variables[dset2_lon_var] )
  tdC_prs = np.asarray( data_hdl.variables[dset2_prs_var] )
  tdC_yr  = np.asarray( data_hdl.variables[dset2_yr_var]  )
  tdC_mm  = np.asarray( data_hdl.variables[dset2_mm_var]  )
  tdC_dy  = np.asarray( data_hdl.variables[dset2_dy_var]  )
  tdC_hr  = np.asarray( data_hdl.variables[dset2_hr_var]  )
  tdC_mn  = np.asarray( data_hdl.variables[dset2_mn_var]  )
  tdC_pccf = np.asarray( data_hdl.variables[dset2_pccf_var]  )

  tdC_wcm = np.asarray( data_hdl.variables[dset2_wcm_var]  )
  tdC_ham = np.asarray( data_hdl.variables[dset2_ham_var]  )
  tdC_spd = np.asarray( data_hdl.variables[dset2_spd_var]  )
  tdC_dir = np.asarray( data_hdl.variables[dset2_dir_var]  )

  tdC_satname = np.asarray( data_hdl.variables[dset2_satname_var]  )

  tdC_wcm = np.asarray( data_hdl.variables[dset2_wcm_var]  )
  tdC_ham = np.asarray( data_hdl.variables[dset2_ham_var]  )
  tdC_spd = np.asarray( data_hdl.variables[dset2_spd_var]  )
  tdC_dir = np.asarray( data_hdl.variables[dset2_dir_var]  )

  data_hdl.close()

	# check pressure units and convert to hPa
  if max(tdC_prs) > 10000.:
    tdC_prs = tdC_prs/100.
    
  if dsetflag == "dep":  
    if existB == 0:
    		#```````````````````````````````````````````````
  		# BEFORE date
      data_hdl = Dataset(dset2_file_B4)

      tdB_lat = np.asarray( data_hdl.variables[dset2_lat_var] )
      tdB_lon = np.asarray( data_hdl.variables[dset2_lon_var] )
      tdB_prs = np.asarray( data_hdl.variables[dset2_prs_var] )
      tdB_yr  = np.asarray( data_hdl.variables[dset2_yr_var]  )
      tdB_mm  = np.asarray( data_hdl.variables[dset2_mm_var]  )
      tdB_dy  = np.asarray( data_hdl.variables[dset2_dy_var]  )
      tdB_hr  = np.asarray( data_hdl.variables[dset2_hr_var]  )
      tdB_mn  = np.asarray( data_hdl.variables[dset2_mn_var]  )
      tdB_pccf = np.asarray( data_hdl.variables[dset2_pccf_var]  )

      tdB_satname = np.asarray( data_hdl.variables[dset2_satname_var]  )

      tdB_wcm = np.asarray( data_hdl.variables[dset2_wcm_var]  )
      tdB_ham = np.asarray( data_hdl.variables[dset2_ham_var]  )
      tdB_spd = np.asarray( data_hdl.variables[dset2_spd_var]  )
      tdB_dir = np.asarray( data_hdl.variables[dset2_dir_var]  )

      data_hdl.close()

	# check pressure units and convert to hPa
      if max(tdB_prs) > 10000.:
        tdB_prs = tdB_prs/100.
    
    if existA == 0:
    		#```````````````````````````````````````````````
  		# AFTER date
      data_hdl = Dataset(dset2_file_A)

      tdA_lat = np.asarray( data_hdl.variables[dset2_lat_var] )
      tdA_lon = np.asarray( data_hdl.variables[dset2_lon_var] )
      tdA_prs = np.asarray( data_hdl.variables[dset2_prs_var] )
      tdA_yr  = np.asarray( data_hdl.variables[dset2_yr_var]  )
      tdA_mm  = np.asarray( data_hdl.variables[dset2_mm_var]  )
      tdA_dy  = np.asarray( data_hdl.variables[dset2_dy_var]  )
      tdA_hr  = np.asarray( data_hdl.variables[dset2_hr_var]  )
      tdA_mn  = np.asarray( data_hdl.variables[dset2_mn_var]  )
      tdA_pccf = np.asarray( data_hdl.variables[dset2_pccf_var]  )

      tdA_satname = np.asarray( data_hdl.variables[dset2_satname_var]  )
 
      tdA_wcm = np.asarray( data_hdl.variables[dset2_wcm_var]  )
      tdA_ham = np.asarray( data_hdl.variables[dset2_ham_var]  )
      tdA_spd = np.asarray( data_hdl.variables[dset2_spd_var]  )
      tdA_dir = np.asarray( data_hdl.variables[dset2_dir_var]  )

      data_hdl.close()

	# check pressure units and convert to hPa
      if max(tdA_prs) > 10000.:
        tdA_prs = tdA_prs/100.

	#-------------------------------------------------
      	# Append current arrays to before date arrays

		# get BEFORE and AFTER obs
    if hour=="00":
    	# keep ihrB hrs before CURRENT and ihrA hrs after CURRENT 6hrs
      tiHHA  = np.where(((tdC_hr >= 3) * (tdC_hr < (3+ihrA))))
      siHHA  = np.asarray(tiHHA)
      iHHA   = siHHA.flatten()

      stdA_lat = tdC_lat[iHHA]
      stdA_lon = tdC_lon[iHHA]
      stdA_yr  = tdC_yr [iHHA]
      stdA_mm  = tdC_mm [iHHA]
      stdA_dy  = tdC_dy [iHHA]
      stdA_hr  = tdC_hr [iHHA]
      stdA_mn  = tdC_mn [iHHA]
      stdA_spd = tdC_spd[iHHA]
      stdA_dir = tdC_dir[iHHA]
      stdA_prs = tdC_prs[iHHA]
      stdA_pccf = tdC_pccf[iHHA]
      stdA_satname = tdC_satname[iHHA]
      stdA_wcm = tdC_wcm[iHHA]
      stdA_ham = tdC_ham[iHHA]
	
      if existB == 0:
        tiHHB4 = np.where(((tdB_hr >= (21-ihrB)) * (tdB_hr < 21)))
        siHHB4 = np.asarray(tiHHB4)
        iHHB4  = siHHB4.flatten()
	
        stdB_lat = tdB_lat[iHHB4]
        stdB_lon = tdB_lon[iHHB4]
        stdB_yr  = tdB_yr [iHHB4]
        stdB_mm  = tdB_mm [iHHB4]
        stdB_dy  = tdB_dy [iHHB4]
        stdB_hr  = tdB_hr [iHHB4]
        stdB_mn  = tdB_mn [iHHB4]
        stdB_spd = tdB_spd[iHHB4]
        stdB_dir = tdB_dir[iHHB4]
        stdB_prs = tdB_prs[iHHB4]
        stdB_pccf = tdB_pccf[iHHB4]
        stdB_satname = tdB_satname[iHHB4]
        stdB_wcm = tdB_wcm[iHHB4]
        stdB_ham = tdB_ham[iHHB4]
	
    elif hour=="06":
    	# keep ihrB before CURRENT and ihrA after CURRENT 6hrs
      tiHHA   = np.where(((tdC_hr >= 9) * (tdC_hr < (9+ihrA))))
      siHHA   = np.asarray(tiHHA)
      iHHA    = siHHA.flatten()
      
      stdA_lat = tdC_lat[iHHA]
      stdA_lon = tdC_lon[iHHA]
      stdA_yr  = tdC_yr [iHHA]
      stdA_mm  = tdC_mm [iHHA]
      stdA_dy  = tdC_dy [iHHA]
      stdA_hr  = tdC_hr [iHHA]
      stdA_mn  = tdC_mn [iHHA]
      stdA_spd = tdC_spd[iHHA]
      stdA_dir = tdC_dir[iHHA]
      stdA_prs = tdC_prs[iHHA]
      stdA_pccf = tdC_pccf[iHHA]
      stdA_satname = tdC_satname[iHHA]
      stdA_wcm = tdC_wcm[iHHA]
      stdA_ham = tdC_ham[iHHA]
      
      if (3-ihrB) < 0:
        dhrB = abs(3 - ihrB)
      
        tiHHB4 = np.where(((tdC_hr >= 0) * (tdC_hr < 3)))
        siHHB4 = np.asarray(tiHHB4)
        iHHB4  = siHHB4.flatten()
      
        if existB == 0:
          tiHHB41 = np.where(((tdB_hr >= (24-dhrB)) * (tdB_hr < 24)))  
          siHHB41 = np.asarray(tiHHB41)
          iHHB41  = siHHB41.flatten()
	  
          stdB_lat = np.append(tdB_lat[iHHB41], tdC_lat[iHHB4], axis=0)
          stdB_lon = np.append(tdB_lon[iHHB41], tdC_lon[iHHB4], axis=0)
          stdB_yr  = np.append(tdB_yr [iHHB41], tdC_yr [iHHB4], axis=0)
          stdB_mm  = np.append(tdB_mm [iHHB41], tdC_mm [iHHB4], axis=0)
          stdB_dy  = np.append(tdB_dy [iHHB41], tdC_dy [iHHB4], axis=0)
          stdB_hr  = np.append(tdB_hr [iHHB41], tdC_hr [iHHB4], axis=0)
          stdB_mn  = np.append(tdB_mn [iHHB41], tdC_mn [iHHB4], axis=0)
          stdB_spd = np.append(tdB_spd[iHHB41], tdC_spd[iHHB4], axis=0)
          stdB_dir = np.append(tdB_dir[iHHB41], tdC_dir[iHHB4], axis=0)
          stdB_prs = np.append(tdB_prs[iHHB41], tdC_prs[iHHB4], axis=0)
          stdB_pccf = np.append(tdB_pccf[iHHB41], tdC_pccf[iHHB4], axis=0)
          stdB_satname = np.append(tdB_satname[iHHB41], tdC_satname[iHHB4], axis=0)
          stdB_wcm = np.append(tdB_wcm[iHHB41], tdC_wcm[iHHB4], axis=0)
          stdB_ham = np.append(tdB_ham[iHHB41], tdC_ham[iHHB4], axis=0)
        else:
          stdB_lat = tdC_lat[iHHB4]
          stdB_lon = tdC_lon[iHHB4]
          stdB_yr  = tdC_yr [iHHB4]
          stdB_mm  = tdC_mm [iHHB4]
          stdB_dy  = tdC_dy [iHHB4]
          stdB_hr  = tdC_hr [iHHB4]
          stdB_mn  = tdC_mn [iHHB4]
          stdB_spd = tdC_spd[iHHB4]
          stdB_dir = tdC_dir[iHHB4]
          stdB_prs = tdC_prs[iHHB4]
          stdB_pccf = tdC_pccf[iHHB4]
          stdB_satname = tdC_satname[iHHB4]
          stdB_wcm = tdC_wcm[iHHB4]
          stdB_ham = tdC_ham[iHHB4]

        del dhrB
	
      else:
        tiHHB4 = np.where(((tdC_hr >= (3-ihrB)) * (tdC_hr < 3)))
        siHHB4 = np.asarray(tiHHB4)
        iHHB4  = siHHB4.flatten()
	
        stdB_lat = tdC_lat[iHHB4]
        stdB_lon = tdC_lon[iHHB4]
        stdB_yr  = tdC_yr [iHHB4]
        stdB_mm  = tdC_mm [iHHB4]
        stdB_dy  = tdC_dy [iHHB4]
        stdB_hr  = tdC_hr [iHHB4]
        stdB_mn  = tdC_mn [iHHB4]
        stdB_spd = tdC_spd[iHHB4]
        stdB_dir = tdC_dir[iHHB4]
        stdB_prs = tdC_prs[iHHB4]
        stdB_pccf = tdC_pccf[iHHB4]
        stdB_satname = tdC_satname[iHHB4]
        stdB_wcm = tdC_wcm[iHHB4]
        stdB_ham = tdC_ham[iHHB4]
	
    elif hour=="12":
    	# keep ihrB before CURRENT and ihrA after CURRENT 6hrs
      tiHHA  = np.where(((tdC_hr >= 15) * (tdC_hr < (15+ihrA)))) 
      tiHHB4 = np.where(((tdC_hr >= (9-ihrB)) * (tdC_hr < 9)))
    
      siHHA  = np.asarray(tiHHA)
      iHHA   = siHHA.flatten()
      siHHB4 = np.asarray(tiHHB4)
      iHHB4  = siHHB4.flatten()
      
      stdA_lat = tdC_lat[iHHA]
      stdA_lon = tdC_lon[iHHA]
      stdA_yr  = tdC_yr [iHHA]
      stdA_mm  = tdC_mm [iHHA]
      stdA_dy  = tdC_dy [iHHA]
      stdA_hr  = tdC_hr [iHHA]
      stdA_mn  = tdC_mn [iHHA]
      stdA_spd = tdC_spd[iHHA]
      stdA_dir = tdC_dir[iHHA]
      stdA_prs = tdC_prs[iHHA]
      stdA_pccf = tdC_pccf[iHHA]
      stdA_satname = tdC_satname[iHHA]
      stdA_wcm = tdC_wcm[iHHA]
      stdA_ham = tdC_ham[iHHA]
      
      stdB_lat = tdC_lat[iHHB4]
      stdB_lon = tdC_lon[iHHB4]
      stdB_yr  = tdC_yr [iHHB4]
      stdB_mm  = tdC_mm [iHHB4]
      stdB_dy  = tdC_dy [iHHB4]
      stdB_hr  = tdC_hr [iHHB4]
      stdB_mn  = tdC_mn [iHHB4]
      stdB_spd = tdC_spd[iHHB4]
      stdB_dir = tdC_dir[iHHB4]
      stdB_prs = tdC_prs[iHHB4]
      stdB_pccf = tdC_pccf[iHHB4]
      stdB_satname = tdC_satname[iHHB4]
      stdB_wcm = tdC_wcm[iHHB4]
      stdB_ham = tdC_ham[iHHB4]
    
    elif hour=="18":
    	# keep ihrB before CURRENT and ihrA after CURRENT 6hrs
      if (21+ihrA) >= 24:
        dhrA = abs(24 - (21 + ihrA))

        tiHHA = np.where(((tdC_hr >= 21) * (tdC_hr < 24)))
        siHHA = np.asarray(tiHHA)
        iHHA  = siHHA.flatten()
      
        if existA == 0:
          tiHHA2 = np.where(((tdA_hr >= 0) * (tdA_hr < dhrA)))  
          siHHA2 = np.asarray(tiHHA2)
          iHHA2  = siHHA2.flatten()
	  
          stdA_lat = np.append(tdC_lat[iHHA], tdA_lat[iHHA2], axis=0)
          stdA_lon = np.append(tdC_lon[iHHA], tdA_lon[iHHA2], axis=0)
          stdA_yr  = np.append(tdC_yr [iHHA], tdA_yr [iHHA2], axis=0)
          stdA_mm  = np.append(tdC_mm [iHHA], tdA_mm [iHHA2], axis=0)
          stdA_dy  = np.append(tdC_dy [iHHA], tdA_dy [iHHA2], axis=0)
          stdA_hr  = np.append(tdC_hr [iHHA], tdA_hr [iHHA2], axis=0)
          stdA_mn  = np.append(tdC_mn [iHHA], tdA_mn [iHHA2], axis=0)
          stdA_spd = np.append(tdC_spd[iHHA], tdA_spd[iHHA2], axis=0)
          stdA_dir = np.append(tdC_dir[iHHA], tdA_dir[iHHA2], axis=0)
          stdA_prs = np.append(tdC_prs[iHHA], tdA_prs[iHHA2], axis=0)
          stdA_pccf = np.append(tdC_pccf[iHHA], tdA_pccf[iHHA2], axis=0)
          stdA_satname = np.append(tdC_satname[iHHA], tdA_satname[iHHA2], axis=0)
          stdA_wcm = np.append(tdC_wcm[iHHA], tdA_wcm[iHHA2], axis=0)
          stdA_ham = np.append(tdC_ham[iHHA], tdA_ham[iHHA2], axis=0)
        else:
          stdA_lat = tdC_lat[iHHA]
          stdA_lon = tdC_lon[iHHA]
          stdA_yr  = tdC_yr [iHHA]
          stdA_mm  = tdC_mm [iHHA]
          stdA_dy  = tdC_dy [iHHA]
          stdA_hr  = tdC_hr [iHHA]
          stdA_mn  = tdC_mn [iHHA]
          stdA_spd = tdC_spd[iHHA]
          stdA_dir = tdC_dir[iHHA]
          stdA_prs = tdC_prs[iHHA]
          stdA_pccf = tdC_pccf[iHHA]
          stdA_satname = tdC_satname[iHHA]
          stdA_wcm = tdC_wcm[iHHA]
          stdA_ham = tdC_ham[iHHA]
      
      tiHHB4 = np.where(((tdC_hr >= (15-ihrB)) * (tdC_hr < 15)))
      siHHB4 = np.asarray(tiHHB4)
      iHHB4  = siHHB4.flatten()
    
      stdB_lat = tdC_lat[iHHB4]
      stdB_lon = tdC_lon[iHHB4]
      stdB_yr  = tdC_yr [iHHB4]
      stdB_mm  = tdC_mm [iHHB4]
      stdB_dy  = tdC_dy [iHHB4]
      stdB_hr  = tdC_hr [iHHB4]
      stdB_mn  = tdC_mn [iHHB4]
      stdB_spd = tdC_spd[iHHB4]
      stdB_dir = tdC_dir[iHHB4]
      stdB_prs = tdC_prs[iHHB4]
      stdB_pccf = tdC_pccf[iHHB4]
      stdB_satname = tdC_satname[iHHB4]
      stdB_wcm = tdC_wcm[iHHB4]
      stdB_ham = tdC_ham[iHHB4]

      		# append CURRENT to BEFORE
    if existB == 0:
      tdset2_lat  = np.append(stdB_lat ,tdC_lat ,axis=0)
      tdset2_lon  = np.append(stdB_lon ,tdC_lon ,axis=0)
      tdset2_prs  = np.append(stdB_prs ,tdC_prs ,axis=0)
      tdset2_yr   = np.append(stdB_yr  ,tdC_yr  ,axis=0)
      tdset2_mm   = np.append(stdB_mm  ,tdC_mm  ,axis=0)
      tdset2_dy   = np.append(stdB_dy  ,tdC_dy  ,axis=0)
      tdset2_hr   = np.append(stdB_hr  ,tdC_hr  ,axis=0)
      tdset2_mn   = np.append(stdB_mn  ,tdC_mn  ,axis=0)
      tdset2_pccf = np.append(stdB_pccf,tdC_pccf,axis=0)
      tdset2_satname = np.append(stdB_satname,tdC_satname,axis=0)
      tdset2_wcm = np.append(stdB_wcm,tdC_wcm,axis=0)
      tdset2_ham = np.append(stdB_ham,tdC_ham,axis=0)
      tdset2_spd = np.append(stdB_spd,tdC_spd,axis=0)
      tdset2_dir = np.append(stdB_dir,tdC_dir,axis=0)
    else:
      tdset2_lat  = tdC_lat 
      tdset2_lon  = tdC_lon 
      tdset2_prs  = tdC_prs 
      tdset2_yr   = tdC_yr  
      tdset2_mm   = tdC_mm  
      tdset2_dy   = tdC_dy  
      tdset2_hr   = tdC_hr  
      tdset2_mn   = tdC_mn  
      tdset2_pccf = tdC_pccf
      tdset2_satname = tdC_satname
      tdset2_wcm = tdC_wcm
      tdset2_ham = tdC_ham
      tdset2_spd = tdC_spd
      tdset2_dir = tdC_dir

  		# append AFTER to CURRENT
    if existA == 0:
      tdset2_lat  = np.append(tdset2_lat ,stdA_lat ,axis=0)
      tdset2_lon  = np.append(tdset2_lon ,stdA_lon ,axis=0)
      tdset2_prs  = np.append(tdset2_prs ,stdA_prs ,axis=0)
      tdset2_yr   = np.append(tdset2_yr  ,stdA_yr  ,axis=0)
      tdset2_mm   = np.append(tdset2_mm  ,stdA_mm  ,axis=0)
      tdset2_dy   = np.append(tdset2_dy  ,stdA_dy  ,axis=0)
      tdset2_hr   = np.append(tdset2_hr  ,stdA_hr  ,axis=0)
      tdset2_mn   = np.append(tdset2_mn  ,stdA_mn  ,axis=0)
      tdset2_pccf = np.append(tdset2_pccf,stdA_pccf,axis=0)
      tdset2_satname = np.append(tdset2_satname,stdA_satname,axis=0)
      tdset2_wcm = np.append(tdset2_wcm,stdA_wcm,axis=0)
      tdset2_ham = np.append(tdset2_ham,stdA_ham,axis=0)
      tdset2_spd = np.append(tdset2_spd,stdA_spd,axis=0)
      tdset2_dir = np.append(tdset2_dir,stdA_dir,axis=0)

  else:
  	# dsetflag = "drv"
    tdset2_lat  = tdC_lat 
    tdset2_lon  = tdC_lon 
    tdset2_prs  = tdC_prs 
    tdset2_yr   = tdC_yr  
    tdset2_mm   = tdC_mm  
    tdset2_dy   = tdC_dy  
    tdset2_hr   = tdC_hr  
    tdset2_mn   = tdC_mn  
    tdset2_pccf = tdC_pccf
    tdset2_satname = tdC_satname
    tdset2_wcm = tdC_wcm
    tdset2_ham = tdC_ham
    tdset2_spd = tdC_spd
    tdset2_dir = tdC_dir

	#----------------------------------------
	# if 'runtype' = 'match', get indices of matches and apply QC if bool_drv_qc=True.
	# if 'runtype' = 'plot', SKIP bool if-block
	
  if runtype == "match":
    if bool_qc:
      tindexes2,qc_list = qc_amv(tdset2_pccf,pct)
    
      sindexes2 = np.asarray(tindexes2)
      indexes2  = sindexes2.flatten()
 
      d_lat = tdset2_lat[indexes2]
      d_lon = tdset2_lon[indexes2]
      d_prs = tdset2_prs[indexes2]
      d_yr  = tdset2_yr [indexes2]
      d_mm  = tdset2_mm [indexes2]
      d_dy  = tdset2_dy [indexes2]
      d_hr  = tdset2_hr [indexes2]
      d_mn  = tdset2_mn [indexes2]
    
      d_satname = tdset2_satname[indexes2]
      d_wcm = tdset2_wcm[indexes2]
      d_ham = tdset2_ham[indexes2]
      d_spd = tdset2_spd[indexes2]
      d_dir = tdset2_dir[indexes2]

    elif not bool_qc:
      # Do not apply QC to AMV_NCEP dataset

      sindexesDC = np.asarray(np.where(tdset2_lat==tdset2_lat))           #get all indices
      indexes2   = sindexesDC.flatten()

      qc_list = "No QC applied"

      d_lat = tdset2_lat 
      d_lon = tdset2_lon
      d_yr  = tdset2_yr 
      d_mm  = tdset2_mm 
      d_dy  = tdset2_dy 
      d_hr  = tdset2_hr 
      d_mn  = tdset2_mn  
      d_prs = tdset2_prs
    
      d_satname = tdset2_satname
      d_wcm = tdset2_wcm
      d_ham = tdset2_ham
      d_spd = tdset2_spd
      d_dir = tdset2_dir

  elif runtype == "plot":

    indexes2 = [0]

    qc_list = "No QC applied"

    d_lat = tdset2_lat[idxs]
    d_lon = tdset2_lon[idxs]
    d_yr  = tdset2_yr [idxs]
    d_mm  = tdset2_mm [idxs]
    d_dy  = tdset2_dy [idxs]
    d_hr  = tdset2_hr [idxs]
    d_mn  = tdset2_mn [idxs]
    d_prs = tdset2_prs[idxs]
    
    d_satname = tdset2_satname[idxs]
    d_wcm = tdset2_wcm[idxs]
    d_ham = tdset2_ham[idxs]
    d_spd = tdset2_spd[idxs]
    d_dir = tdset2_dir[idxs]

	#----------------------------------------

	# create height array but fill with missing -999.
	#	AMV heights not available.
  d_hgt = np.nan * np.ones_like(d_lat)

	# Return variables to MAIN
  return d_lat,d_lon,d_yr,d_mm,d_dy,d_hr,d_mn,d_hgt,d_prs,indexes2,qc_list,dset2_src, d_satname,d_wcm,d_ham,d_spd,d_dir

     
#===========================================================================================
# Read DAWN data
#
#       INPUTS:
#               yyyymmdd .......................... Current date in yyyymmdd format
#               bool_qc ......................... Choice to apply DAWN QC: True=apply QC, False=don't apply QC
#
#       OUTPUTS:
#               d_lat ............................... Latitude in degrees [-90,90]
#               d_lon ............................... Longitude in degrees [0,360]
#               d_prs ............................... Pressure in hPa
#               d_hgt ............................... Height in km
#               d_yr ................................ Year
#               d_mm ................................ Month
#               d_dy ................................ Day
#               d_hr ................................ Hour
#               d_mn ................................ Minute
#               qc_list ............................. List of QC applied (if applicable)
#

def read_DAWN(path_prefix,yyyymmddhh,dateB4,dateA,bool_qc,dsetflag,runtype,time_diff_max,idx):

  print("read_DAWN: bool_qc = "+str(bool_qc))

  qc_list = ""                  #initialize

  yyyy = yyyymmddhh[0:4]
  mm   = yyyymmddhh[4:6]
  dd   = yyyymmddhh[6:8]
  hour = yyyymmddhh[8:10]

  yyyymmdd = yyyy+mm+dd

  yyB4 = dateB4[0:4]
  mmB4 = dateB4[4:6]
  ddB4 = dateB4[6:8]

  yyA  = dateA[0:4]
  mmA  = dateA[4:6]
  ddA  = dateA[6:8]
 

               # convert 'time_max_diff' to rounded integer
  tqc_time = float(time_diff_max)
  if tqc_time%60 != 0.0:
    qc_time = int(math.ceil(tqc_time/60.0))             # round up to nearest or equal integer hour
  else:
    qc_time = int(tqc_time)

                # convert hour to integer
  ihr = int(hour)

                # find integer hours for before/after files for collocation
  if hour == "00":
        # 00
    ihrB = 24 - qc_time
    ihrA = ihr + qc_time
  else:
        # 06,12,18
    ihrB = ihr - qc_time
    ihrA = ihr + qc_time

  del qc_time,tqc_time,ihr

        #-------------------------------------------------
        # Define dataset
  tmp_dset1_path    = '/home/jlocke/wind_datasets/DAWN/CPEX-CV/'+yyyy+'/'+mm+'/'+dd+'/'
  tmp_dset1_path_B4 = '/home/jlocke/wind_datasets/DAWN/CPEX-CV/'+yyB4+'/'+mmB4+'/'+ddB4+'/'
  tmp_dset1_path_A  = '/home/jlocke/wind_datasets/DAWN/CPEX-CV/'+yyA+'/'+mmA+'/'+ddA+'/'
 
                # Path/file
                # ...Current date
  dset1_path            = path_prefix+tmp_dset1_path
  dset1_filename        = 'cpexcv-DAWN_DC8_'+yyyymmdd+'_R0.nc'
  dset1_file            = dset1_path+dset1_filename
                # ...Before date
  dset1_path_B4         = path_prefix+tmp_dset1_path_B4
  dset1_filename_B4     = 'cpexcv-DAWN_DC8_'+dateB4+'_R0.nc'
  dset1_file_B4         = dset1_path_B4+dset1_filename_B4
                # ...After date
  dset1_path_A          = path_prefix+tmp_dset1_path_A
  dset1_filename_A      = 'cpexcv-DAWN_DC8_'+dateA+'_R0.nc'
  dset1_file_A          = dset1_path_A+dset1_filename_A

                        # initialize flag indicating if dataset exists. 0=yes, 1=no
  existB = 0            # ... B = date Before current
  existA = 0            # ... A = date After current

  dset1_exists = exists(dset1_file)
  if dset1_exists==False:
    print("ERROR: file "+dset1_file+" does not exist!")
    sys.exit()                                                  #exit script immediately
  dset1_existsB = exists(dset1_file_B4)
  if dset1_existsB==False:
    print("WARNING: 'before' file "+dset1_file_B4+" does not exist!")
    existB = 1
  dset1_existsA = exists(dset1_file_A)
  if dset1_existsA==False:
    print("WARNING: 'after' file "+dset1_file_A+" does not exist!")
    existA = 1

                # Path on FTP/web archive server (for output NetCDF only)
  str_dset1_path    = dset1_path
  str_dset1_path_B4 = dset1_path_B4
  str_dset1_path_A  = dset1_path_A

  dset1_src  = str_dset1_path_B4+dset1_filename_B4
  dset1_src += ", "+str_dset1_path+dset1_filename
  dset1_src += ", "+str_dset1_path_A+dset1_filename_A

                # Variable names
  dset1_lat_var = 'lat'
  dset1_lon_var = 'lon'
  dset1_t_var = 'time' ## seconds since yyyymmdd 00:00:00
  dset1_hgt_var = 'z'
  dset1_spd_var = 'Wind_Speed'
  dset1_dir_var = 'Wind_Direction'
  dset1_LOSspd_var = 'Line_of_Sight_Wind_Speed'

        #-------------------------------------------------
        # Load dataset

                #```````````````````````````````````````````````
                # CURRENT date

  data_hdl = Dataset(dset1_file)
  
  tdC_lat = np.asarray( data_hdl.variables[dset1_lat_var] )
  tdC_lon = np.asarray( data_hdl.variables[dset1_lon_var] )
  tdC_t = np.asarray( data_hdl.variables[dset1_t_var] )
  tdC_hgt = np.asarray( data_hdl.variables[dset1_hgt_var]  )
  tdC_spdn = np.asarray( data_hdl.variables[dset1_spd_var]  )
  tdC_dirn = np.asarray( data_hdl.variables[dset1_dir_var]  )
  tdC_LOSspd = np.asarray( data_hdl.variables[dset1_LOSspd_var]  )
  
  data_hdl.close()
  
  tdC_spd = tdC_spdn[1]
  tdC_dir = tdC_dirn[1]

###### CONVERT TIME ##### time is written in seconds after 00:00:00 on start date
  tdC_hr = []
  htm = []
  td_hr = []
  for i in tdC_t:              
   hr = (i/3600)             ### Divide seconds by 3600 to get hours not rounded
   hr_r = math.ceil(i/3600)  ### Hours rounded
   left_hr = (hr_r - hr)     ### Whole hours - rounded hours to get leftover -> caluclate min
   td_hr.append(hr_r)
   htm.append(left_hr)  

  for x in td_hr:  
   tdC_hr.append(x-12)
  tdC_hr = np.array(tdC_hr)

  tdC_mn = []
  for k in htm: 
   tdC_m = math.ceil(k*60)  ### using leftover hours * 60 to get min
   tdC_mn.append(tdC_m)
  tdC_mn = np.array(tdC_mn) 
   
###### Create Year, Day, and Months arrays
  tdC_yr = []
  tdC_dy = []
  tdC_mm = []
 
  for n in tdC_t: 
   tdC_y = yyyy
   tdC_d = dd
   tdC_m = mm 
   tdC_yr.append(tdC_y)
   tdC_mm.append(tdC_m)
   tdC_dy.append(tdC_d)
  
  tdC_yr = np.array(tdC_yr)
  tdC_mm = np.array(tdC_mm)
  tdC_dy = np.array(tdC_dy) 
   
##### Create pressure variable 
#  tdC_prs = [] 
#  for m in tdC_hgt:
#   Pr = 101325*(1-2.25577*10**-5*(abs(m)))**5.25588
#   tdC_prs.append(Pr)
#  tdC_prs = np.array(tdC_prs, dtype =np.int) 
 
##### check pressure units and convert to hPa
#  if max(tdC_prs) > 10000.:
#    tdC_prs = tdC_prs/100.
 
  if dsetflag == "dep":
    if existB == 0:
                #```````````````````````````````````````````````
                # BEFORE date
      data_hdl = Dataset(dset1_file_B4)

      tdB_lat = np.asarray( data_hdl.variables[dset1_lat_var] )
      tdB_lon = np.asarray( data_hdl.variables[dset1_lon_var] )
      tdB_t = np.asarray( data_hdl.variables[dset1_t_var] )
      tdB_hgt = np.asarray( data_hdl.variables[dset1_hgt_var]  )
      tdB_spdn = np.asarray( data_hdl.variables[dset1_spd_var]  )
      tdB_dirn= np.asarray( data_hdl.variables[dset1_dir_var]  )
      tdB_LOSspd = np.asarray( data_hdl.variables[dset1_LOSspd_var]  )

      data_hdl.close()

      tdB_spd = tdB_spdn[1]
      tdB_dir = tdB_dirn[1]
###### CONVERT TIME ##### time is written in seconds after 00:00:00 on start date
      tdB_hr = []
      htmB = []
      td_hrB = []
      for i in tdB_t:
       hrB = (i/3600)             ### Divide seconds by 3600 to get hours not rounded
       hr_rB = math.ceil(i/3600)  ### Hours rounded
       left_hrB = (hr_rB - hrB)     ### Whole hours - rounded hours to get leftover -> caluclate min
       td_hrB.append(hr_rB)
       htmB.append(left_hrB)
 
      for x in td_hrB:
       tdB_hr.append(x-12)
      tdB_hr = np.array(tdB_hr)

      tdB_mn = []
      for k in htmB:
       tdB_m = math.ceil(k*60)  ### using leftover hours * 60 to get min
       tdB_mn.append(tdB_m)
      tdB_mn = np.array(tdB_mn)

###### Create Year, Day, and Months arrays
      tdB_yr = []
      tdB_dy = []
      tdB_mm = []

      for n in tdB_t:
       tdB_y = yyB4
       tdB_d = ddB4
       tdB_m = mmB4
       tdB_yr.append(tdB_y)
       tdB_mm.append(tdB_m)
       tdB_dy.append(tdB_d)

      tdB_yr = np.array(tdB_yr)
      tdB_mm = np.array(tdB_mm)
      tdB_dy = np.array(tdB_dy)

##### Create pressure variable
#      tdB_prs = []
#      for m in tdB_hgt:
#       PrB = 101325*(1-2.25577*10**-5*(abs(m)))**5.25588
#       tdB_prs.append(PrB)
#      tdB_prs = np.array(tdB_prs, dtype =np.int)
#
##### check pressure units and convert to hPa
#      if max(tdB_prs) > 10000.:
#       tdB_prs = tdB_prs/100.
  
    if existA == 0:

                #```````````````````````````````````````````````
                # AFTER date
      data_hdl = Dataset(dset1_file_A)

 
      tdA_lat = np.asarray( data_hdl.variables[dset1_lat_var] )
      tdA_lon = np.asarray( data_hdl.variables[dset1_lon_var] )
      tdA_t = np.asarray( data_hdl.variables[dset1_t_var] )
      tdA_hgt = np.asarray( data_hdl.variables[dset1_hgt_var]  )
      tdA_spdn = np.asarray( data_hdl.variables[dset1_spd_var]  )
      tdA_dirn = np.asarray( data_hdl.variables[dset1_dir_var]  )
      tdA_LOSspd = np.asarray( data_hdl.variables[dset1_LOSspd_var]  )

      data_hdl.close()
    
      tdA_spd = tdA_spdn[1]
      tdA_dir = tdA_dirn[1]  
###### CONVERT TIME ##### time is written in seconds after 00:00:00 on start date
      tdA_hr = []
      htmA = []
      td_hrA = []
      for i in tdA_t:
       hrA = (i/3600)             ### Divide seconds by 3600 to get hours not rounded
       hr_rA = math.ceil(i/3600)  ### Hours rounded
       left_hrA = (hr_rA - hrA)     ### Whole hours - rounded hours to get leftover -> caluclate min
       td_hrA.append(hr_rA)
       htmA.append(left_hrA)
     

      for x in td_hrA:
       tdA_hr.append(x-12)
      tdA_hr = np.array(tdA_hr)
  
      tdA_mn = []
      for k in htmA:
       tdA_m = math.ceil(k*60)  ### using leftover hours * 60 to get min
       tdA_mn.append(tdA_m)
      tdA_mn = np.array(tdA_mn)

###### Create Year, Day, and Months arrays
      tdA_yr = []
      tdA_dy = []
      tdA_mm = []

      for n in tdA_t:
       tdA_y = yyA
       tdA_d = ddA
       tdA_m = mmA
       tdA_yr.append(tdA_y)
       tdA_mm.append(tdA_m)
       tdA_dy.append(tdA_d)

      tdA_yr = np.array(tdA_yr)
      tdA_mm = np.array(tdA_mm)
      tdA_dy = np.array(tdA_dy)
      
##### Create pressure variable
#      tdA_prs = []
#      #tdA_hgt.resize((0,g))
#      for m in tdA_hgt:
#       PrA = 101325*(1-2.25577*10**-5*(abs(m)))**5.25588
#       tdA_prs.append(PrA)
#      tdA_prs = np.array(tdA_prs, dtype =np.int)
   
##### check pressure units and convert to hPa
#      if max(tdA_prs) > 10000.:
#       tdA_prs = tdA_prs/100.
  
        #-------------------------------------------------
        # Append current arrays to before date arrays

                # get BEFORE and AFTER obs
    if hour=="00":
        # keep ihrB hrs before CURRENT and ihrA hrs after CURRENT 6hrs
      tiHHA  = np.where(((tdC_hr >= 3) * (tdC_hr < (3+ihrA))))
      siHHA  = np.asarray(tiHHA)
      iHHA   = siHHA.flatten()
      iHHA.resize((434))

      ihhA = np.copy(iHHA)
      ihhA = np.array(ihhA)
      ihhA.resize((401))

      stdA_lat = tdC_lat[iHHA]
      stdA_lon = tdC_lon[iHHA]
      stdA_yr  = tdC_yr [iHHA]
      stdA_mm  = tdC_mm [iHHA]
      stdA_dy  = tdC_dy [iHHA]
      stdA_hr  = tdC_hr [iHHA]
      stdA_mn  = tdC_mn [iHHA]
      stdA_spd = tdC_spd[ihhA]
      stdA_dir = tdC_dir[ihhA]
#      stdA_prs = tdC_prs[ihhA]
      stdA_hgt = tdC_hgt[ihhA]
#      stdA_azm = tdC_azm[iHHA]
#      stdA_elv = tdC_elv[iHHA]
 
      if existB == 0:
        tiHHB4 = np.where(((tdB_hr >= (21-ihrB)) * (tdB_hr < 21)))
        siHHB4 = np.asarray(tiHHB4)
        iHHB4  = siHHB4.flatten()
        iHHB4.resize((434))
         
        ihhB4 = np.copy(iHHB4)
        ihhB4 = np.array(ihhB4)
        ihhB4.resize((401))
     
 
        stdB_lat = tdB_lat[iHHB4]
        stdB_lon = tdB_lon[iHHB4]
        stdB_yr  = tdB_yr [iHHB4]
        stdB_mm  = tdB_mm [iHHB4]
        stdB_dy  = tdB_dy [iHHB4]
        stdB_hr  = tdB_hr [iHHB4]
        stdB_mn  = tdB_mn [iHHB4]
        stdB_spd = tdB_spd[ihhB4]
        stdB_dir = tdB_dir[iHHB4]
#        stdB_prs = tdB_prs[ihhB4]
        stdB_hgt = tdB_hgt[ihhB4]
#        stdB_azm = tdB_azm[iHHB4]
#        stdB_elv = tdB_elv[iHHB4]
       
    elif hour=="06":
        # keep ihrB before CURRENT and ihrA after CURRENT 6hrs
      tiHHA   = np.where(((tdC_hr >= 9) * (tdC_hr < (9+ihrA))))
      siHHA   = np.asarray(tiHHA)
      iHHA    = siHHA.flatten()
      iHHA.resize((434))

      ihhA = np.copy(iHHA)
      ihhA = np.array(ihhA)
      ihhA.resize((401))

      stdA_lat = tdC_lat[iHHA]
      stdA_lon = tdC_lon[iHHA]
      stdA_yr  = tdC_yr [iHHA]
      stdA_mm  = tdC_mm [iHHA]
      stdA_dy  = tdC_dy [iHHA]
      stdA_hr  = tdC_hr [iHHA]
      stdA_mn  = tdC_mn [iHHA]
      stdA_spd = tdC_spd[ihhA]
      stdA_dir = tdC_dir[iHHA]
#      stdA_prs = tdC_prs[ihhA]
      stdA_hgt = tdC_hgt[ihhA]
#      stdA_azm = tdC_azm[iHHA]
#      stdA_elv = tdC_elv[iHHA]

      if (3-ihrB) < 0:
        dhrB = abs(3 - ihrB)

        tiHHB4 = np.where(((tdC_hr >= 0) * (tdC_hr < 3)))
        siHHB4 = np.asarray(tiHHB4)
        iHHB4  = siHHB4.flatten()
        iHHB4.resize((434))

        ihhB4 = np.copy(iHHB4)
        ihhB4 = np.array(ihhB4)
        ihhB4.resize((401))

        if existB == 0:
          tiHHB41 = np.where(((tdB_hr >= (24-dhrB)) * (tdB_hr < 24)))
          siHHB41 = np.asarray(tiHHB41)
          iHHB41  = siHHB41.flatten()
          iHHB41.resize((434))
          
          ihhB41 = np.copy(iHHB41)
          ihhB41 = np.array(ihhB41)
          ihhB41.resize((401))

          stdB_lat = np.append(tdB_lat[iHHB41], tdC_lat[iHHB4], axis=0)
          stdB_lon = np.append(tdB_lon[iHHB41], tdC_lon[iHHB4], axis=0)
          stdB_yr  = np.append(tdB_yr [iHHB41], tdC_yr [iHHB4], axis=0)
          stdB_mm  = np.append(tdB_mm [iHHB41], tdC_mm [iHHB4], axis=0)
          stdB_dy  = np.append(tdB_dy [iHHB41], tdC_dy [iHHB4], axis=0)
          stdB_hr  = np.append(tdB_hr [iHHB41], tdC_hr [iHHB4], axis=0)
          stdB_mn  = np.append(tdB_mn [iHHB41], tdC_mn [iHHB4], axis=0)
          stdB_spd = np.append(tdB_spd[ihhB41], tdC_spd[ihhB4], axis=0)
          stdB_dir = np.append(tdB_dir[iHHB41], tdC_dir[iHHB4], axis=0)
#          stdB_prs = np.append(tdB_prs[ihhB41], tdC_prs[ihhB4], axis=0)
          stdB_hgt = np.append(tdB_hgt[ihhB41], tdC_hgt[ihhB4], axis=0)
#         stdB_azm = np.append(tdB_azm[iHHB41], tdC_azm[iHHB4], axis=0)
#         stdB_elv = np.append(tdB_elv[iHHB41], tdC_elv[iHHB4], axis=0)
        else:
          stdB_lat = tdC_lat[iHHB4]
          stdB_lon = tdC_lon[iHHB4]
          stdB_yr  = tdC_yr [iHHB4]
          stdB_mm  = tdC_mm [iHHB4]
          stdB_dy  = tdC_dy [iHHB4]
          stdB_hr  = tdC_hr [iHHB4]
          stdB_mn  = tdC_mn [iHHB4]
          stdB_spd = tdC_spd[ihhB4]
          stdB_dir = tdC_dir[iHHB4]
#          stdB_prs = tdC_prs[ihhB4]
          stdB_hgt = tdC_hgt[ihhB4]
#          stdB_azm = tdC_azm[iHHB4]
#          stdB_elv = tdC_elv[iHHB4]

        del dhrB

      else:
        tiHHB4 = np.where(((tdC_hr >= (3-ihrB)) * (tdC_hr < 3)))
        siHHB4 = np.asarray(tiHHB4)
        iHHB4  = siHHB4.flatten()
        iHHB4.resize((434))
      
        ihhB4 = np.copy(iHHB4)
        ihhB4 = np.array(ihhB4)
        ihhB4.resize((401))

        stdB_lat = tdC_lat[iHHB4]
        stdB_lon = tdC_lon[iHHB4]
        stdB_yr  = tdC_yr [iHHB4]
        stdB_mm  = tdC_mm [iHHB4]
        stdB_dy  = tdC_dy [iHHB4]
        stdB_hr  = tdC_hr [iHHB4]
        stdB_mn  = tdC_mn [iHHB4]
        stdB_spd = tdC_spd[ihhB4]
        stdB_dir = tdC_dir[iHHB4]
 #       stdB_prs = tdC_prs[ihhB4]
        stdB_hgt = tdC_hgt[ihhB4]
#        stdB_azm = tdC_azm[iHHB4]
#        stdB_elv = tdC_elv[iHHB4]

    elif hour=="12":
        # keep ihrB before CURRENT and ihrA after CURRENT 6hrs
      tiHHA  = np.where(((tdC_hr >= 15) * (tdC_hr < (15+ihrA))))
      tiHHB4 = np.where(((tdC_hr >= (9-ihrB)) * (tdC_hr < 9)))
      siHHA  = np.asarray(tiHHA)
      iHHA   = siHHA.flatten()
      siHHB4 = np.asarray(tiHHB4)
      iHHB4  = siHHB4.flatten()

      iHHB4.resize((434))
      iHHA.resize((434))

      ihhB4 = np.copy(iHHB4)
      ihhB4 = np.array(ihhB4)
      ihhB4.resize((401))

      ihhA = np.copy(iHHA)
      ihhA = np.array(ihhA)
      ihhA.resize((401))
    
      stdA_lat = tdC_lat[iHHA]
      stdA_lon = tdC_lon[iHHA]
      stdA_yr  = tdC_yr [iHHA]
      stdA_mm  = tdC_mm [iHHA]
      stdA_dy  = tdC_dy [iHHA]
      stdA_hr  = tdC_hr [iHHA]
      stdA_mn  = tdC_mn [iHHA]
      stdA_spd = tdC_spd[ihhA]
      stdA_dir = tdC_dir[iHHA]
#      stdA_prs = tdC_prs[ihhA]
      stdA_hgt = tdC_hgt[ihhA]
#      stdA_azm = tdC_azm[iHHA]
#      stdA_elv = tdC_elv[iHHA]

      stdB_lat = tdC_lat[iHHB4]
      stdB_lon = tdC_lon[iHHB4]
      stdB_yr  = tdC_yr [iHHB4]
      stdB_mm  = tdC_mm [iHHB4]
      stdB_dy  = tdC_dy [iHHB4]
      stdB_hr  = tdC_hr [iHHB4]
      stdB_mn  = tdC_mn [iHHB4]
      stdB_spd = tdC_spd[ihhB4]
      stdB_dir = tdC_dir[iHHB4]
#      stdB_prs = tdC_prs[ihhB4]
      stdB_hgt = tdC_hgt[ihhB4]
#      stdB_azm = tdC_azm[iHHB4]
#      stdB_elv = tdC_elv[iHHB4]
    elif hour=="18":
        # keep ihrB before CURRENT and ihrA after CURRENT 6hrs
      if (21+ihrA) >= 24:
        dhrA = abs(24 - (21 + ihrA))
 
        tiHHA = np.where(((tdC_hr >= 21) * (tdC_hr < 24)))
        siHHA = np.asarray(tiHHA)
        iHHA  = siHHA.flatten()
        iHHA.resize((434))
 
        ihhA = np.copy(iHHA)
        ihhA = np.array(ihhA)
        ihhA.resize((401))

        if existA == 0:
         tiHHA2 = np.where(((tdA_hr >= 0) * (tdA_hr < dhrA)))
         siHHA2 = np.asarray(tiHHA2)
         iHHA2  = siHHA2.flatten()
         iHHA2.resize((434))
  
         ihhA2 = np.copy(iHHA2)
         ihhA2 = np.array(ihhA2)
         ihhA2.resize((401))
  
         stdA_lat = np.append(tdC_lat[iHHA], tdA_lat[iHHA2], axis=0) 
         stdA_lon = np.append(tdC_lon[iHHA], tdA_lon[iHHA2], axis=0)
         stdA_yr  = np.append(tdC_yr [iHHA], tdA_yr [iHHA2], axis=0)
         stdA_mm  = np.append(tdC_mm [iHHA], tdA_mm [iHHA2], axis=0)
         stdA_dy  = np.append(tdC_dy [iHHA], tdA_dy [iHHA2], axis=0)
         stdA_hr  = np.append(tdC_hr [iHHA], tdA_hr [iHHA2], axis=0)
         stdA_mn  = np.append(tdC_mn [iHHA], tdA_mn [iHHA2], axis=0)
         stdA_spd = np.append(tdC_spd[ihhA], tdA_spd[ihhA2], axis=0)
         stdA_dir = np.append(tdC_dir[iHHA], tdA_dir[iHHA2], axis=0)
 #        stdA_prs = np.append(tdC_prs[ihhA], tdA_prs[ihhA2], axis=0)
         stdA_hgt = np.append(tdC_hgt[ihhA], tdA_hgt[ihhA2], axis=0)
#        stdA_azm = np.append(tdC_azm[iHHA], tdA_azm[iHHA2], axis=0)
#        stdA_elv = np.append(tdC_elv[iHHA], tdA_elv[iHHA2], axis=0)
      else:
          stdA_lat = tdC_lat[iHHA]
          stdA_lon = tdC_lon[iHHA] 
          stdA_yr  = tdC_yr [iHHA]
          stdA_mm  = tdC_mm [iHHA]
          stdA_dy  = tdC_dy [iHHA]
          stdA_hr  = tdC_hr [iHHA]
          stdA_mn  = tdC_mn [iHHA]
          stdA_spd = tdC_spd[ihhA]
          stdA_dir = tdC_dir[iHHA]
#          stdA_prs = tdC_prs[ihhA]
          stdA_hgt = tdC_hgt[ihhA]
#         stdA_azm = tdC_azm[iHHA]
#         stdA_elv = tdC_elv[iHHA]

      tiHHB4 = np.where(((tdC_hr >= (15-ihrB)) * (tdC_hr < 15)))   
      siHHB4 = np.asarray(tiHHB4)
      iHHB4  = siHHB4.flatten()
      iHHB4.resize((434))

      ihhB4 = np.copy(iHHB4)
      ihhB4 = np.array(ihhB4)
      ihhB4.resize((401))

      stdB_lat = tdC_lat[iHHB4]
      stdB_lon = tdC_lon[iHHB4]
      stdB_yr  = tdC_yr [iHHB4]
      stdB_mm  = tdC_mm [iHHB4]
      stdB_dy  = tdC_dy [iHHB4]
      stdB_hr  = tdC_hr [iHHB4]
      stdB_mn  = tdC_mn [iHHB4]
      stdB_spd = tdC_spd[ihhB4]
      stdB_dir = tdC_dir[iHHB4]
#     stdB_prs = tdC_prs[ihhB4]
      stdB_hgt = tdC_hgt[ihhB4]
#     stdB_azm = tdC_azm[iHHB4]
#     stdB_elv = tdC_elv[iHHB4]
  
              # append CURRENT to BEFORE
    if existB == 0:
      tdset1_lat = np.append(stdB_lat,tdC_lat,axis=0)
      tdset1_lon = np.append(stdB_lon,tdC_lon,axis=0)
      tdset1_yr  = np.append(stdB_yr ,tdC_yr ,axis=0)
      tdset1_mm  = np.append(stdB_mm ,tdC_mm ,axis=0)
      tdset1_dy  = np.append(stdB_dy ,tdC_dy ,axis=0)
      tdset1_hr  = np.append(stdB_hr ,tdC_hr ,axis=0)
      tdset1_mn  = np.append(stdB_mn ,tdC_mn ,axis=0)
      tdset1_hgt = np.append(stdB_hgt,tdC_hgt,axis=0)
#      tdset1_prs = np.append(stdB_prs,tdC_prs,axis=0)
#      tdset1_azm = np.append(stdB_azm,tdC_azm,axis=0)
#      tdset1_elv = np.append(stdB_elv,tdC_elv,axis=0)
      tdset1_spd = np.append(stdB_spd,tdC_spd,axis=0)
      tdset1_dir = np.append(stdB_dir,tdC_dir,axis=0)
    else:
      tdset1_lat = tdC_lat
      tdset1_lon = tdC_lon
      tdset1_yr  = tdC_yr
      tdset1_mm  = tdC_mm
      tdset1_dy  = tdC_dy
      tdset1_hr  = tdC_hr
      tdset1_mn  = tdC_mn
      tdset1_hgt = tdC_hgt
#      tdset1_prs = tdC_prs
#      tdset1_azm = tdC_azm
#      tdset1_elv = tdC_elv
      tdset1_spd = tdC_spd
      tdset1_dir = tdC_dir

                # append AFTER to CURRENT
    if existA == 0:
      tdset1_lat = np.append(tdset1_lat,stdA_lat,axis=0)
      tdset1_lon = np.append(tdset1_lon,stdA_lon,axis=0)
      tdset1_yr  = np.append(tdset1_yr ,stdA_yr ,axis=0)
      tdset1_mm  = np.append(tdset1_mm ,stdA_mm ,axis=0)
      tdset1_dy  = np.append(tdset1_dy ,stdA_dy ,axis=0)
      tdset1_hr  = np.append(tdset1_hr ,stdA_hr ,axis=0)
      tdset1_mn  = np.append(tdset1_mn ,stdA_mn ,axis=0)
      tdset1_hgt = np.append(tdset1_hgt,stdA_hgt,axis=0)
 #     tdset1_prs = np.append(tdset1_prs,stdA_prs,axis=0)
#      tdset1_azm = np.append(tdset1_azm,stdA_azm,axis=0)
#      tdset1_elv = np.append(tdset1_elv,stdA_elv,axis=0)
      tdset1_spd = np.append(tdset1_spd,stdA_spd,axis=0)
      tdset1_dir = np.append(tdset1_dir,stdA_dir,axis=0)

    else:
        # dsetflag = "drv"
      tdset1_lat = tdC_lat
      tdset1_lon = tdC_lon
      tdset1_yr  = tdC_yr
      tdset1_mm  = tdC_mm
      tdset1_dy  = tdC_dy
      tdset1_hr  = tdC_hr
      tdset1_mn  = tdC_mn
      tdset1_hgt = tdC_hgt
  #    tdset1_prs = tdC_prs
#     tdset1_azm = tdC_azm
#     tdset1_elv = tdC_elv
      tdset1_spd = tdC_spd
      tdset1_dir = tdC_dir
      
        #----------------------------------------
        # if 'runtype' = 'match', get indices of matches and apply QC if bool_drv_qc=True.
        # if 'runtype' = 'plot', SKIP bool if-block

  if runtype == "match":
    if bool_qc:
      print("DAWN QC NOT APPLIED: TBA")
    elif not bool_qc:
      # Do not apply LOON QC

    #  sindexesDC = np.asarray(np.where(tdset1_hr==tdset1_hr))           #get all indices
      indexes1   = [0]#sindexesDC.flatten()

      qc_list = "No QC applied"

      d_lat = tdset1_lat
      d_lon = tdset1_lon
      d_yr  = tdset1_yr
      d_mm  = tdset1_mm
      d_dy  = tdset1_dy
      d_hr  = tdset1_hr
      d_mn  = tdset1_mn
   
      d_hgt = tdset1_hgt
#      d_azm = tdset1_azm
#      d_elv = tdset1_elv
      d_spd = tdset1_spd
      d_dir = tdset1_dir

  elif runtype == "plot":

    indexes1 = [0]

    qc_list = "No QC applied"

    d_lat = tdset1_lat[idxs]
    d_lon = tdset1_lon[idxs]
    d_yr  = tdset1_yr [idxs]
    d_mm  = tdset1_mm [idxs]
    d_dy  = tdset1_dy [idxs]
    d_hr  = tdset1_hr [idxs]
    d_mn  = tdset1_mn [idxs]
    d_hgt = tdset1_hgt[idxs]
#    d_azm = tdset1_azm[idxs]
#    d_elv = tdset1_elv[idxs]
    d_spd = tdset1_spd[idxs]
    d_dir = tdset1_dir[idxs]

        #----------------------------------------
        #----------------------------------------

        # create height array but fill with missing -999.
        #       AMV heights not available.
  d_prs = np.nan * np.ones_like(d_hr)
  
  
          # Return variables to MAIN
  return d_lat,d_lon,d_mm,d_yr,d_prs,d_dy,d_hr,d_mn,d_hgt,indexes1,qc_list,dset1_src,d_spd,d_dir

#===============================================================================================
# Read RADIOSONDE (from NCEP)
#
#	INPUTS:
#		path_prefix ......................... Full path up to directories /aeolus... or /atmos... : where this archive is located on local machine
#		yyyymmddhh .......................... Current date in yyyymmddhh format
#		bool_qc ............................. Choice to apply Aeolus QC: True=apply QC, False=don't apply QC
#
#	OUTPUTS:
#		d_lat ............................... Latitude in degrees [-90,90]
#		d_lon ............................... Longitude in degrees [0,360]
#		d_prs ............................... Pressure in hPa
#		d_hgt ............................... Height in km
#		d_yr ................................ Year
#		d_mm ................................ Month
#		d_dy ................................ Day
#		d_hr ................................ Hour
#		d_mn ................................ Minute
#		qc_list ............................. List of QC applied (if applicable)
#
def read_sonde(path_prefix,yyyymmddhh,dateB4,dateA,bool_qc,dsetflag,runtype,time_diff_max,idxs):

  qc_list = ""			#initialize

  yyyy = yyyymmddhh[0:4]
  mm   = yyyymmddhh[4:6]
  dd   = yyyymmddhh[6:8]
  hour = yyyymmddhh[8:10]

	#-------------------------------------
        # Find hour limits ... for Aeolus datafiles

  if hour == "00":
    hB4 = "18"
    hA  = "06"
  elif hour == "06":
    hB4 = "00"
    hA  = "12"
  elif hour == "12":
    hB4 = "06"
    hA  = "18"
  elif hour == "18":
    hB4 = "12"
    hA  = "00"

  if hour == "00":
    yyyymmddhhB4 = dateB4+hB4
    yyB4 = dateB4[0:4]
    mmB4 = dateB4[4:6]
    ddB4 = dateB4[6:8]
  else:
    yyyymmddhhB4 = yyyy+mm+dd+hB4
    yyB4 = yyyy
    mmB4 = mm
    ddB4 = dd
    
  if hour == "18":
    yyyymmddhhA = dateA+hA
    yyA  = dateA[0:4]
    mmA  = dateA[4:6]
    ddA  = dateA[6:8]
  else:
    yyyymmddhhA = yyyy+mm+dd+hA
    yyA = yyyy
    mmA = mm
    ddA = dd

		# convert 'time_max_diff' to rounded integer
  tqc_time = float(time_diff_max)
  if tqc_time%60 != 0.0:
    qc_time = int(math.ceil(tqc_time/60.0))             # round up to nearest or equal integer hour
  else:
    qc_time = int(tqc_time)

		# convert hour to integer
  ihr = int(hour)

		# find integer hours for before/after files for collocation
  if hour == "00":
    	# 00
    ihrB = 24 - qc_time
    ihrA = ihr + qc_time
  else:
	# 06,12,18
    ihrB = ihr - qc_time
    ihrA = ihr + qc_time

  del qc_time,tqc_time,ihr

	#-------------------------------------------------
    	# Define dataset

  tmp_dset1_path    = '/scratch/atmos-nc-dataset/sonde/'+yyyy+'/'+mm+'/'+dd+'/'
  tmp_dset1_path_B4 = '/scratch/atmos-nc-dataset/sonde/'+yyB4+'/'+mmB4+'/'+ddB4+'/'
  tmp_dset1_path_A  = '/scratch/atmos-nc-dataset/sonde/'+yyA+'/'+mmA+'/'+ddA+'/'

		# Path/file
		# ...Current date
  dset1_path   		= path_prefix+tmp_dset1_path
  dset1_filename 	= 'gdas.'+yyyymmddhh+'.adpupa.tm00.bufr_d.nc4'
  dset1_file   		= dset1_path+dset1_filename
  		# ...Before date
  dset1_path_B4		= path_prefix+tmp_dset1_path_B4
  dset1_filename_B4 	= 'gdas.'+yyyymmddhhB4+'.adpupa.tm00.bufr_d.nc4'
  dset1_file_B4		= dset1_path_B4+dset1_filename_B4
  		# ...After date
  dset1_path_A		= path_prefix+tmp_dset1_path_A
  dset1_filename_A 	= 'gdas.'+yyyymmddhhA+'.adpupa.tm00.bufr_d.nc4'
  dset1_file_A		= dset1_path_A+dset1_filename_A

    			# initialize flag indicating if dataset exists. 0=yes, 1=no
  existB = 0		# ... B = date Before current
  existA = 0		# ... A = date After current
  dset1_exists = exists(dset1_file)
  if dset1_exists==False:
    print("ERROR: file "+dset1_file+" does not exist!")
    sys.exit()							#exit script immediately
  dset1_existsB = exists(dset1_file_B4)
  if dset1_existsB==False:
    print("WARNING: 'before' file "+dset1_file_B4+" does not exist!")
    existB = 1
  dset1_existsA = exists(dset1_file_A)
  if dset1_existsA==False:
    print("WARNING: 'after' file "+dset1_file_A+" does not exist!")
    existA = 1

		# Path on FTP/web archive server (for output NetCDF only)
  str_dset1_path    = dset1_path
  str_dset1_path_B4 = dset1_path_B4
  str_dset1_path_A  = dset1_path_A

  dset1_src  = str_dset1_path_B4+dset1_filename_B4
  dset1_src += str_dset1_path+dset1_filename
  dset1_src += str_dset1_path_A+dset1_filename_A

		# Variable names
  dset1_lat_var = 'latitude'
  dset1_lon_var = 'longitude'
  dset1_yr_var  = 'year'
  dset1_mm_var  = 'month'
  dset1_dy_var  = 'day'
  dset1_hr_var  = 'hour'
  dset1_mn_var  = 'minutes'
  dset1_prs_var = 'pressure'
  dset1_ht_var  = 'height'
  dset1_spd_var = 'wind_speed'
  dset1_dir_var = 'wind_direction'
		
    	#-------------------------------------------------
  	# Load dataset
	#	RADIOSONDE data is divided into Groups
 
 		#```````````````````````````````````````````````
  		# CURRENT date
  nsondesC = []
  nlevelsC = []
  ngrpsC   = []
 
  data_hdl = Dataset(dset1_file)

  grps = list(data_hdl.groups)
  ngrpsC.append(np.size(grps))

	# populate full arrays with Group1 (grps(0))
  td_lat = np.asarray( data_hdl.groups[grps[0]].variables[dset1_lat_var] )
  td_lon = np.asarray( data_hdl.groups[grps[0]].variables[dset1_lon_var] )
  td_yr  = np.asarray( data_hdl.groups[grps[0]].variables[dset1_yr_var] )
  td_mm  = np.asarray( data_hdl.groups[grps[0]].variables[dset1_mm_var] )
  td_dy  = np.asarray( data_hdl.groups[grps[0]].variables[dset1_dy_var] )
  td_hr  = np.asarray( data_hdl.groups[grps[0]].variables[dset1_hr_var] )
  td_spd = np.asarray( data_hdl.groups[grps[0]].variables[dset1_spd_var] )
  td_dir = np.asarray( data_hdl.groups[grps[0]].variables[dset1_dir_var] )

  if grps[0]=='NC002001' or grps[0]=='NC002002' or grps[0]=='NC002003' or grps[0]=='NC002004' or grps[0]=='NC002005' or grps[0]=='NC002009':
    td_mn = np.asarray( data_hdl.groups[grps[0]].variables[dset1_hr_var]  )
    td_mn[:] = 0.0
  else:
    td_mn = np.asarray( data_hdl.groups[grps[0]].variables[dset1_mn_var] )
  
  if grps[0]=='NC002001' or grps[0]=='NC002002' or grps[0]=='NC002003' or grps[0]=='NC002004' or grps[0]=='NC002005' or grps[0]=='NC002006' or grps[0]=='NC002009' or grps[0]=='NC002015':
    	# groups with pressure variable but not height
    td_prs = np.asarray( data_hdl.groups[grps[0]].variables[dset1_prs_var]  )
    td_hgt = np.nan * np.ones_like(td_prs)
  else:
    	# groups with height variable but not pressure
    td_hgt = np.asarray( data_hdl.groups[grps[0]].variables[dset1_ht_var]  )
    td_prs = np.nan * np.ones_like(td_hgt)

	# find shape of array (nsondes, nlevels)
  grp_shape = np.shape(td_prs)
  nsondesC.append(grp_shape[0])

	# assign same lat,lon,yr,mm,dy,hr,mn to all levels per sonde (so that all vars have same dimensions)
  sttd_lat = np.nan * np.ones_like(td_prs)
  sttd_lon = np.nan * np.ones_like(td_prs)
  sttd_yr  = np.nan * np.ones_like(td_prs)
  sttd_mm  = np.nan * np.ones_like(td_prs)
  sttd_dy  = np.nan * np.ones_like(td_prs)
  sttd_hr  = np.nan * np.ones_like(td_prs)
  sttd_mn  = np.nan * np.ones_like(td_prs)
  sttd_prs = np.nan * np.ones_like(td_prs)
  sttd_hgt = np.nan * np.ones_like(td_prs)
  sttd_spd = np.nan * np.ones_like(td_prs)
  sttd_dir = np.nan * np.ones_like(td_prs)

  ttd_slv = np.nan * np.ones_like(td_lat)

  for isonde in range(grp_shape[0]):
    tmp = td_prs[isonde]
    ilev = np.where(tmp!=0.0)
    len_ilev = np.size(ilev)

    sttd_lat[isonde,ilev] = td_lat[isonde]
    sttd_lon[isonde,ilev] = td_lon[isonde]
    sttd_yr [isonde,ilev] = td_yr [isonde]
    sttd_mm [isonde,ilev] = td_mm [isonde]
    sttd_dy [isonde,ilev] = td_dy [isonde]
    sttd_hr [isonde,ilev] = td_hr [isonde]
    sttd_mn [isonde,ilev] = td_mn [isonde]
    sttd_prs[isonde,ilev] = td_prs[isonde,ilev]
    sttd_hgt[isonde,ilev] = td_hgt[isonde,ilev]
    sttd_spd[isonde,ilev] = td_spd[isonde,ilev]
    sttd_dir[isonde,ilev] = td_dir[isonde,ilev]
    
    ttd_slv[isonde] = len_ilev
    
    del tmp,ilev,len_ilev

	# find where arrays = nan and reassign as -999 (missing)
  ttd_lat = np.where(np.isnan(sttd_lat),-999.0,sttd_lat)
  ttd_lon = np.where(np.isnan(sttd_lon),-999.0,sttd_lon)
  ttd_yr  = np.where(np.isnan(sttd_yr ),-999.0,sttd_yr )
  ttd_mm  = np.where(np.isnan(sttd_mm ),-999.0,sttd_mm )
  ttd_dy  = np.where(np.isnan(sttd_dy ),-999.0,sttd_dy )
  ttd_hr  = np.where(np.isnan(sttd_hr ),-999.0,sttd_hr )
  ttd_mn  = np.where(np.isnan(sttd_mn ),-999.0,sttd_mn )
  ttd_prs = np.where(np.isnan(sttd_prs),-999.0,sttd_prs)
  ttd_hgt = np.where(np.isnan(sttd_hgt),-999.0,sttd_hgt)
  ttd_spd = np.where(np.isnan(sttd_spd),-999.0,sttd_spd)
  ttd_dir = np.where(np.isnan(sttd_dir),-999.0,sttd_dir)
    
  del td_lat,td_lon,td_yr,td_mm,td_dy,td_hr,td_mn,td_prs,td_hgt,td_spd,td_dir
  del sttd_lat,sttd_lon,sttd_yr,sttd_mm,sttd_dy,sttd_hr,sttd_mn,sttd_prs,sttd_hgt,sttd_spd,sttd_dir

	# conform arrays to 1D
	#	.flatten() COPIES the original array and conforms its dimensions to 1D
	#	Ex) a = [[1,2,3], [4,5,6]]
	#	    b = a.flatten()
	#	    print(b) --> prints b = [1,2,3,4,5,6]
  fdC_lat = ttd_lat.flatten()
  fdC_lon = ttd_lon.flatten()
  fdC_yr  = ttd_yr.flatten()
  fdC_mm  = ttd_mm.flatten()
  fdC_dy  = ttd_dy.flatten()
  fdC_hr  = ttd_hr.flatten()
  fdC_mn  = ttd_mn.flatten()
  fdC_prs = ttd_prs.flatten()
  fdC_hgt = ttd_hgt.flatten()
  fdC_spd = ttd_spd.flatten()
  fdC_dir = ttd_dir.flatten()
  fdC_slv = ttd_slv.flatten()
 
  del ttd_lat,ttd_lon,ttd_yr,ttd_mm,ttd_dy,ttd_hr,ttd_mn,ttd_prs,ttd_hgt,ttd_spd,ttd_dir,ttd_slv

	# append data from remaining Groups
  for x in range(len(grps)-1):
      tdset1_lat = np.asarray( data_hdl.groups[grps[x+1]].variables[dset1_lat_var] )
      tdset1_lon = np.asarray( data_hdl.groups[grps[x+1]].variables[dset1_lon_var] )
      tdset1_yr  = np.asarray( data_hdl.groups[grps[x+1]].variables[dset1_yr_var]  )
      tdset1_mm  = np.asarray( data_hdl.groups[grps[x+1]].variables[dset1_mm_var]  )
      tdset1_dy  = np.asarray( data_hdl.groups[grps[x+1]].variables[dset1_dy_var]  )
      tdset1_hr  = np.asarray( data_hdl.groups[grps[x+1]].variables[dset1_hr_var]  )
      tdset1_spd = np.asarray( data_hdl.groups[grps[x+1]].variables[dset1_spd_var]  )
      tdset1_dir = np.asarray( data_hdl.groups[grps[x+1]].variables[dset1_dir_var]  )

      if grps[x+1]=='NC002001' or grps[x+1]=='NC002002' or grps[x+1]=='NC002003' or grps[x+1]=='NC002004' or grps[x+1]=='NC002005' or grps[x+1]=='NC002009':
      	# groups without minutes variable
        tdset1_mn = np.asarray( data_hdl.groups[grps[x+1]].variables[dset1_hr_var]  )
        tdset1_mn[:] = 0.0
      else:
      	# groups with minutes variable
        tdset1_mn = np.asarray( data_hdl.groups[grps[x+1]].variables[dset1_mn_var]  )

      if grps[x+1]=='NC002001' or grps[x+1]=='NC002002' or grps[x+1]=='NC002003' or grps[x+1]=='NC002004' or grps[x+1]=='NC002005' or grps[x+1]=='NC002006' or grps[x+1]=='NC002009' or grps[x+1]=='NC002015':
      	# groups with pressure variable but not height
        p = np.asarray( data_hdl.groups[grps[x+1]].variables[dset1_prs_var]  )
        z = np.nan * np.ones_like(p)
        tgrp_shape = np.shape(p)
      else:
      	# groups with height variable but not pressure
        z = np.asarray( data_hdl.groups[grps[x+1]].variables[dset1_ht_var]  )
        p = np.nan * np.ones_like(z)
        tgrp_shape = np.shape(z)

	# find shape of array (nsondes, nlevels)
      tnsondes   = tgrp_shape[0]
      nsondesC.append(tnsondes)

	# assign same lat,lon,yr,mm,dy,hr,mn to all levels per sonde (so that all vars have same dimensions)
      sttdset1_lat = np.nan * np.ones_like(p)
      sttdset1_lon = np.nan * np.ones_like(p)
      sttdset1_yr  = np.nan * np.ones_like(p)
      sttdset1_mm  = np.nan * np.ones_like(p)
      sttdset1_dy  = np.nan * np.ones_like(p)
      sttdset1_hr  = np.nan * np.ones_like(p)
      sttdset1_mn  = np.nan * np.ones_like(p)
      sttdset1_prs = np.nan * np.ones_like(p)
      sttdset1_hgt = np.nan * np.ones_like(p)
      sttdset1_spd = np.nan * np.ones_like(p)
      sttdset1_dir = np.nan * np.ones_like(p)
      
      ttdset1_slv = np.nan * np.ones_like(tdset1_lat)

      for isonde in range(tgrp_shape[0]):
        tmp = p[isonde]
        ilev = np.where(tmp!=0.0)
        len_ilev = np.size(ilev)
      
        sttdset1_lat[isonde,ilev] = tdset1_lat[isonde]
        sttdset1_lon[isonde,ilev] = tdset1_lon[isonde]
        sttdset1_yr [isonde,ilev] = tdset1_yr [isonde]
        sttdset1_mm [isonde,ilev] = tdset1_mm [isonde]
        sttdset1_dy [isonde,ilev] = tdset1_dy [isonde]
        sttdset1_hr [isonde,ilev] = tdset1_hr [isonde]
        sttdset1_mn [isonde,ilev] = tdset1_mn [isonde]
        sttdset1_prs[isonde,ilev] = p[isonde,ilev]
        sttdset1_hgt[isonde,ilev] = z[isonde,ilev]
        sttdset1_spd[isonde,ilev] = tdset1_spd[isonde,ilev]
        sttdset1_dir[isonde,ilev] = tdset1_dir[isonde,ilev]
	
        ttdset1_slv[isonde] = len_ilev
    
        del tmp,ilev,len_ilev
 
	# find where arrays = nan and reassign as -999 (missing)
      ttdset1_lat = np.where(np.isnan(sttdset1_lat),-999.0,sttdset1_lat)
      ttdset1_lon = np.where(np.isnan(sttdset1_lon),-999.0,sttdset1_lon)
      ttdset1_yr  = np.where(np.isnan(sttdset1_yr ),-999.0,sttdset1_yr )
      ttdset1_mm  = np.where(np.isnan(sttdset1_mm ),-999.0,sttdset1_mm )
      ttdset1_dy  = np.where(np.isnan(sttdset1_dy ),-999.0,sttdset1_dy )
      ttdset1_hr  = np.where(np.isnan(sttdset1_hr ),-999.0,sttdset1_hr )
      ttdset1_mn  = np.where(np.isnan(sttdset1_mn ),-999.0,sttdset1_mn )
      ttdset1_prs = np.where(np.isnan(sttdset1_prs),-999.0,sttdset1_prs)
      ttdset1_hgt = np.where(np.isnan(sttdset1_hgt),-999.0,sttdset1_hgt)
      ttdset1_spd = np.where(np.isnan(sttdset1_spd),-999.0,sttdset1_spd)
      ttdset1_dir = np.where(np.isnan(sttdset1_dir),-999.0,sttdset1_dir) 
      
      del sttdset1_lat,sttdset1_lon,sttdset1_yr,sttdset1_mm,sttdset1_dy,sttdset1_hr,sttdset1_mn,sttdset1_prs,sttdset1_hgt,sttdset1_spd,sttdset1_dir
 
	# conform arrays to 1D
        #       .flatten() COPIES the original array and conforms its dimensions to 1D
        #       Ex) a = [[1,2,3], [4,5,6]]
        #           b = a.flatten()
        #           print(b) --> prints b = [1,2,3,4,5,6]
      fdset1_lat = ttdset1_lat.flatten()
      fdset1_lon = ttdset1_lon.flatten()
      fdset1_yr  = ttdset1_yr.flatten()
      fdset1_mm  = ttdset1_mm.flatten()
      fdset1_dy  = ttdset1_dy.flatten()
      fdset1_hr  = ttdset1_hr.flatten()
      fdset1_mn  = ttdset1_mn.flatten()
      fdset1_prs = ttdset1_prs.flatten()
      fdset1_hgt = ttdset1_hgt.flatten()
      fdset1_spd = ttdset1_spd.flatten()
      fdset1_dir = ttdset1_dir.flatten()
      fdset1_slv = ttdset1_slv.flatten()

	# append 1D arrays to first group's 1D arrays
      fdC_lat = np.append(fdC_lat,fdset1_lat,axis=0)
      fdC_lon = np.append(fdC_lon,fdset1_lon,axis=0)
      fdC_yr  = np.append(fdC_yr, fdset1_yr, axis=0)
      fdC_mm  = np.append(fdC_mm, fdset1_mm, axis=0)
      fdC_dy  = np.append(fdC_dy, fdset1_dy, axis=0)
      fdC_hr  = np.append(fdC_hr, fdset1_hr, axis=0)
      fdC_mn  = np.append(fdC_mn, fdset1_mn, axis=0)
      fdC_prs = np.append(fdC_prs,fdset1_prs,axis=0)
      fdC_hgt = np.append(fdC_hgt,fdset1_hgt,axis=0)
      fdC_spd = np.append(fdC_spd,fdset1_spd,axis=0)
      fdC_dir = np.append(fdC_dir,fdset1_dir,axis=0)
      fdC_slv = np.append(fdC_slv,fdset1_slv,axis=0)
      
      del tnsondes,tgrp_shape
      del tdset1_lat,tdset1_lon,tdset1_yr,tdset1_mm,tdset1_dy,tdset1_hr,tdset1_mn,tdset1_spd,tdset1_dir,p,z
      del ttdset1_lat,ttdset1_lon,ttdset1_yr,ttdset1_mm,ttdset1_dy,ttdset1_hr,ttdset1_mn,ttdset1_prs,ttdset1_hgt,ttdset1_spd,ttdset1_dir,ttdset1_slv
      del fdset1_lat,fdset1_lon,fdset1_yr,fdset1_mm,fdset1_dy,fdset1_hr,fdset1_mn,fdset1_prs,fdset1_hgt,fdset1_spd,fdset1_dir,fdset1_slv

  data_hdl.close()
  del grps,grp_shape
  
  if dsetflag == "dep":
    if existB == 0:
  		#```````````````````````````````````````````````
  		# BEFORE date
      nsondesB = []
      nlevelsB = []
      ngrpsB   = []
 
      data_hdl = Dataset(dset1_file_B4)

      grps = list(data_hdl.groups)
      ngrpsB.append(np.size(grps))

	# populate full arrays with Group1 (grps(0))
      td_lat = np.asarray( data_hdl.groups[grps[0]].variables[dset1_lat_var] )
      td_lon = np.asarray( data_hdl.groups[grps[0]].variables[dset1_lon_var] )
      td_yr  = np.asarray( data_hdl.groups[grps[0]].variables[dset1_yr_var] )
      td_mm  = np.asarray( data_hdl.groups[grps[0]].variables[dset1_mm_var] )
      td_dy  = np.asarray( data_hdl.groups[grps[0]].variables[dset1_dy_var] )
      td_hr  = np.asarray( data_hdl.groups[grps[0]].variables[dset1_hr_var] )
  
      td_spd = np.asarray( data_hdl.groups[grps[0]].variables[dset1_spd_var] )
      td_dir = np.asarray( data_hdl.groups[grps[0]].variables[dset1_dir_var] )

      if grps[0]=='NC002001' or grps[0]=='NC002002' or grps[0]=='NC002003' or grps[0]=='NC002004' or grps[0]=='NC002005' or grps[0]=='NC002009':
        td_mn = np.asarray( data_hdl.groups[grps[0]].variables[dset1_hr_var]  )
        td_mn[:] = 0.0
      else:
        td_mn = np.asarray( data_hdl.groups[grps[0]].variables[dset1_mn_var] )
  
      if grps[0]=='NC002001' or grps[0]=='NC002002' or grps[0]=='NC002003' or grps[0]=='NC002004' or grps[0]=='NC002005' or grps[0]=='NC002006' or grps[0]=='NC002009' or grps[0]=='NC002015':
    	# groups with pressure variable but not height
        td_prs = np.asarray( data_hdl.groups[grps[0]].variables[dset1_prs_var]  )
        td_hgt = np.nan * np.ones_like(td_prs)
      else:
    	# groups with height variable but not pressure
        td_hgt = np.asarray( data_hdl.groups[grps[0]].variables[dset1_ht_var]  )
        td_prs = np.nan * np.ones_like(td_hgt)

	# find shape of array (nsondes, nlevels)
      grp_shape = np.shape(td_prs)
      nsondesB.append(grp_shape[0])

	# assign same lat,lon,yr,mm,dy,hr,mn to all levels per sonde (so that all vars have same dimensions)
      sttd_lat = np.nan * np.ones_like(td_prs)
      sttd_lon = np.nan * np.ones_like(td_prs)
      sttd_yr  = np.nan * np.ones_like(td_prs)
      sttd_mm  = np.nan * np.ones_like(td_prs)
      sttd_dy  = np.nan * np.ones_like(td_prs)
      sttd_hr  = np.nan * np.ones_like(td_prs)
      sttd_mn  = np.nan * np.ones_like(td_prs)
      sttd_prs = np.nan * np.ones_like(td_prs)
      sttd_hgt = np.nan * np.ones_like(td_prs)
      sttd_spd = np.nan * np.ones_like(td_prs)
      sttd_dir = np.nan * np.ones_like(td_prs)
      
      ttd_slv = np.nan * np.ones_like(td_lat)
      
      for isonde in range(grp_shape[0]):
        tmp = td_prs[isonde]
        ilev = np.where(tmp!=0.0)
        len_ilev = np.size(ilev)

        sttd_lat[isonde,ilev] = td_lat[isonde]
        sttd_lon[isonde,ilev] = td_lon[isonde]
        sttd_yr [isonde,ilev] = td_yr [isonde]
        sttd_mm [isonde,ilev] = td_mm [isonde]
        sttd_dy [isonde,ilev] = td_dy [isonde]
        sttd_hr [isonde,ilev] = td_hr [isonde]
        sttd_mn [isonde,ilev] = td_mn [isonde]
        sttd_prs[isonde,ilev] = td_prs[isonde,ilev]
        sttd_hgt[isonde,ilev] = td_hgt[isonde,ilev]
        sttd_spd[isonde,ilev] = td_spd[isonde,ilev]
        sttd_dir[isonde,ilev] = td_dir[isonde,ilev]
        
        ttd_slv[isonde] = len_ilev
    
        del tmp,ilev,len_ilev	
      
      	  # find where arrays = nan and reassign as -999 (missing)
      ttd_lat = np.where(np.isnan(sttd_lat),-999.0,sttd_lat)
      ttd_lon = np.where(np.isnan(sttd_lon),-999.0,sttd_lon)
      ttd_yr  = np.where(np.isnan(sttd_yr ),-999.0,sttd_yr )
      ttd_mm  = np.where(np.isnan(sttd_mm ),-999.0,sttd_mm )
      ttd_dy  = np.where(np.isnan(sttd_dy ),-999.0,sttd_dy )
      ttd_hr  = np.where(np.isnan(sttd_hr ),-999.0,sttd_hr )
      ttd_mn  = np.where(np.isnan(sttd_mn ),-999.0,sttd_mn )
      ttd_prs = np.where(np.isnan(sttd_prs),-999.0,sttd_prs)
      ttd_hgt = np.where(np.isnan(sttd_hgt),-999.0,sttd_hgt)
      ttd_spd = np.where(np.isnan(sttd_spd),-999.0,sttd_spd)
      ttd_dir = np.where(np.isnan(sttd_dir),-999.0,sttd_dir)
    
      del td_lat,td_lon,td_yr,td_mm,td_dy,td_hr,td_mn,td_prs,td_hgt,td_spd,td_dir
      del sttd_lat,sttd_lon,sttd_yr,sttd_mm,sttd_dy,sttd_hr,sttd_mn,sttd_prs,sttd_hgt,sttd_spd,sttd_dir
  
          # conform arrays to 1D
          #	  .flatten() COPIES the original array and conforms its dimensions to 1D
          #	  Ex) a = [[1,2,3], [4,5,6]]
          #	      b = a.flatten()
          #	      print(b) --> prints b = [1,2,3,4,5,6]
      fdB_lat = ttd_lat.flatten()
      fdB_lon = ttd_lon.flatten()
      fdB_yr  = ttd_yr.flatten()
      fdB_mm  = ttd_mm.flatten()
      fdB_dy  = ttd_dy.flatten()
      fdB_hr  = ttd_hr.flatten()
      fdB_mn  = ttd_mn.flatten()
      fdB_prs = ttd_prs.flatten()
      fdB_hgt = ttd_hgt.flatten()
      fdB_spd = ttd_spd.flatten()
      fdB_dir = ttd_dir.flatten()
      fdB_slv = ttd_slv.flatten()

      del ttd_lat,ttd_lon,ttd_yr,ttd_mm,ttd_dy,ttd_hr,ttd_mn,ttd_prs,ttd_hgt,ttd_spd,ttd_dir

	# append data from remaining Groups
      for x in range(len(grps)-1):
        tdset1_lat = np.asarray( data_hdl.groups[grps[x+1]].variables[dset1_lat_var] )
        tdset1_lon = np.asarray( data_hdl.groups[grps[x+1]].variables[dset1_lon_var] )
        tdset1_yr  = np.asarray( data_hdl.groups[grps[x+1]].variables[dset1_yr_var]  )
        tdset1_mm  = np.asarray( data_hdl.groups[grps[x+1]].variables[dset1_mm_var]  )
        tdset1_dy  = np.asarray( data_hdl.groups[grps[x+1]].variables[dset1_dy_var]  )
        tdset1_hr  = np.asarray( data_hdl.groups[grps[x+1]].variables[dset1_hr_var]  )
        tdset1_spd = np.asarray( data_hdl.groups[grps[x+1]].variables[dset1_spd_var]  )
        tdset1_dir = np.asarray( data_hdl.groups[grps[x+1]].variables[dset1_dir_var]  )

        if grps[x+1]=='NC002001' or grps[x+1]=='NC002002' or grps[x+1]=='NC002003' or grps[x+1]=='NC002004' or grps[x+1]=='NC002005' or grps[x+1]=='NC002009':
  	  # groups without minutes variable
          tdset1_mn = np.asarray( data_hdl.groups[grps[x+1]].variables[dset1_hr_var]  )
          tdset1_mn[:] = 0.0
        else:
  	  # groups with minutes variable
          tdset1_mn = np.asarray( data_hdl.groups[grps[x+1]].variables[dset1_mn_var]  )
    
        if grps[x+1]=='NC002001' or grps[x+1]=='NC002002' or grps[x+1]=='NC002003' or grps[x+1]=='NC002004' or grps[x+1]=='NC002005' or grps[x+1]=='NC002006' or grps[x+1]=='NC002009' or grps[x+1]=='NC002015':
  	  # groups with pressure variable but not height
          p = np.asarray( data_hdl.groups[grps[x+1]].variables[dset1_prs_var]  )
          z = np.nan * np.ones_like(p)
          tgrp_shape = np.shape(p)
        else:
  	  # groups with height variable but not pressure
          z = np.asarray( data_hdl.groups[grps[x+1]].variables[dset1_ht_var]  )
          p = np.nan * np.ones_like(z)
          tgrp_shape = np.shape(z)
    
  	  # find shape of array (nsondes, nlevels)
        tnsondes   = tgrp_shape[0]
        nsondesB.append(tnsondes)

	# assign same lat,lon,yr,mm,dy,hr,mn to all levels per sonde (so that all vars have same dimensions)
        sttdset1_lat = np.nan * np.ones_like(p)
        sttdset1_lon = np.nan * np.ones_like(p)
        sttdset1_yr  = np.nan * np.ones_like(p)
        sttdset1_mm  = np.nan * np.ones_like(p)
        sttdset1_dy  = np.nan * np.ones_like(p)
        sttdset1_hr  = np.nan * np.ones_like(p)
        sttdset1_mn  = np.nan * np.ones_like(p)
        sttdset1_prs = np.nan * np.ones_like(p)
        sttdset1_hgt = np.nan * np.ones_like(p)
        sttdset1_spd = np.nan * np.ones_like(p)
        sttdset1_dir = np.nan * np.ones_like(p)
        
        ttdset1_slv = np.nan * np.ones_like(tdset1_lat)
  	
        for isonde in range(tgrp_shape[0]):
          tmp = p[isonde]
          ilev = np.where(tmp!=0.0)
          len_ilev = np.size(ilev)
	
          sttdset1_lat[isonde,ilev] = tdset1_lat[isonde]
          sttdset1_lon[isonde,ilev] = tdset1_lon[isonde]
          sttdset1_yr [isonde,ilev] = tdset1_yr [isonde]
          sttdset1_mm [isonde,ilev] = tdset1_mm [isonde]
          sttdset1_dy [isonde,ilev] = tdset1_dy [isonde]
          sttdset1_hr [isonde,ilev] = tdset1_hr [isonde]
          sttdset1_mn [isonde,ilev] = tdset1_mn [isonde]
          sttdset1_prs[isonde,ilev] = p[isonde,ilev]
          sttdset1_hgt[isonde,ilev] = z[isonde,ilev]
          sttdset1_spd[isonde,ilev] = tdset1_spd[isonde,ilev]
          sttdset1_dir[isonde,ilev] = tdset1_dir[isonde,ilev]
          
          ttdset1_slv[isonde] = len_ilev
    
          del tmp,ilev,len_ilev	
  	
		# find where arrays = nan and reassign as -999 (missing)
        ttdset1_lat = np.where(np.isnan(sttdset1_lat),-999.0,sttdset1_lat)
        ttdset1_lon = np.where(np.isnan(sttdset1_lon),-999.0,sttdset1_lon)
        ttdset1_yr  = np.where(np.isnan(sttdset1_yr ),-999.0,sttdset1_yr )
        ttdset1_mm  = np.where(np.isnan(sttdset1_mm ),-999.0,sttdset1_mm )
        ttdset1_dy  = np.where(np.isnan(sttdset1_dy ),-999.0,sttdset1_dy )
        ttdset1_hr  = np.where(np.isnan(sttdset1_hr ),-999.0,sttdset1_hr )
        ttdset1_mn  = np.where(np.isnan(sttdset1_mn ),-999.0,sttdset1_mn )
        ttdset1_prs = np.where(np.isnan(sttdset1_prs),-999.0,sttdset1_prs)
        ttdset1_hgt = np.where(np.isnan(sttdset1_hgt),-999.0,sttdset1_hgt)
        ttdset1_spd = np.where(np.isnan(sttdset1_spd),-999.0,sttdset1_spd)
        ttdset1_dir = np.where(np.isnan(sttdset1_dir),-999.0,sttdset1_dir) 
      
        del sttdset1_lat,sttdset1_lon,sttdset1_yr,sttdset1_mm,sttdset1_dy,sttdset1_hr,sttdset1_mn,sttdset1_prs,sttdset1_hgt,sttdset1_spd,sttdset1_dir	
	
  	  # conform arrays to 1D
    	  #	  .flatten() COPIES the original array and conforms its dimensions to 1D
  	  #	  Ex) a = [[1,2,3], [4,5,6]]
      	  #	      b = a.flatten()
  	  #	      print(b) --> prints b = [1,2,3,4,5,6]
        fdset1_lat = ttdset1_lat.flatten()
        fdset1_lon = ttdset1_lon.flatten()
        fdset1_yr  = ttdset1_yr.flatten()
        fdset1_mm  = ttdset1_mm.flatten()
        fdset1_dy  = ttdset1_dy.flatten()
        fdset1_hr  = ttdset1_hr.flatten()
        fdset1_mn  = ttdset1_mn.flatten()
        fdset1_prs = ttdset1_prs.flatten()
        fdset1_hgt = ttdset1_hgt.flatten()
        fdset1_spd = ttdset1_spd.flatten()
        fdset1_dir = ttdset1_dir.flatten()
        fdset1_slv = ttdset1_slv.flatten()

	# append 1D arrays to first group's 1D arrays
        fdB_lat = np.append(fdB_lat,fdset1_lat,axis=0)
        fdB_lon = np.append(fdB_lon,fdset1_lon,axis=0)
        fdB_yr  = np.append(fdB_yr, fdset1_yr, axis=0)
        fdB_mm  = np.append(fdB_mm, fdset1_mm, axis=0)
        fdB_dy  = np.append(fdB_dy, fdset1_dy, axis=0)
        fdB_hr  = np.append(fdB_hr, fdset1_hr, axis=0)
        fdB_mn  = np.append(fdB_mn, fdset1_mn, axis=0)
        fdB_prs = np.append(fdB_prs,fdset1_prs,axis=0)
        fdB_hgt = np.append(fdB_hgt,fdset1_hgt,axis=0)
        fdB_spd = np.append(fdB_spd,fdset1_spd,axis=0)
        fdB_dir = np.append(fdB_dir,fdset1_dir,axis=0)
        fdB_slv = np.append(fdB_slv,fdset1_slv,axis=0)
      
        del tnsondes,tgrp_shape
        del tdset1_lat,tdset1_lon,tdset1_yr,tdset1_mm,tdset1_dy,tdset1_hr,tdset1_mn,tdset1_spd,tdset1_dir,p,z
        del ttdset1_lat,ttdset1_lon,ttdset1_yr,ttdset1_mm,ttdset1_dy,ttdset1_hr,ttdset1_mn,ttdset1_prs,ttdset1_hgt,ttdset1_spd,ttdset1_dir,ttdset1_slv
        del fdset1_lat,fdset1_lon,fdset1_yr,fdset1_mm,fdset1_dy,fdset1_hr,fdset1_mn,fdset1_prs,fdset1_hgt,fdset1_spd,fdset1_dir,fdset1_slv

      data_hdl.close()
      del grps,grp_shape
  
    if existA == 0:
  		#```````````````````````````````````````````````
  		# AFTER date
      nsondesA = []
      nlevelsA = []
      ngrpsA   = []
 
      data_hdl = Dataset(dset1_file_A)

      grps = list(data_hdl.groups)
      ngrpsA.append(np.size(grps))

          # populate full arrays with Group1 (grps(0))
      td_lat = np.asarray( data_hdl.groups[grps[0]].variables[dset1_lat_var] )
      td_lon = np.asarray( data_hdl.groups[grps[0]].variables[dset1_lon_var] )
      td_yr  = np.asarray( data_hdl.groups[grps[0]].variables[dset1_yr_var] )
      td_mm  = np.asarray( data_hdl.groups[grps[0]].variables[dset1_mm_var] )
      td_dy  = np.asarray( data_hdl.groups[grps[0]].variables[dset1_dy_var] )
      td_hr  = np.asarray( data_hdl.groups[grps[0]].variables[dset1_hr_var] )  
      td_spd = np.asarray( data_hdl.groups[grps[0]].variables[dset1_spd_var] )
      td_dir = np.asarray( data_hdl.groups[grps[0]].variables[dset1_dir_var] )

      if grps[0]=='NC002001' or grps[0]=='NC002002' or grps[0]=='NC002003' or grps[0]=='NC002004' or grps[0]=='NC002005' or grps[0]=='NC002009':
        td_mn = np.asarray( data_hdl.groups[grps[0]].variables[dset1_hr_var]  )
        td_mn[:] = 0.0
      else:
        td_mn = np.asarray( data_hdl.groups[grps[0]].variables[dset1_mn_var] )
  
      if grps[0]=='NC002001' or grps[0]=='NC002002' or grps[0]=='NC002003' or grps[0]=='NC002004' or grps[0]=='NC002005' or grps[0]=='NC002006' or grps[0]=='NC002009' or grps[0]=='NC002015':
          # groups with pressure variable but not height
        td_prs = np.asarray( data_hdl.groups[grps[0]].variables[dset1_prs_var]  )
        td_hgt = np.nan * np.ones_like(td_prs)
      else:
          # groups with height variable but not pressure
        td_hgt = np.asarray( data_hdl.groups[grps[0]].variables[dset1_ht_var]  )
        td_prs = np.nan * np.ones_like(td_hgt)

          # find shape of array (nsondes, nlevels)
      grp_shape = np.shape(td_prs)
      nsondesA.append(grp_shape[0])

	# assign same lat,lon,yr,mm,dy,hr,mn to all levels per sonde (so that all vars have same dimensions)
      sttd_lat = np.nan * np.ones_like(td_prs)
      sttd_lon = np.nan * np.ones_like(td_prs)
      sttd_yr  = np.nan * np.ones_like(td_prs)
      sttd_mm  = np.nan * np.ones_like(td_prs)
      sttd_dy  = np.nan * np.ones_like(td_prs)
      sttd_hr  = np.nan * np.ones_like(td_prs)
      sttd_mn  = np.nan * np.ones_like(td_prs)
      sttd_prs = np.nan * np.ones_like(td_prs)
      sttd_hgt = np.nan * np.ones_like(td_prs)
      sttd_spd = np.nan * np.ones_like(td_prs)
      sttd_dir = np.nan * np.ones_like(td_prs)
      
      ttd_slv = np.nan * np.ones_like(td_lat)
      
      for isonde in range(grp_shape[0]):
        tmp = td_prs[isonde]
        ilev = np.where(tmp!=0.0)
        len_ilev = np.size(ilev)
      
        sttd_lat[isonde,ilev] = td_lat[isonde]
        sttd_lon[isonde,ilev] = td_lon[isonde]
        sttd_yr [isonde,ilev] = td_yr [isonde]
        sttd_mm [isonde,ilev] = td_mm [isonde]
        sttd_dy [isonde,ilev] = td_dy [isonde]
        sttd_hr [isonde,ilev] = td_hr [isonde]
        sttd_mn [isonde,ilev] = td_mn [isonde]
        sttd_prs[isonde,ilev] = td_prs[isonde,ilev]
        sttd_hgt[isonde,ilev] = td_hgt[isonde,ilev]
        sttd_spd[isonde,ilev] = td_spd[isonde,ilev]
        sttd_dir[isonde,ilev] = td_dir[isonde,ilev]
        
        ttd_slv[isonde] = len_ilev
    
        del tmp,ilev,len_ilev	
      
          # find where arrays = nan and reassign as -999 (missing)
      ttd_lat = np.where(np.isnan(sttd_lat),-999.0,sttd_lat)
      ttd_lon = np.where(np.isnan(sttd_lon),-999.0,sttd_lon)
      ttd_yr  = np.where(np.isnan(sttd_yr ),-999.0,sttd_yr )
      ttd_mm  = np.where(np.isnan(sttd_mm ),-999.0,sttd_mm )
      ttd_dy  = np.where(np.isnan(sttd_dy ),-999.0,sttd_dy )
      ttd_hr  = np.where(np.isnan(sttd_hr ),-999.0,sttd_hr )
      ttd_mn  = np.where(np.isnan(sttd_mn ),-999.0,sttd_mn )
      ttd_prs = np.where(np.isnan(sttd_prs),-999.0,sttd_prs)
      ttd_hgt = np.where(np.isnan(sttd_hgt),-999.0,sttd_hgt)
      ttd_spd = np.where(np.isnan(sttd_spd),-999.0,sttd_spd)
      ttd_dir = np.where(np.isnan(sttd_dir),-999.0,sttd_dir)
    
      del td_lat,td_lon,td_yr,td_mm,td_dy,td_hr,td_mn,td_prs,td_hgt,td_spd,td_dir
      del sttd_lat,sttd_lon,sttd_yr,sttd_mm,sttd_dy,sttd_hr,sttd_mn,sttd_prs,sttd_hgt,sttd_spd,sttd_dir
  
          # conform arrays to 1D
          #	  .flatten() COPIES the original array and conforms its dimensions to 1D
          #	  Ex) a = [[1,2,3], [4,5,6]]
          #	      b = a.flatten()
          #	      print(b) --> prints b = [1,2,3,4,5,6]
      fdA_lat = ttd_lat.flatten()
      fdA_lon = ttd_lon.flatten()
      fdA_yr  = ttd_yr.flatten()
      fdA_mm  = ttd_mm.flatten()
      fdA_dy  = ttd_dy.flatten()
      fdA_hr  = ttd_hr.flatten()
      fdA_mn  = ttd_mn.flatten()
      fdA_prs = ttd_prs.flatten()
      fdA_hgt = ttd_hgt.flatten()
      fdA_spd = ttd_spd.flatten()
      fdA_dir = ttd_dir.flatten()
      fdA_slv = ttd_slv.flatten()
  
      del ttd_lat,ttd_lon,ttd_yr,ttd_mm,ttd_dy,ttd_hr,ttd_mn,ttd_prs,ttd_hgt,ttd_spd,ttd_dir,ttd_slv

	# append data from remaining Groups
      for x in range(len(grps)-1):
        tdset1_lat = np.asarray( data_hdl.groups[grps[x+1]].variables[dset1_lat_var] )
        tdset1_lon = np.asarray( data_hdl.groups[grps[x+1]].variables[dset1_lon_var] )
        tdset1_yr  = np.asarray( data_hdl.groups[grps[x+1]].variables[dset1_yr_var]  )
        tdset1_mm  = np.asarray( data_hdl.groups[grps[x+1]].variables[dset1_mm_var]  )
        tdset1_dy  = np.asarray( data_hdl.groups[grps[x+1]].variables[dset1_dy_var]  )
        tdset1_hr  = np.asarray( data_hdl.groups[grps[x+1]].variables[dset1_hr_var]  )
        tdset1_spd = np.asarray( data_hdl.groups[grps[x+1]].variables[dset1_spd_var]  )
        tdset1_dir = np.asarray( data_hdl.groups[grps[x+1]].variables[dset1_dir_var]  )

        if grps[x+1]=='NC002001' or grps[x+1]=='NC002002' or grps[x+1]=='NC002003' or grps[x+1]=='NC002004' or grps[x+1]=='NC002005' or grps[x+1]=='NC002009':
  	  # groups without minutes variable
          tdset1_mn = np.asarray( data_hdl.groups[grps[x+1]].variables[dset1_hr_var]  )
          tdset1_mn[:] = 0.0
        else:
  	  # groups with minutes variable
          tdset1_mn = np.asarray( data_hdl.groups[grps[x+1]].variables[dset1_mn_var]  )
    
        if grps[x+1]=='NC002001' or grps[x+1]=='NC002002' or grps[x+1]=='NC002003' or grps[x+1]=='NC002004' or grps[x+1]=='NC002005' or grps[x+1]=='NC002006' or grps[x+1]=='NC002009' or grps[x+1]=='NC002015':
  	  # groups with pressure variable but not height
          p = np.asarray( data_hdl.groups[grps[x+1]].variables[dset1_prs_var]  )
          z = np.nan * np.ones_like(p)
          tgrp_shape = np.shape(p)
        else:
  	  # groups with height variable but not pressure
          z = np.asarray( data_hdl.groups[grps[x+1]].variables[dset1_ht_var]  )
          p = np.nan * np.ones_like(z)
          tgrp_shape = np.shape(z)
    
  	  # find shape of array (nsondes, nlevels)
        tnsondes   = tgrp_shape[0]
        nsondesA.append(tnsondes)

  	  # assign same lat,lon,yr,mm,dy,hr,mn to all levels per sonde (so that all vars have same dimensions)
        sttdset1_lat = np.nan * np.ones_like(p)
        sttdset1_lon = np.nan * np.ones_like(p)
        sttdset1_yr  = np.nan * np.ones_like(p)
        sttdset1_mm  = np.nan * np.ones_like(p)
        sttdset1_dy  = np.nan * np.ones_like(p)
        sttdset1_hr  = np.nan * np.ones_like(p)
        sttdset1_mn  = np.nan * np.ones_like(p)
        sttdset1_prs = np.nan * np.ones_like(p)
        sttdset1_hgt = np.nan * np.ones_like(p)
        sttdset1_spd = np.nan * np.ones_like(p)
        sttdset1_dir = np.nan * np.ones_like(p)
        
        ttdset1_slv = np.nan * np.ones_like(tdset1_lat)
   
        for isonde in range(tgrp_shape[0]):
          tmp = p[isonde]
          ilev = np.where(tmp!=0.0)
          len_ilev = np.size(ilev)
	
          sttdset1_lat[isonde,ilev] = tdset1_lat[isonde]
          sttdset1_lon[isonde,ilev] = tdset1_lon[isonde]
          sttdset1_yr [isonde,ilev] = tdset1_yr [isonde]
          sttdset1_mm [isonde,ilev] = tdset1_mm [isonde]
          sttdset1_dy [isonde,ilev] = tdset1_dy [isonde]
          sttdset1_hr [isonde,ilev] = tdset1_hr [isonde]
          sttdset1_mn [isonde,ilev] = tdset1_mn [isonde]
          sttdset1_prs[isonde,ilev] = p[isonde,ilev]
          sttdset1_hgt[isonde,ilev] = z[isonde,ilev]
          sttdset1_spd[isonde,ilev] = tdset1_spd[isonde,ilev]
          sttdset1_dir[isonde,ilev] = tdset1_dir[isonde,ilev]
          
          ttdset1_slv[isonde] = len_ilev
    
          del tmp,ilev,len_ilev
  	
		# find where arrays = nan and reassign as -999 (missing)
        ttdset1_lat = np.where(np.isnan(sttdset1_lat),-999.0,sttdset1_lat)
        ttdset1_lon = np.where(np.isnan(sttdset1_lon),-999.0,sttdset1_lon)
        ttdset1_yr  = np.where(np.isnan(sttdset1_yr ),-999.0,sttdset1_yr )
        ttdset1_mm  = np.where(np.isnan(sttdset1_mm ),-999.0,sttdset1_mm )
        ttdset1_dy  = np.where(np.isnan(sttdset1_dy ),-999.0,sttdset1_dy )
        ttdset1_hr  = np.where(np.isnan(sttdset1_hr ),-999.0,sttdset1_hr )
        ttdset1_mn  = np.where(np.isnan(sttdset1_mn ),-999.0,sttdset1_mn )
        ttdset1_prs = np.where(np.isnan(sttdset1_prs),-999.0,sttdset1_prs)
        ttdset1_hgt = np.where(np.isnan(sttdset1_hgt),-999.0,sttdset1_hgt)
        ttdset1_spd = np.where(np.isnan(sttdset1_spd),-999.0,sttdset1_spd)
        ttdset1_dir = np.where(np.isnan(sttdset1_dir),-999.0,sttdset1_dir) 
      
        del sttdset1_lat,sttdset1_lon,sttdset1_yr,sttdset1_mm,sttdset1_dy,sttdset1_hr,sttdset1_mn,sttdset1_prs,sttdset1_hgt,sttdset1_spd,sttdset1_dir	
	
    	  # conform arrays to 1D
  	  #	  .flatten() COPIES the original array and conforms its dimensions to 1D
  	  #	  Ex) a = [[1,2,3], [4,5,6]]
  	  #	      b = a.flatten()
  	  #	      print(b) --> prints b = [1,2,3,4,5,6]
        fdset1_lat = ttdset1_lat.flatten()
        fdset1_lon = ttdset1_lon.flatten()
        fdset1_yr  = ttdset1_yr.flatten()
        fdset1_mm  = ttdset1_mm.flatten()
        fdset1_dy  = ttdset1_dy.flatten()
        fdset1_hr  = ttdset1_hr.flatten()
        fdset1_mn  = ttdset1_mn.flatten()
        fdset1_prs = ttdset1_prs.flatten()
        fdset1_hgt = ttdset1_hgt.flatten()
        fdset1_spd = ttdset1_spd.flatten()
        fdset1_dir = ttdset1_dir.flatten()
        fdset1_slv = ttdset1_slv.flatten()
      
  	  # append 1D arrays to first group's 1D arrays
        fdA_lat = np.append(fdA_lat,fdset1_lat,axis=0)
        fdA_lon = np.append(fdA_lon,fdset1_lon,axis=0)
        fdA_yr  = np.append(fdA_yr, fdset1_yr, axis=0)
        fdA_mm  = np.append(fdA_mm, fdset1_mm, axis=0)
        fdA_dy  = np.append(fdA_dy, fdset1_dy, axis=0)
        fdA_hr  = np.append(fdA_hr, fdset1_hr, axis=0)
        fdA_mn  = np.append(fdA_mn, fdset1_mn, axis=0)
        fdA_prs = np.append(fdA_prs,fdset1_prs,axis=0)
        fdA_hgt = np.append(fdA_hgt,fdset1_hgt,axis=0)
        fdA_spd = np.append(fdA_spd,fdset1_spd,axis=0)
        fdA_dir = np.append(fdA_dir,fdset1_dir,axis=0)
        fdA_slv = np.append(fdA_slv,fdset1_slv,axis=0)
  
        del tnsondes,tgrp_shape
        del tdset1_lat,tdset1_lon,tdset1_yr,tdset1_mm,tdset1_dy,tdset1_hr,tdset1_mn,tdset1_spd,tdset1_dir,p,z
        del ttdset1_lat,ttdset1_lon,ttdset1_yr,ttdset1_mm,ttdset1_dy,ttdset1_hr,ttdset1_mn,ttdset1_prs,ttdset1_hgt,ttdset1_spd,ttdset1_dir,ttdset1_slv
        del fdset1_lat,fdset1_lon,fdset1_yr,fdset1_mm,fdset1_dy,fdset1_hr,fdset1_mn,fdset1_prs,fdset1_hgt,fdset1_spd,fdset1_dir,fdset1_slv

      data_hdl.close()
      del grps,grp_shape
  
  	#-------------------------------------------------
      	# Append current arrays to before date arrays

		# get BEFORE and AFTER obs
    if hour=="00":
    	# keep ihrB hrs before CURRENT and ihrA hrs after CURRENT 6hrs
      tiHHA  = np.where(((fdC_hr >= 3) * (fdC_hr < (3+ihrA))))
      siHHA  = np.asarray(tiHHA)
      iHHA   = siHHA.flatten()

      sfdA_lat = fdC_lat[iHHA]
      sfdA_lon = fdC_lon[iHHA]
      sfdA_yr  = fdC_yr [iHHA]
      sfdA_mm  = fdC_mm [iHHA]
      sfdA_dy  = fdC_dy [iHHA]
      sfdA_hr  = fdC_hr [iHHA]
      sfdA_mn  = fdC_mn [iHHA]
      sfdA_spd = fdC_spd[iHHA]
      sfdA_dir = fdC_dir[iHHA]
      sfdA_prs = fdC_prs[iHHA]
      sfdA_hgt = fdC_hgt[iHHA]
      sfdA_slv = fdC_slv
	
      if existB == 0:
        tiHHB4 = np.where(((fdB_hr >= (21-ihrB)) * (fdB_hr < 21)))
        siHHB4 = np.asarray(tiHHB4)
        iHHB4  = siHHB4.flatten()
	
        sfdB_lat = fdB_lat[iHHB4]
        sfdB_lon = fdB_lon[iHHB4]
        sfdB_yr  = fdB_yr [iHHB4]
        sfdB_mm  = fdB_mm [iHHB4]
        sfdB_dy  = fdB_dy [iHHB4]
        sfdB_hr  = fdB_hr [iHHB4]
        sfdB_mn  = fdB_mn [iHHB4]
        sfdB_spd = fdB_spd[iHHB4]
        sfdB_dir = fdB_dir[iHHB4]
        sfdB_prs = fdB_prs[iHHB4]
        sfdB_hgt = fdB_hgt[iHHB4]
        sfdB_slv = fdB_slv
	
    elif hour=="06":
    	# keep ihrB before CURRENT and ihrA after CURRENT 6hrs
      tiHHA   = np.where(((fdC_hr >= 9) * (fdC_hr < (9+ihrA))))
      siHHA   = np.asarray(tiHHA)
      iHHA    = siHHA.flatten()
      
      sfdA_lat = fdC_lat[iHHA]
      sfdA_lon = fdC_lon[iHHA]
      sfdA_yr  = fdC_yr [iHHA]
      sfdA_mm  = fdC_mm [iHHA]
      sfdA_dy  = fdC_dy [iHHA]
      sfdA_hr  = fdC_hr [iHHA]
      sfdA_mn  = fdC_mn [iHHA]
      sfdA_spd = fdC_spd[iHHA]
      sfdA_dir = fdC_dir[iHHA]
      sfdA_prs = fdC_prs[iHHA]
      sfdA_hgt = fdC_hgt[iHHA]
      sfdA_slv = fdC_slv
      
      if (3-ihrB) < 0:
        dhrB = abs(3 - ihrB)
      
        tiHHB4 = np.where(((fdC_hr >= 0) * (fdC_hr < 3)))
        siHHB4 = np.asarray(tiHHB4)
        iHHB4  = siHHB4.flatten()
      
        if existB == 0:
          tiHHB41 = np.where(((fdB_hr >= (24-dhrB)) * (fdB_hr < 24)))  
          siHHB41 = np.asarray(tiHHB41)
          iHHB41  = siHHB41.flatten()
	  
          sfdB_lat = np.append(fdB_lat[iHHB41], fdC_lat[iHHB4], axis=0)
          sfdB_lon = np.append(fdB_lon[iHHB41], fdC_lon[iHHB4], axis=0)
          sfdB_yr  = np.append(fdB_yr [iHHB41], fdC_yr [iHHB4], axis=0)
          sfdB_mm  = np.append(fdB_mm [iHHB41], fdC_mm [iHHB4], axis=0)
          sfdB_dy  = np.append(fdB_dy [iHHB41], fdC_dy [iHHB4], axis=0)
          sfdB_hr  = np.append(fdB_hr [iHHB41], fdC_hr [iHHB4], axis=0)
          sfdB_mn  = np.append(fdB_mn [iHHB41], fdC_mn [iHHB4], axis=0)
          sfdB_spd = np.append(fdB_spd[iHHB41], fdC_spd[iHHB4], axis=0)
          sfdB_dir = np.append(fdB_dir[iHHB41], fdC_dir[iHHB4], axis=0)
          sfdB_prs = np.append(fdB_prs[iHHB41], fdC_prs[iHHB4], axis=0)
          sfdB_hgt = np.append(fdB_hgt[iHHB41], fdC_hgt[iHHB4], axis=0)
          sfdB_slv = np.append(fdB_slv        , fdC_slv       , axis=0)
        else:
          sfdB_lat = fdC_lat[iHHB4]
          sfdB_lon = fdC_lon[iHHB4]
          sfdB_yr  = fdC_yr [iHHB4]
          sfdB_mm  = fdC_mm [iHHB4]
          sfdB_dy  = fdC_dy [iHHB4]
          sfdB_hr  = fdC_hr [iHHB4]
          sfdB_mn  = fdC_mn [iHHB4]
          sfdB_spd = fdC_spd[iHHB4]
          sfdB_dir = fdC_dir[iHHB4]
          sfdB_prs = fdC_prs[iHHB4]
          sfdB_hgt = fdC_hgt[iHHB4]
          sfdB_slv = fdC_slv

        del dhrB
	
      else:
        tiHHB4 = np.where(((fdC_hr >= (3-ihrB)) * (fdC_hr < 3)))
        siHHB4 = np.asarray(tiHHB4)
        iHHB4  = siHHB4.flatten()
	
        sfdB_lat = fdC_lat[iHHB4]
        sfdB_lon = fdC_lon[iHHB4]
        sfdB_yr  = fdC_yr [iHHB4]
        sfdB_mm  = fdC_mm [iHHB4]
        sfdB_dy  = fdC_dy [iHHB4]
        sfdB_hr  = fdC_hr [iHHB4]
        sfdB_mn  = fdC_mn [iHHB4]
        sfdB_spd = fdC_spd[iHHB4]
        sfdB_dir = fdC_dir[iHHB4]
        sfdB_prs = fdC_prs[iHHB4]
        sfdB_hgt = fdC_hgt[iHHB4]
        sfdB_slv = fdC_slv
	
    elif hour=="12":
    	# keep ihrB before CURRENT and ihrA after CURRENT 6hrs
      tiHHA  = np.where(((fdC_hr >= 15) * (fdC_hr < (15+ihrA)))) 
      tiHHB4 = np.where(((fdC_hr >= (9-ihrB)) * (fdC_hr < 9)))
    
      siHHA  = np.asarray(tiHHA)
      iHHA   = siHHA.flatten()
      siHHB4 = np.asarray(tiHHB4)
      iHHB4  = siHHB4.flatten()
      
      sfdA_lat = fdC_lat[iHHA]
      sfdA_lon = fdC_lon[iHHA]
      sfdA_yr  = fdC_yr [iHHA]
      sfdA_mm  = fdC_mm [iHHA]
      sfdA_dy  = fdC_dy [iHHA]
      sfdA_hr  = fdC_hr [iHHA]
      sfdA_mn  = fdC_mn [iHHA]
      sfdA_spd = fdC_spd[iHHA]
      sfdA_dir = fdC_dir[iHHA]
      sfdA_prs = fdC_prs[iHHA]
      sfdA_hgt = fdC_hgt[iHHA]
      sfdA_slv = fdC_slv
     
      sfdB_lat = fdC_lat[iHHB4]
      sfdB_lon = fdC_lon[iHHB4]
      sfdB_yr  = fdC_yr [iHHB4]
      sfdB_mm  = fdC_mm [iHHB4]
      sfdB_dy  = fdC_dy [iHHB4]
      sfdB_hr  = fdC_hr [iHHB4]
      sfdB_mn  = fdC_mn [iHHB4]
      sfdB_spd = fdC_spd[iHHB4]
      sfdB_dir = fdC_dir[iHHB4]
      sfdB_prs = fdC_prs[iHHB4]
      sfdB_hgt = fdC_hgt[iHHB4]
      sfdB_slv = fdC_slv
    
    elif hour=="18":
    	# keep ihrB before CURRENT and ihrA after CURRENT 6hrs
      if (21+ihrA) >= 24:
        dhrA = abs(24 - (21 + ihrA))

        tiHHA = np.where(((fdC_hr >= 21) * (fdC_hr < 24)))
        siHHA = np.asarray(tiHHA)
        iHHA  = siHHA.flatten()
      
        if existA == 0:
          tiHHA2 = np.where(((fdA_hr >= 0) * (fdA_hr < dhrA)))  
          siHHA2 = np.asarray(tiHHA2)
          iHHA2  = siHHA2.flatten()
	  
          sfdA_lat = np.append(fdC_lat[iHHA], fdA_lat[iHHA2], axis=0)
          sfdA_lon = np.append(fdC_lon[iHHA], fdA_lon[iHHA2], axis=0)
          sfdA_yr  = np.append(fdC_yr [iHHA], fdA_yr [iHHA2], axis=0)
          sfdA_mm  = np.append(fdC_mm [iHHA], fdA_mm [iHHA2], axis=0)
          sfdA_dy  = np.append(fdC_dy [iHHA], fdA_dy [iHHA2], axis=0)
          sfdA_hr  = np.append(fdC_hr [iHHA], fdA_hr [iHHA2], axis=0)
          sfdA_mn  = np.append(fdC_mn [iHHA], fdA_mn [iHHA2], axis=0)
          sfdA_spd = np.append(fdC_spd[iHHA], fdA_spd[iHHA2], axis=0)
          sfdA_dir = np.append(fdC_dir[iHHA], fdA_dir[iHHA2], axis=0)
          sfdA_prs = np.append(fdC_prs[iHHA], fdA_prs[iHHA2], axis=0)
          sfdA_hgt = np.append(fdC_hgt[iHHA], fdA_hgt[iHHA2], axis=0)
          sfdA_slv = np.append(fdC_slv      , fdA_slv       , axis=0)
        else:
          sfdA_lat = fdC_lat[iHHA]
          sfdA_lon = fdC_lon[iHHA]
          sfdA_yr  = fdC_yr [iHHA]
          sfdA_mm  = fdC_mm [iHHA]
          sfdA_dy  = fdC_dy [iHHA]
          sfdA_hr  = fdC_hr [iHHA]
      nsondes = np.append(nsondesB,nsondesC,axis=0)
      nlevels = np.append(sfdB_slv,fdC_slv,axis=0)
      ngrps   = np.append(ngrpsB,ngrpsC,axis=0)
    else:
      fd_lat = fdC_lat
      fd_lon = fdC_lon
      fd_yr  = fdC_yr 
      fd_mm  = fdC_mm 
      fd_dy  = fdC_dy 
      fd_hr  = fdC_hr 
      fd_mn  = fdC_mn 
      fd_prs = fdC_prs
      fd_hgt = fdC_hgt
      fd_spd = fdC_spd
      fd_dir = fdC_dir
      
      nsondes = nsondesC
      nlevels = fdC_slv
      ngrps   = ngrpsC

  		# append AFTER to CURRENT
    if existA == 0:
      fd_lat = np.append(fd_lat,sfdA_lat,axis=0)
      fd_lon = np.append(fd_lon,sfdA_lon,axis=0)
      fd_yr  = np.append(fd_yr ,sfdA_yr ,axis=0)
      fd_mm  = np.append(fd_mm ,sfdA_mm ,axis=0)
      fd_dy  = np.append(fd_dy ,sfdA_dy ,axis=0)
      fd_hr  = np.append(fd_hr ,sfdA_hr ,axis=0)
      fd_mn  = np.append(fd_mn ,sfdA_mn ,axis=0)
      fd_prs = np.append(fd_prs,sfdA_prs,axis=0)
      fd_hgt = np.append(fd_hgt,sfdA_hgt,axis=0)
      fd_spd = np.append(fd_spd,sfdA_spd,axis=0)
      fd_dir = np.append(fd_dir,sfdA_dir,axis=0)
  
      nsondes = np.append(nsondes,nsondesA,axis=0)
      nlevels = np.append(nlevels,sfdA_slv,axis=0)
      ngrps   = np.append(ngrps,ngrpsA,axis=0)
  
  else:
  	# dsetflag = "drv"
    fd_lat = fdC_lat
    fd_lon = fdC_lon
    fd_yr  = fdC_yr 
    fd_mm  = fdC_mm 
    fd_dy  = fdC_dy 
    fd_hr  = fdC_hr 
    fd_mn  = fdC_mn 
    fd_prs = fdC_prs
    fd_hgt = fdC_hgt
    fd_spd = fdC_spd
    fd_dir = fdC_dir
    
    nsondes = nsondesC
    nlevels = fdC_slv
    ngrps   = ngrpsC
 
	#----------------------------------------
	# if 'runtype' = 'match', get indices of matches and apply QC if bool_drv_qc=True.
	# if 'runtype' = 'plot', SKIP bool if-block
	
  if runtype == "match":
    if bool_qc:
      print("RADIOSONDE QC TBA")	
    elif not bool_qc:
      # Do not apply RADIOSONDE QC

      sindexesDC = np.asarray(np.where(fd_prs==fd_prs))           #get all indices
      indexes1   = sindexesDC.flatten()

      qc_list = "No QC applied"

      d_lat = fd_lat
      d_lon = fd_lon
      d_prs = fd_prs
      d_hgt = fd_hgt
      d_yr  = fd_yr
      d_mm  = fd_mm
      d_dy  = fd_dy
      d_hr  = fd_hr
      d_mn  = fd_mn
      d_spd = fd_spd
      d_dir = fd_dir

  elif runtype == "plot":

    indexes1 = [0]

    qc_list = "No QC applied"

    d_lat = fd_lat[idxs]
    d_lon = fd_lon[idxs]
    d_prs = fd_prs[idxs]
    d_hgt = fd_hgt[idxs]
    d_yr  = fd_yr [idxs]
    d_mm  = fd_mm [idxs]
    d_dy  = fd_dy [idxs]
    d_hr  = fd_hr [idxs]
    d_mn  = fd_mn [idxs]
    d_spd = fd_spd[idxs]
    d_dir = fd_dir[idxs]

	#----------------------------------------
	
	# check pressure units and convert to hPa
  if max(d_prs) > 10000.:
    d_prs = d_prs/100.

	# check height units and convert to km
  if max(d_hgt) > 1000.:
    d_hgt = d_hgt/1000.
 
# Return variables to MAIN
  return d_lat,d_lon,d_yr,d_mm,d_dy,d_hr,d_mn,d_hgt,d_prs,indexes1,qc_list,dset1_src,nsondes,nlevels,ngrps,d_spd,d_dir
  
#===============================================================================================
