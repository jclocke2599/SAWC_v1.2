#########################################################################################################
# MAIN SCRIPT FOR Collocation of Winds from Multiple Datasets
# 	-- This script has been specifically modified to collocate Aircraft winds (dataset 1) and AMVs (dataset 2) with Aeolus winds (driver dataset)
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
#	        David Santek
#
# Modified by:  Katherine E. Lukens
#
# HISTORY: 
#	2021		B. Hoover, D. Santek    Created/developed program.
#	2021-10-25	K.E. Lukens		Added aircraft, AMV input configs. Added quality controls (QC).
#	2021-11-22	K.E. Lukens		Added date before consideration.
#	2021-11-23	K.E. Lukens		Created QC module and added relevant calls.
#	2022-03-08	K.E. Lukens		Added metadata to output files.
#	2023-03-29	K.E. Lukens		Finalized for upload to wind-datasets archive.
#	2023-07-15	D. Huber, K.E. Lukens	Optimized program to improve performance.
#
# Usage: python3 ThisScriptName.py $driver_wind_type $TF_drv_qc $TF_dset2_qc $pct $dateIN $driver_dset_type_str
#
#	driver_wind_type ..................... Aeolus wind type: RayClear = Rayleigh-clear winds, MieCloud = Mie-cloudy winds
#	TF_drv_qc ............................ Choose to apply Aeolus QC as listed in ./quality_controls.py: 0 = no QC, 1 = apply QC
#	TF_dset2_qc .......................... Choose to apply AMV QC as listed in ./quality_controls.py: 0 = no QC, 1 = apply QC
#	pct .................................. Percent confidence (quality indicator in %) for AMV data. If TF_dset2_qc=0, this value will be ignored
#	dateIN ............................... Date over which to perform collocation. Format is yyyymmddhh (yyyy = year, mm = month, dd = day, hh = hour (can only be 00,06,12, or 18))
#	driver_dset_type_str.................. Choose if Aeolus reprocessed data will be used: orig = not reprocessed (original data), for reprocessed files use baseline number (example: B10)
#
# Input: NetCDF
#	Driver dataset ....................... Daily Aeolus winds archived on S4, in yyyymmdd format
#	  -- Driver is independent dataset: All other datasets are collocated to Driver obs.
#	Dependent datasets ............................ 6-hourly winds archived on S4 (converted to NetCDF), in yyyymmddhh format
#	  -- Datasets are collocated to Driver obs.
#
# Output: NetCDF
#	Files contain:
#		1. Indices of collocated observations. Indices pertain to arrays in each source data file.
#		2. Metadata listing collocation criteria, QC criteria (if applicable), paths to source data
#
#########################################################################################################
# Import python modules
#########################################################################################################
print("***** Import Python Modules *****")

import os
os.environ['OPENBLAS_NUM_THREADS'] = '1'

from os.path import exists

import sys
import numpy as np
import matplotlib.pyplot as plt

from netCDF4 import Dataset

from match_dependencies import fix_lon
from match_dependencies import fix_yr
from match_dependencies import fix_mn
from match_dependencies import fix_dy
from match_dependencies import fix_hr
from match_dependencies import fix_mn
from match_dependencies import match_obs
from match_dependencies import obs_match_3d
from match_dependencies import match_lists_to_netCDF

from read_data import read_aeolus
from read_data import read_aircraft
from read_data import read_amv_ncep
from read_data import read_loon
from read_data import read_sonde

##############################################
print("========================================")
print("========== BEGIN MAIN PROGRAM ==========")

runtype = "match"	# type of program coded here: MATCH code

#=============================================
# Read raw user input from command line

aeolus_wind_type        = sys.argv[1]           # Aeolus wind type
qc_flags                = sys.argv[2]           # Apply QC to datasets? 0=no, 1=yes
pct_str                 = sys.argv[3]      	# AMV Quality Indicator (QI) lower limits in % (per AMV dataset listed)
qi_choice_str           = sys.argv[4]      	# AMV Quality Indicator (QI) choice: NO_FC (default) is QI without forecast; YES_FC is QI with forecast
dateIN                  = sys.argv[5]           # Date in question: YYYYMMDDHH
aeolus_dset_type_str    = sys.argv[6]           # Aeolus dataset type: orig = original, for reprocessed files use baseline number (example: B10)
output_path             = sys.argv[7]           # Output (archive) directory
dst_max_str 		= sys.argv[8]	        # Collocation distance maximum in km ... horizontal distance
prs_max_str 		= sys.argv[9]    	# Collocation pressure difference maximum in log10(hPa) ... for AMV with given pressures
tim_max_str 		= sys.argv[10]    	# Colloation time difference maximum in minutes
hgt_max_str 		= sys.argv[11]   	# Collocation height difference maximum in km ... vertical distance ... for AIRCRAFT with height data
dataset_names		= sys.argv[12]		# Names of all datasets to use for collocation, separated by comma "," delimiter.
archive_parent		= sys.argv[13]		# Archive parent path: path where home archive directory is located
n_max			= int(sys.argv[14])	# Maximum number of matches allowed per data point
nproc			= int(sys.argv[15])	# Number of processors to use during parallelization
AMV_center		= sys.argv[16]		# Abbreviation of wind-producing center

	#strip quotes from string inputs
qc_flags 	= qc_flags.strip('\"')
pct_str		= pct_str.strip('\"')
dst_max_str	= dst_max_str.strip('\"')
prs_max_str	= prs_max_str.strip('\"')
tim_max_str	= tim_max_str.strip('\"')
hgt_max_str	= hgt_max_str.strip('\"')
dataset_names 	= dataset_names.strip('\"')

print("CHECK USER INPUT ARGUMENTS:")
print("... Aeolus wind type = "+str(aeolus_wind_type))
print("... qc flags = "+str(qc_flags))
print("... AMV QI(s) = "+str(pct_str))
print("... AMV QI choice(s) = "+str(qi_choice_str))
print("... date IN = "+str(dateIN))
print("... aeolus dataset type = "+str(aeolus_dset_type_str))
print("... output path = "+str(output_path))
print("... DIST = "+str(dst_max_str))
print("... PRES = "+str(prs_max_str))
print("... TIME = "+str(tim_max_str))
print("... HGT  = "+str(hgt_max_str))
print("... all dataset names = "+str(dataset_names))
print("... archive parent path = "+str(archive_parent))

#=============================================
# Initial Setup

	#-------------------------------------
	# Get Input Dataset Names
dependent_names = []

dset_array = dataset_names.split(",")
iname=0
for name in range(len(dset_array)):
  iname=iname+1
  if name == 0:
    driver_name = dset_array[name]
  else:
    dependent_names.append(dset_array[name])
print("DRIVER dataset name: "+str(driver_name))
for x in range(len(dependent_names)):
  print("DEPENDENT dataset name(s): "+str(dependent_names[x]))

	# Number of dependent datasets
ndependent_datasets = np.size(dependent_names)

	#-------------------------------------
	# Get Input Dataset QC Flags
qc_dependent_flag = []

qc_array = qc_flags.split(",")
iqc=0
for qc in qc_array:
  iqc = iqc + 1
  if iqc == 1:
    qc_drv_flag = int(qc)
  elif iqc > 1:
    qc_dependent_flag.append(int(qc))

if iname!=iqc:
  print("ERROR: number of datasets does not equal number of QC dataset flags!")
  sys.exit()

	#-------------------------------------
        # Get AMV QI Value(s)

		#split QI string into list by delimiter ","
pct_array       = pct_str.split(",")
qi_choice_array = pct_str.split(",")

	#if driver is an AMV dataset
if ((driver_name.find('AMV') != -1) or (driver_name.find('amv') != -1)):
  pct_drv       = float(pct_array[0])
  qi_choice_drv = qi_choice_array
  namv=1	#indicates which QI value to assign to each dependent dataset (see next for-loop)
else:
  pct_drv       = np.nan
  qi_choice_drv = str(np.nan)
  namv=0
	#if 1+ dependent datasets is an AMV dataset,
	#  assign AMV QI value to array space of corresponding AMV dataset(s)

		#create 'pct' array (same dimension size as 'dependent_names'
pct       = np.nan * np.ones(ndependent_datasets, dtype = float)
qi_choice = np.nan * np.ones(ndependent_datasets, dtype = float)

for i in range(ndependent_datasets):
  if (dependent_names[i].find('AMV') != -1) or (dependent_names[i].find('amv') != -1):
    pct[i]       = pct_array[namv]
    qi_choice[i] = qi_choice_array[namv]
    namv = namv + 1

	#-------------------------------------
        # Get Collocation Criteria Value(s)

		#split string into list by delimiter ","
tdst = dst_max_str.split(",")
tprs = prs_max_str.split(",")
ttim = tim_max_str.split(",")
thgt = hgt_max_str.split(",")

dst_max = np.nan * np.ones(ndependent_datasets, dtype = float)
prs_max = np.nan * np.ones(ndependent_datasets, dtype = float)
tim_max = np.nan * np.ones(ndependent_datasets, dtype = float)
hgt_max = np.nan * np.ones(ndependent_datasets, dtype = float)
for i in range(ndependent_datasets):
  dst_max[i] = tdst[i]
  prs_max[i] = tprs[i]
  tim_max[i] = ttim[i]
  hgt_max[i] = thgt[i]

	#-------------------------------------
	# Set list of collocation criteria (for output NetCDF metadata)

colloc_list = "max distance = "+str(dst_max)+" km, max pressure difference = "+str(prs_max)+" log10(hPa), max time difference = "+str(tim_max)+" minutes, max height difference = "+str(hgt_max)+" km"
print("Collocation Criteria: "+str(colloc_list))

#=============================================
# Set QC booleans based on raw user input

	#-------------------------------------
	# Driver dataset
if qc_drv_flag == 1:
  bool_drv_qc = True
elif qc_drv_flag == 0:
  bool_drv_qc = False
  
  	#-------------------------------------
  	# Dependent datasets
bool_dset_qc = []
for i in range(ndependent_datasets):
  if qc_dependent_flag[i] == 1:
    bool_dset_qc.append(True)
  elif qc_dependent_flag[i] == 0:
    bool_dset_qc.append(False)
  
#=============================================
# Full date arrays

mmARR 		= [ "01","02","03","04","05","06","07","08","09","10","11","12" ]
ddARRend 	= [ "31","28","31","30","31","30","31","31","30","31","30","31" ]
ddARR 		= [ "01","02","03","04","05","06","07","08","09","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24","25","26","27","28","29","30","31" ]
hhARR		= [ "00","06","12","18" ]

#=============================================
# Find date before (B4) and date after (A) current date

print("----- DATE being processed = "+dateIN)

yyyy = dateIN[0:4]
mm   = dateIN[4:6]
dd   = dateIN[6:8]
hour = dateIN[8:10]

iyy = int(yyyy)
imm = int(mm)
idd = int(dd)

	#`````````````````````````````````````
	# Date BEFORE
if dd == "01":
  if (mm == "01") | (mm == "02") | (mm == "04") | (mm == "06") | (mm == "08") | (mm == "09") | (mm == "11"):
    ddB4 = "31"
    if mm == "01":
      yyB4 = str(iyy-1)
      mmB4 = "12"
    else:
      yyB4 = yyyy
      immB4 = mmARR.index(mm)-1
      mmB4  = mmARR[immB4] 
  elif mm == "03":
    yyB4  = yyyy
    immB4 = mmARR.index(mm)-1
    mmB4  = mmARR[immB4]
    ddB4  = "28"
    if (iyy % 4) == 0:
      ddB4 = "29"
  elif (mm == "05") | (mm == "07") | (mm == "10") | (mm == "12"):
    yyB4  = yyyy
    immB4 = mmARR.index(mm)-1
    mmB4  = mmARR[immB4]
    ddB4  = "30"    
else:
  yyB4  = yyyy
  mmB4  = mm
  iddB4	= ddARR.index(dd)-1
  ddB4	= ddARR[iddB4]
  
  	#`````````````````````````````````````
  	# Date AFTER
if (mm == "02" and (dd == "29" or (dd == "28" and iyy%4 != 0))) |\
   ((mm == "04" or mm == "06" or mm == "09" or mm == "11") and dd == "30") |\
   ((mm == "01" or mm == "03" or mm == "05" or mm == "07" or mm == "08" or mm == "10") and dd == "31"):
  yyA  = yyyy
  immA = mmARR.index(mm)+1
  mmA  = mmARR[immA]
  ddA  = "01"
elif mm == "12" and dd == "31":
  yyA  = str(iyy+1)
  mmA  = "01"
  ddA  = "01"
else:
  yyA  = yyyy
  mmA  = mm
  iddA = ddARR.index(dd)+1
  ddA  = ddARR[iddA]

 	#`````````````````````````````````````
yyyymmdd = yyyy+mm+dd		#Current date
dateB4   = yyB4+mmB4+ddB4	#Day before current date
dateA    = yyA+mmA+ddA		#Day after current date
 
#print("Current date to be processed = "+yyyymmdd)
print("...Date before = "+dateB4) 
print("...Date after = "+dateA)

yyyymmddhh = yyyymmdd+hour

#=============================================
# Define Datasets

idxs = [0]	# this is a bogus variable that is only used for plotting (not in this script)

	#---------------------------------------------
	# Define and load DRIVER dataset
	#---------------------------------------------
print("Define driver dataset")

if driver_name=='Aeolus':
  drv_lat,drv_lon,drv_prs,drv_hgt,drv_yr,drv_mm,drv_dy,drv_hr,drv_mn,indexesD,qc_drv_list,drv_src,drv_err,drv_len,drv_spd,drv_dir = read_aeolus(archive_parent,yyyymmddhh,dateB4,dateA,aeolus_dset_type_str,aeolus_wind_type,bool_drv_qc,"drv",runtype,tim_max_str,idxs)
elif driver_name=='Aircraft':
  drv_lat,drv_lon,drv_yr,drv_mm,drv_dy,drv_hr,drv_mn,drv_hgt,drv_prs,indexesD,qc_drv_list,drv_src,drv_spd,drv_dir = read_aircraft(archive_parent,yyyymmddhh,dateB4,dateA,bool_drv_qc,"drv",runtype,tim_max_str,idxs)
elif driver_name=='AMV_NCEP':
  drv_lat,drv_lon,drv_yr,drv_mm,drv_dy,drv_hr,drv_mn,drv_hgt,drv_prs,indexesD,qc_drv_list,drv_src,drv_satname,drv_wcm,drv_ham,drv_spd,drv_dir = read_amv_ncep(archive_parent,yyyymmddhh,dateB4,dateA,bool_drv_qc,pct_drv,qi_choice_drv,"drv",runtype,tim_max_str,idxs)
elif driver_name=='Loon':
  drv_lat,drv_lon,drv_yr,drv_mm,drv_dy,drv_hr,drv_mn,drv_hgt,drv_prs,indexesD,qc_drv_list,drv_src,drv_hgt,drv_azm,drv_elv,drv_spd,drv_dir = read_loon(archive_parent,yyyymmddhh,dateB4,dateA,bool_drv_qc,"drv",runtype,tim_max_str,idxs)
elif driver_name=='Sonde':
  drv_lat,drv_lon,drv_yr,drv_mm,drv_dy,drv_hr,drv_mn,drv_hgt,drv_prs,indexesD,qc_drv_list,drv_src,nsondes,nlevels,ngroups,drv_spd,drv_dir = read_sonde(archive_parent,yyyymmddhh,dateB4,dateA,bool_drv_qc,"drv",runtype,tim_max_str,idxs)
else:
  print("ERROR: Driver dataset not defined!")
  sys.exit()

	# longitude check: set range to 0=360
for iP in range(np.size(drv_lon)):
  if drv_lon[iP] < 0.0:
    drv_lon[iP] = drv_lon[iP] + 360.

        # pressure check: set pressure to NaN if p<=0
for iP in range(np.size(drv_prs)):
  if drv_prs[iP] <= 0.0:
    drv_prs[iP] = np.nan

        # height check: set height to NaN if hgt<0
for iP in range(np.size(drv_hgt)):
  if drv_hgt[iP] < 0.0:
    drv_hgt[iP] = np.nan

	#---------------------------------------------
	# Define dependent datasets and append to each other to create total matching array for collocation
	#---------------------------------------------
print("Define dependent datasets")

dset_len = np.nan * np.ones(ndependent_datasets, dtype = float)

qc_dset_list = []
dset_src     = []

for i in range(ndependent_datasets):
  print("i="+str(i)+" | "+dependent_names[i])
	# read dependent dataset
  if dependent_names[i]=='Aeolus':
    t_lat,t_lon,t_prs,t_hgt,t_yr,t_mm,t_dy,t_hr,t_mn,t_indexes,t_qc_dset_list,t_src,t_err,t_len,t_spd,t_dir = read_aeolus(archive_parent,yyyymmddhh,dateB4,dateA,aeolus_dset_type_str,aeolus_wind_type,bool_dset_qc[i],"dep",runtype,tim_max_str,idxs)
  elif dependent_names[i]=='Aircraft':
    t_lat,t_lon,t_yr,t_mm,t_dy,t_hr,t_mn,t_hgt,t_prs,t_indexes,t_qc_dset_list,t_src,t_spd,t_dir = read_aircraft(archive_parent,yyyymmddhh,dateB4,dateA,bool_dset_qc[i],"dep",runtype,tim_max_str,idxs)
  elif dependent_names[i]=='AMV_NCEP':
    t_lat,t_lon,t_yr,t_mm,t_dy,t_hr,t_mn,t_hgt,t_prs,t_indexes,t_qc_dset_list,t_src,t_satname,t_wcm,t_ham,t_spd,t_dir = read_amv_ncep(archive_parent,yyyymmddhh,dateB4,dateA,bool_dset_qc[i],pct[i],qi_choice[i],"dep",runtype,tim_max_str,idxs)
  elif dependent_names[i]=='Loon':
    t_lat,t_lon,t_yr,t_mm,t_dy,t_hr,t_mn,t_hgt,t_prs,t_indexes,t_qc_dset_list,t_src,t_azm,t_elv,t_spd,t_dir = read_loon(archive_parent,yyyymmddhh,dateB4,dateA,bool_dset_qc[i],"dep",runtype,tim_max_str,idxs)
  elif dependent_names[i]=='Sonde':
	# RADIOSONDE data are 2D (except for lat, lon, dates/times): [nsondes, nlevels]
	# 	nsondes = number of sondes
	# 	nlevels = number of vertical levels per sonde
	# nsondes, nlevels are used to find collocated indices that pertain to original data file
    t_lat,t_lon,t_yr,t_mm,t_dy,t_hr,t_mn,t_hgt,t_prs,t_indexes,t_qc_dset_list,t_src,nsondes,nlevels,ngroups,t_spd,t_dir = read_sonde(archive_parent,yyyymmddhh,dateB4,dateA,bool_dset_qc[i],"dep",runtype,tim_max_str,idxs)
  else:
    print("ERROR: Dependent dataset "+str(i)+" ("+dependent_names[i]+") cannot be read by this program! Please add function to read_input_datasets module and try again.")
    sys.exit()

	# longitude check: set range to 0=360
  for iP in range(np.size(t_lon)):
    if t_lon[iP] < 0.0:
      t_lon[iP] = t_lon[iP] + 360.

	# pressure check: set pressure to NaN if p<=0
  for iP in range(np.size(t_prs)):
    if t_prs[iP] <= 0.0:
      t_prs[iP] = np.nan

	# height check: set height to NaN if hgt<0
  for iP in range(np.size(t_hgt)):
    if t_hgt[iP] < 0.0:
      t_hgt[iP] = np.nan

	# assign integers to 'usePH' to determine if should use height or pressure as vertical collocation criterion
	#	0 = use pressure (P)
	#	1 = use height (H)
  if np.isnan(drv_prs).all() or np.isnan(t_prs).all():
		# use height
    tmp = 1
  elif (np.isnan(drv_prs).all() and np.isnan(drv_hgt).all()) or (np.isnan(t_prs).all() and np.isnan(t_hgt).all()):
    if np.isnan(drv_prs).all() and np.isnan(drv_hgt).all():
      print("ERROR: Driver dataset has neither pressure nor height variable. One of these must be present to continue.")
    if np.isnan(t_prs).all() and np.isnan(t_hgt).all():
      print("ERROR: Dependent dataset "+str(dependent_names[i])+" has neither pressure nor height variable. One of these must be present to continue.")
    sys.exit()
  else:
    tmp = 0	# use pressure

  if i==0:
	# initialize all_match* arrays with first dependent dataset listed
    all_match_lat = t_lat
    all_match_lon = t_lon
    all_match_prs = t_prs     
    all_match_hgt = t_hgt    
    all_match_yr  = t_yr 
    all_match_mm  = t_mm 
    all_match_dy  = t_dy 
    all_match_hr  = t_hr 
    all_match_mn  = t_mn 
    
    	# collocation criteria to be applied to each obs
    all_dst_max    = np.nan * np.ones(np.size(t_prs),dtype=float)
    all_prs_max    = np.nan * np.ones(np.size(t_prs),dtype=float)
    all_tim_max    = np.nan * np.ones(np.size(t_prs),dtype=float)
    all_hgt_max    = np.nan * np.ones(np.size(t_prs),dtype=float)

    all_dst_max[:] = dst_max[i]
    all_prs_max[:] = prs_max[i]
    all_tim_max[:] = tim_max[i]
    all_hgt_max[:] = hgt_max[i]

	# indices
    indexes       = t_indexes
    qc_dset_list  = [t_qc_dset_list]
    dset_src      = [t_src]

    dset_ind	  = np.ones(np.size(t_prs),dtype=int)
    dset_ind[:]	  = 0

	# flag determining if pressure or height should be used for collocation
    usePH         = np.nan * np.ones(np.size(t_prs),dtype=int)
    usePH[:]      = tmp

  else:
  	# append remaining dependent datasets to all_match* arrays
    all_match_lat = np.append(all_match_lat,t_lat,axis=0)
    all_match_lon = np.append(all_match_lon,t_lon,axis=0)
    all_match_prs = np.append(all_match_prs,t_prs,axis=0)     
    all_match_hgt = np.append(all_match_hgt,t_hgt,axis=0)    
    all_match_yr  = np.append(all_match_yr ,t_yr ,axis=0)
    all_match_mm  = np.append(all_match_mm ,t_mm ,axis=0)
    all_match_dy  = np.append(all_match_dy ,t_dy ,axis=0)
    all_match_hr  = np.append(all_match_hr ,t_hr ,axis=0)
    all_match_mn  = np.append(all_match_mn ,t_mn ,axis=0)

	# collocation criteria to be applied to each obs
    tall_dst_max  = np.nan * np.ones(np.size(t_prs),dtype=float)
    tall_prs_max  = np.nan * np.ones(np.size(t_prs),dtype=float)
    tall_tim_max  = np.nan * np.ones(np.size(t_prs),dtype=float)
    tall_hgt_max  = np.nan * np.ones(np.size(t_prs),dtype=float)
    
    tall_dst_max[:] = dst_max[i]
    tall_prs_max[:] = prs_max[i]
    tall_tim_max[:] = tim_max[i]
    tall_hgt_max[:] = hgt_max[i]
    
    all_dst_max   = np.append(all_dst_max, tall_dst_max, axis=0)
    all_prs_max   = np.append(all_prs_max, tall_prs_max, axis=0)
    all_tim_max   = np.append(all_tim_max, tall_tim_max, axis=0)
    all_hgt_max   = np.append(all_hgt_max, tall_hgt_max, axis=0)

	# indices
    indexes 	  = np.append(indexes,t_indexes,axis=0)
    qc_dset_list.append(t_qc_dset_list)
    dset_src.append(t_src)

    tdset_ind	  = np.ones(np.size(t_prs),dtype=int)
    tdset_ind[:]  = int(i)
    dset_ind	  = np.append(dset_ind,tdset_ind,axis=0)

	# flag determining if pressure or height should be used for collocation
    tusePH        = np.nan * np.ones(np.size(t_prs),dtype=int)
    tusePH[:]     = tmp
    usePH         = np.append(usePH,tusePH,axis=0)

  	# length of each datasets
  dset_len[i] = np.size(t_lat)

  print("Size of dataset "+str(i)+"/"+str(ndependent_datasets-1)+" ("+str(dependent_names[i])+") = "+str(dset_len[i])+" = "+str(np.size(usePH)))

# Total size of all_match* arrays
all_match_size = np.size(all_match_prs)
  	
#=============================================
# Create output directory
#	If the program made it this far without exiting, then all datasets are confirmed to exist.
#	Thus, this program is able to read all datasets, and index files will be generated.

print("Create archive directory")
os.system('mkdir -m 775 -p '+output_path)

#=============================================
# Set output NetCDF filename

output_name = 'index.'+str(dateIN)

	# add driver name
if driver_name == 'Aeolus':
  if aeolus_dset_type_str == 'orig':
    drv_dset_type_str = 'Orig'
  else:
    drv_dset_type_str = 'Reproc'+str(aeolus_dset_type_str)
  outname_driver = driver_name+'_'+aeolus_wind_type+'_'+drv_dset_type_str
elif driver_name != 'Aeolus':
  if driver_name == 'AMV_4th_Int':
    outname_driver = new_driver_name 
  else:
    outname_driver = driver_name

output_name += '.drv_'+outname_driver

		# add 'QC' if qc_flag=True
if bool_drv_qc:
  output_name+='__QC'
else:
  output_name+='__NoQC'

	# add dependent dataset names
for ndsets in range(ndependent_datasets):
  if dependent_names[ndsets] == 'Aeolus':
    if aeolus_dset_type_str == 'orig':
      a_dset_type_str = ''
    else:
      a_dset_type_str = 'Reproc'+str(aeolus_dset_type_str)
    outname_dep = str(dependent_names[ndsets])+'_'+aeolus_wind_type+'_'+a_dset_type_str
  elif dependent_names[ndsets] != 'Aeolus':
    if dependent_names[ndsets] == 'AMV_4th_Int':
      outname_dep = str(dependent_names[ndsets])+'_'+amv4th_type_str
    else:
      outname_dep = str(dependent_names[ndsets])

  output_name+='.dset'+str(ndsets+1)+'_'+outname_dep

		# add 'QC' if qc_flag=True
  if bool_dset_qc[ndsets]:
    output_name+='__QC'
  else:
    output_name+='__NoQC'

# Full output path and filename
output_string = output_path+output_name

#=============================================
# Run obs_match_3d(): Match amv_match* data to amv_drv* data

if __name__ ==  '__main__':
    print("CALL obs_match_3d")
    a_match_idx,a_drv_idx,GCD_match,DP_match,DPlog_match,HT_match,DT_match = obs_match_3d(
                                         all_match_lat ,
                                         all_match_lon ,
                                         all_match_prs ,
                                         all_match_yr  ,
                                         all_match_mm  ,
                                         all_match_dy  ,
                                         all_match_hr  ,
                                         all_match_mn  ,
                                         drv_lat ,
                                         drv_lon ,
                                         drv_prs ,
                                         drv_yr  ,
                                         drv_mm  ,
                                         drv_dy  ,
                                         drv_hr  ,
                                         drv_mn  ,
                                         all_dst_max  ,
                                         all_prs_max  ,
                                         all_tim_max  ,
                                         n_max    ,
                                         nproc    ,
					 all_match_hgt ,
					 drv_hgt ,
					 all_hgt_max ,
                                         usePH
                                        )
    print("Number of matches = "+str(np.size(GCD_match))) # str(np.shape(GCD_match)))
    if np.size(GCD_match)==0:
      print("Zero matches found. EXIT PROGRAM")
      sys.exit()		# if no matches between DRIVER and DEPENDENT, exit MAIN program

    # Get match indices from original files (before QC) - Driver dataset (AEOLUS)
    tmp = a_drv_idx
    a_drv_idx_out = indexesD[tmp]
    del tmp

    # Get indices of dependent dataset indicators from matched indices
    if np.size(dset_ind) == np.size(drv_prs):
      tmp = a_drv_idx
    elif np.size(dset_ind) == np.size(all_match_prs):
      tmp = a_match_idx
    dset_ind_match = dset_ind[tmp]
    del tmp

    #==============================================================
    # Get variables to output to NetCDF

    fill = -999
    
	#index matches, for output NetCDF
    var_vals  = []
    	#variable names, for output NetCDF
    var_names = []
        #variable short names, for output NetCDF
    var_shortnames = []
	#variable long names, for output NetCDF
    var_longnames = []
    	#variable units, for output NetCDF
    var_units = []
    	#variable missing value, for output NetCDF
    var_miss = []
	#wind type (for Aeolus only)
    var_aeol_wind_type = []
        #dataset type (for Aeolus only)
    var_aeol_dset_type = []
    	#extra notes
    var_notes = []

	#for dataset variables only
		#path
    var_path   = []
    var_qcflag = []
    var_qclist = []
    var_qcnotes = []
    
    	# Add dataset names as variables
		# DRIVER
    var_vals  = [driver_name]	#data for variable
    var_names = ['drv']		#name of variable
    var_shortnames = [str(outname_driver)]
    var_longnames = ['DRIVER dataset with which all dependent datasets are collocated']
    var_units = [str(fill)]
    var_miss  = [str(fill)]
    var_path  = [drv_src]
    var_qcflag = [qc_drv_flag]
    var_qclist = [qc_drv_list]
    var_qcnotes = ['If QC_flag=0, QC is not applied; if QC_flag=1, QC listed in QC_list is applied.']
    if (driver_name.find('Aeolus') != -1):
      var_aeol_wind_type = [aeolus_wind_type]
      var_aeol_dset_type = [drv_dset_type_str]
    else:
      var_aeol_wind_type = [str(fill)]    
      var_aeol_dset_type = [str(fill)]    

    		# DEPENDENT(S)
    for x in range(len(dependent_names)):
      var_vals  += [dependent_names[x]]
      var_names += ['dset'+str(x+1)]
      var_shortnames += [str(dependent_names[x])]
      var_longnames += ['DEPENDENT DATASET '+str(x+1)]
      var_units += [str(fill)]
      var_miss  += [str(fill)]
      var_path  += [dset_src[x]]
      var_qcflag += [qc_dependent_flag[x]]
      var_qclist += [qc_dset_list[x]]
      var_qcnotes += ['If QC_flag=0, QC is not applied; if QC_flag=1, QC listed in QC_list is applied.']
      if (dependent_names[x].find('Aeolus') != -1):
        var_aeol_wind_type += [aeolus_wind_type]
        var_aeol_dset_type += [a_dset_type_str]
      else:
        var_aeol_wind_type += [str(fill)]	
        var_aeol_dset_type += [str(fill)]	

	# Add number of dependent datasets as variable
    var_vals  += [ndependent_datasets]
    var_names += ['ndset']
    var_longnames += ['Number of DEPENDENT DATASETs collocated with DRIVER']
    var_units += [str(fill)]
    var_miss  += [str(fill)]

	# Add collocation criteria as variables
    var_vals  += [tim_max]
    var_names += ['time_max']
    var_longnames += ['Collocation criterion: Maximum absolute time difference']
    var_units += ['minutes']
    var_notes += ['each value corresponds to each dependent dataset, in order (e.g., *_max at dimension 0 corresponds to dset1, *_max at dimension 1 corresponds to dset2, etc.)']
    var_miss  += [str(fill)]
    
    var_vals  += [prs_max]
    var_names += ['pres_max']
    var_longnames += ['Collocation criterion: Maximum absolute log10(pressure) difference']
    var_units += ['hPa']
    var_notes += ['each value corresponds to each dependent dataset, in order (e.g., *_max at dimension 0 corresponds to dset1, *_max at dimension 1 corresponds to dset2, etc.)']
    var_miss  += [str(fill)]
    
    var_vals  += [hgt_max]
    var_names += ['hgt_max']
    var_longnames += ['Collocation criterion: Maximum absolute height difference']
    var_units += ['km']
    var_notes += ['each value corresponds to each dependent dataset, in order (e.g., *_max at dimension 0 corresponds to dset1, *_max at dimension 1 corresponds to dset2, etc.)']
    var_miss  += [str(fill)]
    
    var_vals  += [dst_max]
    var_names += ['dist_max']
    var_longnames += ['Collocation criterion: Maximum great circle distance']
    var_units += ['km']
    var_notes += ['each value corresponds to each dependent dataset, in order (e.g., *_max at dimension 0 corresponds to dset1, *_max at dimension 1 corresponds to dset2, etc.)']
    var_miss  += [str(fill)]
    
        # Split a_match_idx and a_drv_idx into two lists each, anything with an index value greater than dset1_len
    	# is in the second list, and corresponds to dset2 data, etc.
    print("Find matched indices from source data files")

    drv_size = np.nan * np.ones(ndependent_datasets,dtype=int)	# used to determine whether output NetCDF file should be created: if size(drv_flg)=0, do not create.
    
    for ndsets in range(ndependent_datasets):
      tmpmatch = []
      tmpdrv   = []
      tmpDT    = []
      tmpDP    = []
      tmpDPlog = []
      tmpHT    = []
      tmpGCD   = []

      PorH = "N"		# indicates whether pressure or height collocation criterion was used in collocation; initialize "N" for None
      for i in range(len(a_match_idx)):
        idrv 	 = a_drv_idx_out[i]         #driver dataset index matches
        itmp 	 = a_match_idx[i]           #dependent dataset index matches
        idx_orig = indexes[itmp]
	
        tDT 	 = DT_match[i]
        tDP	 = DP_match[i]
        tDPlog	 = DPlog_match[i]
        tHT	 = HT_match[i]
        tGCD	 = GCD_match[i] 

        if dset_ind_match[i] == int(ndsets):
                #matching indices per dependent dataset
          tmpmatch.append(idx_orig)
          tmpdrv.append(idrv)
	  	#collocation dimensions for matching pairs
          tmpDT.append(tDT)
          tmpGCD.append(tGCD)
          if not np.isnan(tDP):
            PorH = "P"
            tmpDP.append(tDP)
            tmpDPlog.append(tDPlog)
          elif not np.isnan(tHT) and np.isnan(tDP):
            PorH = "H"
            tmpHT.append(tHT)

		#matching driver dataset indices
      var_vals  += [tmpdrv]					#appends driver indices that match dependent dataset to 'var_vals'
      var_names += ['idx_drv_dset'+str(ndsets+1)]		#appends dependent dataset index variable name to 'var_names'
      var_longnames += ['Indices of DRIVER that match DATASET '+str(ndsets+1)] 
      var_units += [str(fill)]
      var_miss  += [str(fill)]
 
		#matching dependent dataset indices
      var_vals  += [tmpmatch]					#appends dependent dataset indices that match driver to 'var_vals'
      var_names += ['idx_dset'+str(ndsets+1)]			#appends dependent dataset index variable name to 'var_names'
      var_longnames += ['Indices of DATASET '+str(ndsets+1)+' that match DRIVER'] 
      var_units += [str(fill)]
      var_miss  += [str(fill)]

      	#----------------------------
      	# Fill drv_size (for check later)
      
      drv_size[ndsets] = np.size(tmpdrv)
      #print("... ndsets = "+str(ndsets)+" | DRV_SIZE = "+str(drv_size[ndsets]))

      	#----------------------------	
      	# Collocation dimensions (GCD, DP, DT) for matched pairs

      var_vals  += [tmpDT]
      var_names += ['DT_match_drv_dset'+str(ndsets+1)]
      var_longnames += ['Time difference for each collocated pair between DRIVER and DATASET '+str(ndsets+1)]
      var_units += ['minutes']
      var_miss  += [str(fill)]
	
      if PorH == "P":
        var_vals  += [tmpDP]
        var_names += ['DP_match_drv_dset'+str(ndsets+1)]
        var_longnames += ['Difference in pressure for each collocated pair between DRIVER and DATASET '+str(ndsets+1)]
        var_units += ['hPa']
        var_miss  += [str(fill)]

        var_vals  += [tmpDPlog]
        var_names += ['DPlog_match_drv_dset'+str(ndsets+1)]
        var_longnames += ['Difference in log10(pressure) for each collocated pair between DRIVER and DATASET '+str(ndsets+1)]
        var_units += ['log10(hPa)']
        var_miss  += [str(fill)]
      elif PorH == "H":
        var_vals  += [tmpHT]
        var_names += ['HT_match_drv_dset'+str(ndsets+1)]
        var_longnames += ['Height difference for each collocated pair between DRIVER and DATASET '+str(ndsets+1)]
        var_units += ['km']
        var_miss  += [str(fill)]
	
      var_vals  += [tmpGCD]
      var_names += ['GCD_match_drv_dset'+str(ndsets+1)]
      var_longnames += ['Great circle distance for each collocated pair between DRIVER and DATASET '+str(ndsets+1)]
      var_units += ['km']
      var_miss  += [str(fill)]
		
	#------------------------------------------
	# This code block is for profile datasets with 2D variables (e.g., Sonde)

      if driver_name == "Sonde" or dependent_names[ndsets] == "Sonde":
        var_vals  += [ngroups]  #append number of sondes per group to NetCDF output var
        var_names += ['ntypes_dset'+str(ndsets+1)]
        var_longnames += ['Number of Sonde types']
        var_units += [str(fill)]
        var_miss  += [str(fill)]

        var_vals  += [nsondes]	#append number of sondes per group to NetCDF output var
        var_names += ['nprofiles_dset'+str(ndsets+1)]
        var_longnames += ['Number of Sondes per type'] 
        var_units += [str(fill)]
        var_miss  += [str(fill)]

        var_vals  += [nlevels]	#append number of levels per sonde per group to NetCDF output var
        var_names += ['nlevels_dset'+str(ndsets+1)]
        var_longnames += ['Number of max vertical levels for each Sonde per type'] 
        var_units += [str(fill)]
        var_miss  += [str(fill)]

		#NOTE: plotting script will use these variables to extract indices of Sonde obs collocated with DRIVER

	#------------------------------------------

    	# matches exist! create output file
    print("Matches exist for "+str(dateIN)+"! Create output index file:")
	    
    	# Send index values to netCDF file
    nc_out_filenm = output_string+".nc4"
    print("... OUTPATH = "+output_path)
    print("... OUTFILE = "+nc_out_filenm)
    
      #==============================================================
      # Output variables and metadata to NetCDF

    match_lists_to_netCDF(var_names,var_vals,nc_out_filenm,var_longnames,var_units,var_miss,var_path,var_qcflag,var_qclist,var_qcnotes,var_shortnames,var_aeol_wind_type,var_aeol_dset_type)
    
      #==============================================================

print("========== END MAIN PROGRAM ==========")
print("======================================")    
#########################################################################################################
# END SCRIPT
#########################################################################################################
