###################################################################################################################
###################################################################################################################
# MAIN SCRIPT FOR Plotting Matched Winds from Multiple Datasets
###################################################################################################################
###################################################################################################################
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
# History: 
#	2021		B. Hoover, D. Santek    Created/developed program.
#	2021-10-25	K.E. Lukens		Added aircraft, AMV input configs. Added quality controls (QC).
#	2021-11-22	K.E. Lukens		Added date before consideration.
#	2021-11-23	K.E. Lukens		Created QC module and added relevant calls.
#	2022-03-08	K.E. Lukens		Added metadata to output files.
#	2022-10-05	K.E. Lukens		Updated script with additional plotting capabilities.
#	2022-10-11	K.E. Lukens		Added capability to ingest multiple collocation files.
#       2023-04-14      K.E. Lukens             Finalized for upload to SAWC archive.
#	2023-04-21	K.E. Lukens		Bug fix: Now uses correct input path for archived index files. str(dset_name) now used with .find function
#	2023-07-15	D. Huber, K.E. Lukens	Optimized program to improve performance. Added new Pressure/Height vs Longitude plot type. Various bug fixes.
#
# Output: Figures (.png) and Videos (.gif)
#
###################################################################################################################
# Import python modules
###################################################################################################################
print("***** Import Python Modules *****")

import os
os.environ['OPENBLAS_NUM_THREADS'] = '1'

from os.path import exists

from datetime import datetime

import sys
import math
import statistics as stats
import numpy as np
import warnings

import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import cartopy.crs as ccrs
import cartopy.mpl.ticker as cticker
from scipy.interpolate import griddata
from scipy.interpolate import interpn
from scipy.interpolate import RegularGridInterpolator

from netCDF4 import Dataset

from read_data import read_index_file
from read_data import read_aeolus
from read_data import read_aircraft
from read_data import read_amv_ncep
#from read_data import read_loon
from read_data import read_sonde

from quality_controls import qc_winds

from tools_analysis import compute_hlos
from tools_analysis import to_nan
from tools_analysis import superob_matches
from tools_analysis import prs_to_hgt
from tools_analysis import var_to_2d
from tools_analysis import var_to_2d_counts

from tools_plotting import density_scatter
from tools_plotting import hist_diffs
from tools_plotting import map_contour2d_ce
from tools_plotting import map_points2d_ce
from tools_plotting import map_locations_ce
from tools_plotting import map_locations_ortho
from tools_plotting import map_locations_ortho_rotate
from tools_plotting import map_3d_profile
from tools_plotting import scatter_matches
from tools_plotting import time_series
from tools_plotting import z_distr
from tools_plotting import stats_windbins
from tools_plotting import contour2d
from tools_plotting import contour2d_orthomap
from tools_plotting import contour2d_orthomap_rotate

###################################################################################################################
print("========================================")
print("========== BEGIN MAIN PROGRAM ==========")

#=============================================
# Set global parameters

fill = -999.0

runtype = "plot"			# type of program coded here: PLOT code

	# set horizontal resolution
#latlonstride = 0.5				# 0.5-degree resolution
latlonstride = 1.0				# 1.0-degree resolution
flats = list(np.arange(-90,91,latlonstride))
flons = list(np.arange(0,361,latlonstride)) 
flons180 = list(np.arange(-180,181,latlonstride)) 

	# set vertical resolution
plevstride = 25.0
fpres = np.arange(0,1050,plevstride)		# pressure array
hlevstride = 1.0
fhgts = list(np.arange(0,21,hlevstride))	# height array

	# convert scattered data to regular grid
xgrid,ygrid = np.meshgrid(flons,flats)

#=============================================
# Read raw user input from command line

aeolus_wind_type        = sys.argv[1]           # Aeolus wind type
dateSTART               = sys.argv[2]           # Start date: YYYYMMDDHH
dateEND			= sys.argv[3]		# End date: YYYYMMDDHH
input_path		= sys.argv[4]		# Input directory for index files
dependent_names_str	= sys.argv[5]		# Dependent dataset names
input_file_suffix	= sys.argv[6]		# Input index file suffixes
output_path             = sys.argv[7]           # Output (archive) directory
archive_parent		= sys.argv[8]		# Archive parent path: path where home archive directory is located
avgthin_choice		= int(sys.argv[9])	# Choice to super-ob, thin, or plot all matches

# temporary assignment
qi_choice = "NO_FC"

# strip quotes from string inputs
dependent_names_str = dependent_names_str.strip('\"')
input_file_suffix = input_file_suffix.strip('\"')

print("CHECK USER INPUT ARGUMENTS:")
print("... Aeolus wind type = "+str(aeolus_wind_type))
print("... date START = "+str(dateSTART))
print("... date END = "+str(dateEND))
print("... input path = "+str(input_path))
print("... dependent dataset names = "+str(dependent_names_str))
print("... input index file suffix = "+str(input_file_suffix))
print("... output path = "+str(output_path))
print("... archive parent path = "+str(archive_parent))
print("... avg/thin choice = "+str(avgthin_choice))

        #-------------------------------------
        # Get each index file suffix

# split string into list by delimiter ","
tsuffix = input_file_suffix.split(",")

dependent_names = dependent_names_str.split(",")

tndset = np.size(dependent_names)		# number of dependent datasets 
print("number of dependent datasets = "+str(tndset))

suffix = []
for i in range(tndset):
  suffix.append(tsuffix[i])

#=============================================
# Full date arrays

mmARR 		= [ "01","02","03","04","05","06","07","08","09","10","11","12" ]
ddARRend 	= [ "31","28","31","30","31","30","31","31","30","31","30","31" ]
ddARR 		= [ "01","02","03","04","05","06","07","08","09","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24","25","26","27","28","29","30","31" ]
hhARR		= [ "00","06","12","18" ]

#=============================================
# Find date before (B4) ... used for Aeolus input (driver dataset)

print("----- DATE range: "+dateSTART+"-"+dateEND)

yyyyS = dateSTART[0:4]
mmS   = dateSTART[4:6]
ddS   = dateSTART[6:8]
hhS   = dateSTART[8:10]
yyyyE = dateEND[0:4]
mmE   = dateEND[4:6]
ddE   = dateEND[6:8]
hhE   = dateEND[8:10]

print("START yyyy="+yyyyS+",mm="+mmS+",dd="+ddS+",hour="+hhS)
print("END   yyyy="+yyyyE+",mm="+mmE+",dd="+ddE+",hour="+hhE)

#---------------------------------------------------------------
# LOOP through dates

idx_file_vertstr = [[] for i in range(tndset)]

tim_max   = np.nan * np.ones(tndset, dtype=float)
prs_max   = np.nan * np.ones(tndset, dtype=float)
hgt_max   = np.nan * np.ones(tndset, dtype=float)
dst_max   = np.nan * np.ones(tndset, dtype=float)

dset_name   = [[] for i in range(tndset)]
DT_match    = [[] for i in range(tndset)]
GCD_match   = [[] for i in range(tndset)]
Vert_match  = [[] for i in range(tndset)]

Dyyyymmdd   = [[] for i in range(tndset)]
Dyyyymmddhh = [[] for i in range(tndset)]
drv_lat   = [[] for i in range(tndset)]
drv_lon   = [[] for i in range(tndset)]
drv_prs   = [[] for i in range(tndset)]
drv_hgt   = [[] for i in range(tndset)]
drv_spd   = [[] for i in range(tndset)]
drv_dir   = [[] for i in range(tndset)]
dset_lat  = [[] for i in range(tndset)]
dset_lon  = [[] for i in range(tndset)]
dset_prs  = [[] for i in range(tndset)]
dset_hgt  = [[] for i in range(tndset)]
dset_spd  = [[] for i in range(tndset)]
dset_dir  = [[] for i in range(tndset)]

drv_wcm   = [[] for i in range(tndset)]
dset_wcm  = [[] for i in range(tndset)]

yyyymmddARR = []
yyyymmddhhARR = []

yyyy = yyyyS
iyy  = int(yyyy)
iyyE = int(yyyyE)
while iyy <= int(yyyyE):
  if yyyyS==yyyyE:
    Simm = mmARR.index(mmS)
    Eimm = mmARR.index(mmE)
  elif yyyyS!=yyyyE:
    if iyy==int(yyyyS):
      Simm = mmARR.index(mmS)
      Eimm = 12-1
    elif iyy>int(yyyyS) and iyy<iyyE:
      Simm = 0
      Eimm = 12-1
    elif iyy==iyyE:
      Simm = 0
      Eimm = mmARR.index(mmE)

  imm=Simm
  while imm <= Eimm:
    mm = mmARR[imm]
    if mm=="02":
      if int(yyyy)%4==0:
        # leap year
        nDD = 29
      else:
        nDD = 28
    elif mm=="04" or mm=="06" or mm=="09" or mm=="11":
      nDD = 30
    else:
      nDD = 31
 
    if yyyyS==yyyyE:
      if mmS==mmE:
        Sidd = ddARR.index(ddS)
        Eidd = ddARR.index(ddE)
      else:
        if mm==mmS:
          Sidd = ddARR.index(ddS)
          Eidd = nDD-1
        elif mm==mmE:
          Sidd = 0
          Eidd = ddARR.index(ddE)
        else:
          Sidd = 0
          Eidd = nDD-1
    elif yyyyS!=yyyyE:
      if yyyy==yyyyS and mm==mmS:
        Sidd = ddARR.index(ddS)
        Eidd = nDD-1
      elif iyy==int(yyyyE) and mm==mmE:
        Sidd = 0
        Eidd = ddARR.index(ddE)
      else:
        Sidd = 0
        Eidd = nDD-1

    idd=Sidd
    while idd <= Eidd:
      dd = ddARR[idd]
      if yyyyS==yyyyE and mmS==mmE and ddS==ddE:
        Sihh = hhARR.index(hhS)
        Eihh = hhARR.index(hhE)
      else:
        if yyyy==yyyyS and mm==mmS and dd==ddS:
          Sihh = hhARR.index(hhS)
          Eihh = 4-1
        elif yyyy==yyyyE and mm==mmE and dd==ddE:
          Sihh = 0
          Eihh = hhARR.index(hhE)
        else:
          Sihh = 0
          Eihh = 4-1
   
      ihh=Sihh
      while ihh <= Eihh:
        hh = hhARR[ihh]
        
  		#`````````````````````````````````````
		# Date BEFORE
        if dd == "01":
          if (mm == "01") | (mm == "02") | (mm == "04") | (mm == "06") | (mm == "08") | (mm == "09") | (mm == "11"):
            ddB4 = "31"
            if mm == "01":
              yyB4 = str(iyy-1)
              mmB4 = "12"
              prevyy = iyy-1
            else:
              yyB4 = yyyy
              immB4 = mmARR.index(mm)-1
              mmB4  = mmARR[immB4] 
              prevyy = iyy
          elif mm == "03":
            yyB4  = yyyy
            immB4 = mmARR.index(mm)-1
            mmB4  = mmARR[immB4]
            ddB4  = "28"
            if (iyy % 4) == 0:
              ddB4 = "29"
            prevyy = iyy
          elif (mm == "05") | (mm == "07") | (mm == "10") | (mm == "12"):
            yyB4  = yyyy
            immB4 = mmARR.index(mm)-1
            mmB4  = mmARR[immB4]
            ddB4  = "30"    
            prevyy = iyy
        else:
          yyB4  = yyyy
          mmB4  = mm
          iddB4	= ddARR.index(dd)-1
          ddB4	= ddARR[iddB4]
          prevyy = iyy
 
 		#`````````````````````````````````````
  		# Date AFTER
        if (mm == "02" and (dd == "29" or (dd == "28" and iyy%4 != 0))) |\
           ((mm == "04" or mm == "06" or mm == "09" or mm == "11") and dd == "30") |\
           ((mm == "01" or mm == "03" or mm == "05" or mm == "07" or mm == "08" or mm == "10") and dd == "31"):
          yyA  = yyyy
          immA = mmARR.index(mm)+1
          mmA  = mmARR[immA]
          ddA  = "01"
          nextyy = iyy
        elif mm == "12" and dd == "31":
          yyA  = str(iyy+1)
          mmA  = "01"
          ddA  = "01"
          if hh == "18":
            nextyy = iyy+1
          else:
            nextyy = iyy
        else:
          yyA  = yyyy
          mmA  = mm
          iddA = ddARR.index(dd)+1
          ddA  = ddARR[iddA]
          nextyy = iyy
 
 	#`````````````````````````````````````
        yyyymmdd = str(iyy)+mm+dd		#Current date
        dateB4   = str(prevyy)+mmB4+ddB4	#Day before current date
        dateA    = str(nextyy)+mmA+ddA		#Day after current date
 
        print("Date to be processed in LOOP   = "+yyyymmdd)
        print("...Date before = "+dateB4) 
        print("...Date after  = "+dateA)
        
        yyyymmddARR.append(dateB4)
        yyyymmddARR.append(yyyymmdd)
        yyyymmddARR.append(dateA)

        yyyymmddhh = yyyymmdd+hh
        yyyymmddhhARR.append(yyyymmddhh)

	#=============================================
	# Define Datasets

		#---------------------------------------------
		# Define and load INDEX FILES
		#---------------------------------------------
        print("*****Define INDEX FILES")

		# read file
        if input_path == str(archive_parent)+"collocation/index_files":	# index files are located within 'archive_parent'
          input_path_str = input_path+"/"+str(iyy)+"/"+str(mm)+"/"+str(dd)+"/"
        else:								# index files are somewhere else
          input_path_str = input_path
        dset_path	      = input_path_str
        dset_filename	      = 'index.'+yyyymmddhh 
        del input_path_str

        tidx_file_vertstr = [[] for i in range(tndset)]
        qc_dset_list      = []
        dset_src          = []
        
        for j in range(int(tndset)):
          tfile = dset_path+dset_filename+suffix[j]
          dset_exists = exists(tfile)
          if dset_exists==False:
            print("WARNING: index file "+tfile+" does not exist!")
            tidx_file_vertstr[j] = "NO_MATCHES"
            d_names = str(fill)
            continue
	      	# load dataset
          data_hdl = Dataset(tfile)

	   	# extract data
			# collocation criteria
          ttim_max   = np.asarray( data_hdl.variables['time_max'] )
          tprs_max   = np.asarray( data_hdl.variables['pres_max'] )
          thgt_max   = np.asarray( data_hdl.variables['hgt_max' ] )
          tdst_max   = np.asarray( data_hdl.variables['dist_max'] )
          tim_max[j] = ttim_max
          prs_max[j] = tprs_max
          hgt_max[j] = thgt_max
          dst_max[j] = tdst_max
          print("Collocation Criteria for "+str(tfile)+":")
          print("... Max time = "+str(tim_max[j])+" | distance = "+str(dst_max[j])+" | pressure = "+str(prs_max[j])+" | height = "+str(hgt_max[j]))
          del ttim_max,tprs_max,thgt_max,tdst_max
	  
			# DRIVER dataset attributes
          drv_key = 'drv'
          for attname in data_hdl.variables[drv_key].ncattrs():
            att = data_hdl.variables[drv_key].getncattr(attname)
            if attname == 'short_name':
              driver_name = att
            if attname == 'wind_type':
              aeolus_wind_type = att
            if attname == 'dataset_type':
              tatt = att
              if tatt.find("Reproc") != -1:
                tatt = tatt.strip("Reproc")
              else:
                tatt = "orig"
              aeolus_dset_type_str = tatt
              del tatt
            if attname == 'QC_flag':
              bool_drv_qc = False				# don't apply any QC. Need to extract the FULL dataset.
            if attname == 'QC_list':
              if ((driver_name.find('AMV') != -1) or (driver_name.find('amv') != -1)):
                pct_drv = att
                pct_drv = pct_drv.strip(' %')
            del att
          del tfile,dset_exists

          idx_file_str = suffix[j]
	  
          print("Read index file")
          dset1_name,idx_drv_dset1,idx_dset1,DT_match1,GCD_match1,Vert_match1,Vert_str1 = read_index_file(archive_parent,yyyymmddhh,idx_file_str,dset_path)
          if dset1_name=="NoMatches": continue 		# double check that index file exists
	  
          print("Set collocation criteria arrays")
          dset_name[j]         = dset1_name
          DT_match[j]          = np.append(DT_match[j], DT_match1, axis=0)
          GCD_match[j]         = np.append(GCD_match[j], GCD_match1, axis=0)
          Vert_match[j]        = np.append(Vert_match[j], Vert_match1, axis=0)
          tidx_file_vertstr[j] = Vert_str1
 
          del dset1_name,DT_match1,GCD_match1,Vert_match1,Vert_str1
          
 	  #---------------------------------------------
	  # Define and load DRIVER dataset
	  #---------------------------------------------
          print("*****Define DRIVER dataset")
	  
          idxs = idx_drv_dset1

          if (driver_name.find('Aeolus') != -1):
            tdrv_lat,tdrv_lon,tdrv_prs,tdrv_hgt,tdrv_yr,tdrv_mm,tdrv_dy,tdrv_hr,tdrv_mn,tindexesD,tqc_drv_list,tdrv_src,tdrv_err,tdrv_len,tdrv_spd,tdrv_dir = read_aeolus(archive_parent,yyyymmddhh,dateB4,dateA,aeolus_dset_type_str,aeolus_wind_type,bool_drv_qc,"drv",runtype,tim_max[j],idxs)
            del tindexesD,tqc_drv_list,tdrv_src,tdrv_err,tdrv_len
          elif driver_name=='Aircraft':
            tdrv_lat,tdrv_lon,tdrv_yr,tdrv_mm,tdrv_dy,tdrv_hr,tdrv_mn,tdrv_hgt,tdrv_prs,tindexesD,tqc_drv_list,tdrv_src,tdrv_spd,tdrv_dir = read_aircraft(archive_parent,yyyymmddhh,dateB4,dateA,bool_drv_qc,"drv",runtype,tim_max[j],idxs)
            del tindexesD,tqc_drv_list,tdrv_src
          elif driver_name=='AMV_NCEP':
            tdrv_lat,tdrv_lon,tdrv_yr,tdrv_mm,tdrv_dy,tdrv_hr,tdrv_mn,tdrv_hgt,tdrv_prs,tindexesD,tqc_drv_list,tdrv_src,tdrv_satname,tdrv_wcm,tdrv_ham,tdrv_spd,tdrv_dir = read_amv_ncep(archive_parent,yyyymmddhh,dateB4,dateA,bool_drv_qc,pct_drv,qi_choice,"drv",runtype,tim_max[j],idxs)
            del tindexesD,tqc_drv_list,tdrv_src,tdrv_satname,tdrv_ham
          elif driver_name=='Loon':
            tdrv_lat,tdrv_lon,tdrv_yr,tdrv_mm,tdrv_dy,tdrv_hr,tdrv_mn,tdrv_hgt,tdrv_prs,tindexesD,tqc_drv_list,tdrv_src,tdrv_hgt,tdrv_azm,tdrv_elv,tdrv_spd,tdrv_dir = read_loon(archive_parent,yyyymmddhh,dateB4,dateA,bool_drv_qc,"drv",runtype,tim_max[j],idxs)
            del tindexesD,tqc_drv_list,tdrv_src,tdrv_azm,tdrv_elv
          elif driver_name=='Sonde':
            tdrv_lat,tdrv_lon,tdrv_yr,tdrv_mm,tdrv_dy,tdrv_hr,tdrv_mn,tdrv_hgt,tdrv_prs,tindexesD,tqc_drv_list,tdrv_src,tnsondes,tnlevels,tngroups,tdrv_spd,tdrv_dir = read_sonde(archive_parent,yyyymmddhh,dateB4,dateA,bool_drv_qc,"drv",runtype,tim_max[j],idxs)
            del tindexesD,tqc_drv_list,tdrv_src,tnsondes,tnlevels,tngroups
          elif driver_name=='DAWN':
            tdrv_lat,tdrv_lon,tdrv_yr,tdrv_mm,tdrv_dy,tdrv_hr,tdrv_mn,tdrv_hgt,tdrv_prs,tindexesD,tqc_drv_list,tdrv_src,tdrv_hgt,tdrv_spd,tdrv_dir = read_DAWN(archive_parent,yyyymmddhh,dateB4,dateA,bool_drv_qc,"drv",runtype,tim_max[j],idxs)
            del tindexesD,tqc_drv_list,tdrv_src,tdrv_azm,tdrv_elv
          else:
            print("ERROR: Driver dataset not defined!")
            sys.exit()
	    
          del idxs
          	  
  	  #---------------------------------------------
  	  # Define dependent datasets and append to each other to create total matching array for collocation
  	  #---------------------------------------------
          print("*****Define DEPENDENT dataset(s)")

          dset_key = "dset1"

		# get attributes for each dataset
          for attname in data_hdl.variables[dset_key].ncattrs():
            att = data_hdl.variables[dset_key].getncattr(attname)
            if attname == 'short_name':
              d_names = att
            if attname == 'wind_type':
              aeolus_wind_type = att
            if attname == 'dataset_type':
              tatt = att
              if tatt.find("Reproc") != -1:
                tatt = tatt.strip("Reproc")
              else:
                tatt = "orig"
              aeolus_dset_type_str = tatt
              del tatt
            if attname == 'QC_flag':
              bool_dset_qc = False
            if attname == 'QC_list':
              if ((d_names.find('AMV') != -1) or (d_names.find('amv') != -1)):
                pct = att
                pct = pct.strip(' %')
            del att

          data_hdl.close()

          if tidx_file_vertstr[j] != "NO_MATCHES":			# check if DEPENDENT dataset has any matches to DRIVER. If yes, proceed.

            idx_file_vertstr[j] = tidx_file_vertstr[j]
	    
            idxs = idx_dset1

		# read dependent dataset
            if (d_names.find('Aeolus') != -1):
              t_lat,t_lon,t_prs,t_hgt,t_yr,t_mm,t_dy,t_hr,t_mn,t_indexes,t_qc_dset_list,t_src,t_err,t_len,t_spd,t_dir = read_aeolus(archive_parent,yyyymmddhh,dateB4,dateA,aeolus_dset_type_str,aeolus_wind_type,bool_dset_qc,"dep",runtype,tim_max[j],idxs)
              del t_err,t_len
            elif d_names=='Aircraft':
              t_lat,t_lon,t_yr,t_mm,t_dy,t_hr,t_mn,t_hgt,t_prs,t_indexes,t_qc_dset_list,t_src,t_spd,t_dir = read_aircraft(archive_parent,yyyymmddhh,dateB4,dateA,bool_dset_qc,"dep",runtype,tim_max[j],idxs)
            elif d_names=='AMV_NCEP':
              t_lat,t_lon,t_yr,t_mm,t_dy,t_hr,t_mn,t_hgt,t_prs,t_indexes,t_qc_dset_list,t_src,t_satname,t_wcm,t_ham,t_spd,t_dir = read_amv_ncep(archive_parent,yyyymmddhh,dateB4,dateA,bool_dset_qc,pct,qi_choice,"dep",runtype,tim_max[j],idxs)
              del t_satname,t_ham
#            elif d_names=='Loon':
#              t_lat,t_lon,t_yr,t_mm,t_dy,t_hr,t_mn,t_hgt,t_prs,t_indexes,t_qc_dset_list,t_src,t_azm,t_elv,t_spd,t_dir = read_loon(archive_parent,yyyymmddhh,dateB4,dateA,bool_dset_qc,"dep",runtype,tim_max[j],idxs)
#              del t_azm,t_elv
            elif d_names=='Sonde':
		# RADIOSONDE data are 2D (except for lat, lon, dates/times): [nsondes, nlevels]
		# 	nsondes = number of sondes
		# 	nlevels = number of vertical levels per sonde
		# nsondes, nlevels are used to find collocated indices that pertain to original data file
              t_lat,t_lon,t_yr,t_mm,t_dy,t_hr,t_mn,t_hgt,t_prs,t_indexes,t_qc_dset_list,t_src,nsondes,nlevels,ngroups,t_spd,t_dir = read_sonde(archive_parent,yyyymmddhh,dateB4,dateA,bool_dset_qc,"dep",runtype,tim_max[j],idxs)
              del nsondes,nlevels,ngroups
            elif d_names=='DAWN':
              t_lat,t_lon,t_yr,t_mm,t_dy,t_hr,t_mn,t_hgt,t_prs,t_indexes,t_qc_dset_list,t_src,t_spd,t_dir = read_DAWN(archive_parent,yyyymmddhh,dateB4,dateA,bool_dset_qc,"dep",runtype,tim_max[j],idxs)
            else:
              print("ERROR: Dependent dataset "+str(i)+" ("+d_names+") cannot be read by this program! Please add function to read_input_datasets module and try again.")
              sys.exit()

            del idxs
          del t_yr,t_mm,t_dy,t_hr,t_mn,t_indexes,t_qc_dset_list,t_src
	  
	  #---------------------------------------------
          #---------------------------------------------

          	#---------------------------------------------
	  	# reassign all missing values to NaN
          ttdrv_yr  = to_nan(tdrv_yr )
          ttdrv_mm  = to_nan(tdrv_mm )
          ttdrv_dy  = to_nan(tdrv_dy )
          ttdrv_hr  = to_nan(tdrv_hr )
          ttdrv_mn  = to_nan(tdrv_mn )
          ttdrv_lat = to_nan(tdrv_lat)
          ttdrv_lon = to_nan(tdrv_lon)
          ttdrv_prs = to_nan(tdrv_prs)
          ttdrv_hgt = to_nan(tdrv_hgt)
          ttdrv_spd = to_nan(tdrv_spd)
          ttdrv_dir = to_nan(tdrv_dir)
          if driver_name=="AMV_NCEP":
            ttdrv_wcm = to_nan(tdrv_wcm)
            del tdrv_wcm
          tt_lat    = to_nan(t_lat   )
          tt_lon    = to_nan(t_lon   )
          tt_prs    = to_nan(t_prs   )
          tt_hgt    = to_nan(t_hgt   )
          tt_spd    = to_nan(t_spd   )
          tt_dir    = to_nan(t_dir   )
          if d_names=="AMV_NCEP":
            tt_wcm  = to_nan(t_wcm)
            del t_wcm
          del tdrv_yr,tdrv_mm,tdrv_dy,tdrv_hr,tdrv_mn
          del tdrv_lat,tdrv_lon,tdrv_prs,tdrv_hgt,tdrv_spd,tdrv_dir
          del t_lat,t_lon,t_prs,t_hgt,t_spd,t_dir

          	#---------------------------------------------
		# compute HLOS wind velocity, if applicable
          if (d_names.find('Aeolus') != -1):
            # compute DRIVER hlos if Aeolus is DEPENDENT
            tmp_spd = ttdrv_spd
            tmp_dir = ttdrv_dir
            del ttdrv_spd,ttdrv_dir
            ttdrv_spd  = np.nan * np.ones_like(ttdrv_prs)
            ttdrv_dir  = np.nan * np.ones_like(ttdrv_prs)
            for jj in range(np.size(ttdrv_prs)):
              thlos = compute_hlos(tt_dir[jj],tt_spd[jj],tmp_dir[jj],tmp_spd[jj])
              if thlos>fill:
                ttdrv_spd[jj] = compute_hlos(tt_dir[jj],tt_spd[jj],tmp_dir[jj],tmp_spd[jj])
                ttdrv_dir[jj] = tt_dir[jj]
              del thlos
            del tmp_spd,tmp_dir    

          if (driver_name.find('Aeolus') != -1):
            # compute DEPENDENT hlos if Aeolus is DRIVER
            tmp_spd = tt_spd
            tmp_dir = tt_dir
            del tt_spd,tt_dir
            tt_spd = np.nan * np.ones_like(tt_prs)
            tt_dir = np.nan * np.ones_like(tt_prs)
            for jj in range(np.size(tt_prs)):
              thlos = compute_hlos(ttdrv_dir[jj],ttdrv_spd[jj],tmp_dir[jj],tmp_spd[jj])
              if thlos>fill:
                tt_spd[jj] = compute_hlos(ttdrv_dir[jj],ttdrv_spd[jj],tmp_dir[jj],tmp_spd[jj])
                tt_dir[jj] = ttdrv_dir[jj]
              del thlos
            del tmp_spd,tmp_dir

          	#---------------------------------------------
          	# gross check for wind speed: 
		#	Omit missing values and unrealistic winds (e.g., winds larger than 'spd_max' as set in 'qc_winds' function in quality_controls.py)
                #	'ispd' = indices where inputs meet wind speed QC criteria in 'qc_winds' function
          ispd = qc_winds(ttdrv_spd,tt_spd)
	 
			# matched index arrays
          qidx_drv_dset1 = idx_drv_dset1[ispd]
          qidx_dset1     = idx_dset1[ispd] 
			# variables
          qdrv_yr  = ttdrv_yr [ispd]
          qdrv_mm  = ttdrv_mm [ispd]
          qdrv_dy  = ttdrv_dy [ispd]
          qdrv_hr  = ttdrv_hr [ispd]
          qdrv_mn  = ttdrv_mn [ispd]
          qdrv_lat = ttdrv_lat[ispd]
          qdrv_lon = ttdrv_lon[ispd]
          qdrv_prs = ttdrv_prs[ispd]
          qdrv_hgt = ttdrv_hgt[ispd]
          qdrv_spd = ttdrv_spd[ispd]
          qdrv_dir = ttdrv_dir[ispd]
          if driver_name=="AMV_NCEP":
            qdrv_wcm = ttdrv_wcm[ispd]
            del ttdrv_wcm
          qt_lat   = tt_lat   [ispd]
          qt_lon   = tt_lon   [ispd]
          qt_prs   = tt_prs   [ispd]
          qt_hgt   = tt_hgt   [ispd]
          qt_spd   = tt_spd   [ispd]
          qt_dir   = tt_dir   [ispd]
          if d_names=="AMV_NCEP":
            qt_wcm = tt_wcm[ispd]
            del tt_wcm
          del ttdrv_yr,ttdrv_mm,ttdrv_dy,ttdrv_hr,ttdrv_mn
          del ttdrv_lat,ttdrv_lon,ttdrv_prs,ttdrv_hgt,ttdrv_spd,ttdrv_dir
          del tt_lat,tt_lon,tt_prs,tt_hgt,tt_spd,tt_dir
          del ispd

          	#---------------------------------------------
                # super-ob (if option is chosen in MAINSCRIPT_for__Plotting.bash)

          if avgthin_choice==0:
          	# super-ob matches: average DEPENDENT obs for same DRIVER ob.
            ttdrv_yr  = qdrv_yr 
            ttdrv_mm  = qdrv_mm 
            ttdrv_dy  = qdrv_dy 
            ttdrv_hr  = qdrv_hr 
            ttdrv_mn  = qdrv_mn 
            ttdrv_lat = qdrv_lat
            ttdrv_lon = qdrv_lon
            ttdrv_prs = qdrv_prs
            ttdrv_hgt = qdrv_hgt
            ttdrv_spd = qdrv_spd
            ttdrv_dir = qdrv_dir
            if driver_name=="AMV_NCEP":
              ttdrv_wcm = qdrv_wcm
              del qdrv_wcm
            else:
              ttdrv_wcm = qdrv_lat	# this is just to temporarily fill the array. It will not be used in plotting.
            tt_lat = qt_lat
            tt_lon = qt_lon
            tt_prs = qt_prs
            tt_hgt = qt_hgt
            tt_spd = qt_spd
            tt_dir = qt_dir
            if d_names=="AMV_NCEP":
              tt_wcm = qt_wcm
              del qt_wcm
            else:
              tt_wcm = qt_lat		# this is just to temporarily fill the array. It will not be used in plotting.
            del qdrv_yr,qdrv_mm,qdrv_dy,qdrv_hr,qdrv_mn
            del qdrv_lat,qdrv_lon,qdrv_prs,qdrv_hgt,qdrv_spd,qdrv_dir
            del qt_lat,qt_lon,qt_prs,qt_hgt,qt_spd,qt_dir

            qdrv_yr,qdrv_mm,qdrv_dy,qdrv_hr,qdrv_mn,qdrv_lat,qdrv_lon,qdrv_prs,qdrv_hgt,qdrv_spd,qdrv_dir,qt_lat,qt_lon,qt_prs,qt_hgt,qt_spd,qt_dir,qdrv_wcm,qt_wcm = superob_matches(qidx_drv_dset1,ttdrv_yr,ttdrv_mm,ttdrv_dy,ttdrv_hr,ttdrv_mn,ttdrv_lat,ttdrv_lon,ttdrv_prs,ttdrv_hgt,ttdrv_spd,ttdrv_dir,tt_lat,tt_lon,tt_prs,tt_hgt,tt_spd,tt_dir,ttdrv_wcm,tt_wcm)

            del ttdrv_yr,ttdrv_mm,ttdrv_dy,ttdrv_hr,ttdrv_mn
            del ttdrv_lat,ttdrv_lon,ttdrv_prs,ttdrv_hgt,ttdrv_spd,ttdrv_dir
            del tt_lat,tt_lon,tt_prs,tt_hgt,tt_spd,tt_dir

          del qidx_drv_dset1,qidx_dset1

          	#---------------------------------------------
		# make date array
          ttDyyyymmdd   = []
          ttDyyyymmddhh = []
          for iP in range(np.size(qdrv_yr)):
            jyr = int(qdrv_yr[iP])
            jmm = int(qdrv_mm[iP])
            jdy = int(qdrv_dy[iP])
            jhh = int(qdrv_hr[iP])
            if jmm < 10:
              mms = "0"+str(jmm)
            else:
              mms = str(jmm)
            if jdy < 10:
              dds = "0"+str(jdy)
            else:
              dds = str(jdy)
            if jhh < 10:
              hhs = "0"+str(jhh)
            else:
              hhs = str(jhh)
            tmp   = str(jyr)+str(mms)+str(dds)
            tmphh = tmp+str(jhh)
            ttDyyyymmdd.append(tmp)
            ttDyyyymmddhh.append(tmphh)
            del jyr,jmm,jdy,jhh,mms,dds,hhs,tmp
          del qdrv_yr,qdrv_mm,qdrv_dy,qdrv_hr,qdrv_mn
	  
	  		# convert date list to np.array
          tDyyyymmdd   = np.asarray(ttDyyyymmdd)
          tDyyyymmddhh = np.asarray(ttDyyyymmddhh)
          del ttDyyyymmdd,ttDyyyymmddhh

		# lon/prs/hgt checks
		# ... longitude check: set range to 0=360
			# DRIVER
          st_lon = [qdrv_lon[i]+360.0 if qdrv_lon[i]<0.0 else qdrv_lon[i] for i in range(np.size(qdrv_lon))]
          del qdrv_lon
          qdrv_lon = np.asarray(st_lon)
          del st_lon
			# DEPENDENT
          st_lon = [qt_lon[i]+360.0 if qt_lon[i]<0.0 else qt_lon[i] for i in range(np.size(qt_lon))]
          del qt_lon
          qt_lon = np.asarray(st_lon)
          del st_lon
		# ... pressure check: set pressure to NaN if p<=0
			# DRIVER
          st_prs = [np.nan if qdrv_prs[i]<=0.0 else qdrv_prs[i] for i in range(np.size(qdrv_prs))]
          del qdrv_prs
          qdrv_prs = np.asarray(st_prs)
          del st_prs
			# DEPENDENT
          st_prs = [np.nan if qt_prs[i]<=0.0 else qt_prs[i] for i in range(np.size(qt_prs))]
          del qt_prs
          qt_prs = np.asarray(st_prs)
          del st_prs
		# ... height check: set height to NaN if hgt<0
			# DRIVER
          st_hgt = [np.nan if qdrv_hgt[i]<0.0 else qdrv_hgt[i] for i in range(np.size(qdrv_hgt))]
          del qdrv_hgt
          qdrv_hgt = np.asarray(st_hgt)
          del st_hgt
			# DEPENDENT
          st_hgt = [np.nan if qt_hgt[i]<0.0 else qt_hgt[i] for i in range(np.size(qt_hgt))]
          del qt_hgt
          qt_hgt = np.asarray(st_hgt)
          del st_hgt
	  
          del idx_drv_dset1,idx_dset1
	 
          # append matched obs 
          Dyyyymmdd[j]   = np.append(Dyyyymmdd[j],tDyyyymmdd,axis=0)
          Dyyyymmddhh[j] = np.append(Dyyyymmddhh[j],tDyyyymmddhh,axis=0)
          drv_lat[j]   = np.append(drv_lat[j]  ,qdrv_lat  ,axis=0)
          drv_lon[j]   = np.append(drv_lon[j]  ,qdrv_lon  ,axis=0)
          drv_prs[j]   = np.append(drv_prs[j]  ,qdrv_prs  ,axis=0)
          drv_hgt[j]   = np.append(drv_hgt[j]  ,qdrv_hgt  ,axis=0)
          drv_spd[j]   = np.append(drv_spd[j]  ,qdrv_spd  ,axis=0)
          drv_dir[j]   = np.append(drv_dir[j]  ,qdrv_dir  ,axis=0)
          if driver_name=="AMV_NCEP":
            drv_wcm[j] = np.append(drv_wcm[j]  ,qdrv_wcm  ,axis=0)
            del qdrv_wcm
          dset_lat[j]  = np.append(dset_lat[j] ,qt_lat    ,axis=0)
          dset_lon[j]  = np.append(dset_lon[j] ,qt_lon    ,axis=0)
          dset_prs[j]  = np.append(dset_prs[j] ,qt_prs    ,axis=0)
          dset_hgt[j]  = np.append(dset_hgt[j] ,qt_hgt    ,axis=0)
          dset_spd[j]  = np.append(dset_spd[j] ,qt_spd    ,axis=0)
          dset_dir[j]  = np.append(dset_dir[j] ,qt_dir    ,axis=0)
          if d_names=="AMV_NCEP":
            dset_wcm[j] = np.append(dset_wcm[j],qt_wcm    ,axis=0)
            del qt_wcm
	  
          del d_names
          del tDyyyymmdd

          del qdrv_lat,qdrv_lon,qdrv_prs,qdrv_hgt,qdrv_spd,qdrv_dir
          del qt_lat,qt_lon,qt_prs,qt_hgt,qt_spd,qt_dir
	  	  
          #---------------------------------------------

        del dset_path,dset_filename,tidx_file_vertstr,qc_dset_list,dset_src

        ihh+=1
      idd+=1
    imm+=1
  iyy+=1

#----------------------------------------------------------
# Plotting prep.
#----------------------------------------------------------

# Setup for regional plots
#       NH = Northern Hemisphere = lat range [30,90)
#       TR = Tropics             = lat range [-30,30)
#       SH = Southern Hemisphere = lat range [-90,-30)
regions        = ['NH','TR','SH']
nregions       = np.size(regions)
regions_latmax = [90, 30,-30]
regions_latmin = [30,-30,-90]

# Find number of collocations per dataset/region
Dlatlonlist   = []
DlatlonlistNH = []
DlatlonlistTR = []
DlatlonlistSH = []
nDuniq_list   = [[] for i in range(tndset)]	# global obs count
nDuniq_listNH = [[] for i in range(tndset)]	# NH obs count
nDuniq_listTR = [[] for i in range(tndset)]	# TR obs count
nDuniq_listSH = [[] for i in range(tndset)]	# SH obs count
for j in range(tndset):
	# set string per DEPENDENT dataset to 'NO_MATCHES' if dataset has no matches with DRIVER
  if np.size(idx_file_vertstr[j])==0:
    idx_file_vertstr[j] = "NO_MATCHES"

	# get UNIQUE counts for DRIVER dataset (for legend labels)
  x = drv_lon[j]
  y = drv_lat[j]
  xarr = np.asarray(x)
  yarr = np.asarray(y)
  xNH = xarr[np.where(yarr>regions_latmin[0])]
  xTR = xarr[np.where((yarr>regions_latmin[1])*(yarr<=regions_latmax[1]))]
  xSH = xarr[np.where(yarr<=regions_latmax[2])]
  yNH = yarr[np.where(yarr>regions_latmin[0])]
  yTR = yarr[np.where((yarr>regions_latmin[1])*(yarr<=regions_latmax[1]))]
  ySH = yarr[np.where(yarr<=regions_latmax[2])]
  if idx_file_vertstr[j]=="Height":
    # use height
    z = drv_hgt[j]
  elif idx_file_vertstr[j]=="Pressure":
    # use pressure
    z = drv_prs[j]
  else:
    nDuniq_list[j]   = 0
    nDuniq_listNH[j] = 0
    nDuniq_listTR[j] = 0
    nDuniq_listSH[j] = 0
    continue

  if idx_file_vertstr[j]!="NO_MATCHES":
    zarr = np.asarray(z)
    zNH = zarr[np.where(yarr>regions_latmin[0])]
    zTR = zarr[np.where((yarr>regions_latmin[1])*(yarr<=regions_latmax[1]))]
    zSH = zarr[np.where(yarr<=regions_latmax[2])]
    del zarr

    xstr = np.array(x, dtype=str)
    ystr = np.array(y, dtype=str)
    zstr = np.array(z, dtype=str)
    xstrNH = np.array(xNH, dtype=str)
    xstrTR = np.array(xTR, dtype=str)
    xstrSH = np.array(xSH, dtype=str)
    ystrNH = np.array(yNH, dtype=str)
    ystrTR = np.array(yTR, dtype=str)
    ystrSH = np.array(ySH, dtype=str)
    zstrNH = np.array(zNH, dtype=str)
    zstrTR = np.array(zTR, dtype=str)
    zstrSH = np.array(zSH, dtype=str)
    del zNH,zTR,zSH

    xstr = np.char.add(xstr,",")
    ystr = np.char.add(xstr,ystr)
    ystr = np.char.add(ystr,",")
    xstrNH = np.char.add(xstrNH,",")
    xstrTR = np.char.add(xstrTR,",")
    xstrSH = np.char.add(xstrSH,",")
    ystrNH = np.char.add(xstrNH,ystrNH)
    ystrNH = np.char.add(ystrNH,",")
    ystrTR = np.char.add(xstrTR,ystrTR)
    ystrTR = np.char.add(ystrTR,",")
    ystrSH = np.char.add(xstrSH,ystrSH)
    ystrSH = np.char.add(ystrSH,",")    

    ttDlatlonlist   = np.char.add(ystr,zstr)
    ttDlatlonlistNH = np.char.add(ystrNH,zstrNH)
    ttDlatlonlistTR = np.char.add(ystrTR,zstrTR)
    ttDlatlonlistSH = np.char.add(ystrSH,zstrSH)
    del xstr,ystr,zstr,xstrNH,xstrTR,xstrSH,ystrNH,ystrTR,ystrSH,zstrNH,zstrTR,zstrSH
    tDlatlonlist   = list(set(ttDlatlonlist))	# find unique values
    tDlatlonlistNH = list(set(ttDlatlonlistNH))	# find unique values
    tDlatlonlistTR = list(set(ttDlatlonlistTR))	# find unique values
    tDlatlonlistSH = list(set(ttDlatlonlistSH))	# find unique values
    Dlatlonlist   = np.append(Dlatlonlist, tDlatlonlist, axis=0)
    DlatlonlistNH = np.append(DlatlonlistNH, tDlatlonlistNH, axis=0)
    DlatlonlistTR = np.append(DlatlonlistTR, tDlatlonlistTR, axis=0)
    DlatlonlistSH = np.append(DlatlonlistSH, tDlatlonlistSH, axis=0)
    del ttDlatlonlist,tDlatlonlist,ttDlatlonlistNH,tDlatlonlistNH,ttDlatlonlistTR,tDlatlonlistTR,ttDlatlonlistSH,tDlatlonlistSH
		# count unique values per DEPENDENT dataset
    tmp = np.size(Dlatlonlist)
    nDuniq_list[j] = tmp
    del tmp
    tmp = np.size(DlatlonlistNH)
    nDuniq_listNH[j] = tmp
    del tmp
    tmp = np.size(DlatlonlistTR)
    nDuniq_listTR[j] = tmp
    del tmp
    tmp = np.size(DlatlonlistSH)
    nDuniq_listSH[j] = tmp
    del tmp
  del x,y,z,xarr,yarr

# count unique number of DRIVER obs
	# global
Duniq_listALL = list(set(Dlatlonlist))
nDuniq_listALL = np.size(Duniq_listALL)
del Dlatlonlist,Duniq_listALL
print("number of unique DRIVER obs global = "+str(nDuniq_listALL))
        # NH
Duniq_listNH = list(set(DlatlonlistNH))
nDuniq_listNH = np.size(Duniq_listNH)
del DlatlonlistNH,Duniq_listNH
print("number of unique DRIVER obs NH = "+str(nDuniq_listNH))
        # TR
Duniq_listTR = list(set(DlatlonlistTR))
nDuniq_listTR = np.size(Duniq_listTR)
del DlatlonlistTR,Duniq_listTR
print("number of unique DRIVER obs TR = "+str(nDuniq_listTR))
        # SH
Duniq_listSH = list(set(DlatlonlistSH))
nDuniq_listSH = np.size(Duniq_listSH)
del DlatlonlistSH,Duniq_listSH
print("number of unique DRIVER obs SH = "+str(nDuniq_listSH))

	# substring for plot filenames
if avgthin_choice==-1:
  print("PLOT ALL MATCHES")
  avgthin_str = ".AllMatches"
elif avgthin_choice==0:
  print("SUPER-OB (AVERAGE) MATCHES")
  avgthin_str = ".SuperOb"

	# get unique dates (yyyymmdd)
date_uniq = list(set(yyyymmddARR))
date_uniq.sort()
        # get unique dates (yyyymmddhh)
datehh_uniq = list(set(yyyymmddhhARR))
datehh_uniq.sort()

#******************************************************************************************************************
#******************************************************************************************************************
# PLOTS
#******************************************************************************************************************
#******************************************************************************************************************
print("***** PLOTS *****")

#----------------------------------------------------------
# More plotting prep.
#----------------------------------------------------------

dateIN = dateSTART+"-"+dateEND

amvtypestr = ["IR","WVcloud","WVclear","Visible"]			# AMV type
iamvtypes  = [1,3,5,2]							# value of WCM pertaining to each AMV type
amvcolors  = ["firebrick","deepskyblue","green","magenta"]		# color assigned to each AMV type

dcolors   = ["red","limegreen","blue","gold","cyan"]			# colors assigned to each DEPENDENT dataset
imark     = ["^","v","s","P","o"]					# markers assigned to each DEPENDENT dataset

names  = dependent_names
colors = dcolors[0:tndset]
marks  = imark[0:tndset]
print("dset names = "+str(names))
print("colors per dset = "+str(colors))
print("markers per dset = "+str(marks))

# Check if all idx_file* = NO_MATCHES. If true, this dataset is skipped in plotting.
nidx = np.zeros(tndset,dtype=int)
for j in range(tndset):
  if idx_file_vertstr[j]!="NO_MATCHES":
    nidx[j] = np.size(dset_lat[j])					# count how many obs match DRIVER
print("number of obs matching DRIVER = "+str(nidx))

sizes = nidx

# Vertical levels for plotting
levsP = [70.0,200.0,500.0]
levsZ = [5.0,7.0,10.0]

# Longitude ranges for vertical cross sections
# 	First value is min lon, second value is max lon. 
#	Obs with lons between min and max are averaged and plotted on VERT vs LAT plot
# NOTE: *min and *max must have same shape
vertcrossLONSmin = [0., 90., 180., 270.]	
vertcrossLONSmax = [90., 180., 270., 360.]

#==========================================================
#==========================================================
# Histograms of Collocation Difference
#==========================================================
#==========================================================

alphaval = 1.0		# transparency factor (between 0 and 1, with 1=opaque)

x_name  = driver_name
py_name = ""
zy_name = ""
if driver_name == 'Aeolus':
  x_name += "_"+aeolus_wind_type+"_"+aeolus_dset_type_str

# y_name string for plot filename
y_name = ""
for jj in range(np.size(names)):
  y_name += "_"+names[jj]
  if names[jj] == 'Aeolus':
    y_name += "_"+aeolus_wind_type+"_"+aeolus_dset_type_str

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# ... Global
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
print("HISTOGRAMS - Global")

latD        = []
lonD        = []
a_tmatch    = []
a_pmatch    = []
a_zmatch    = []
a_distmatch = []

a_tname     = []
a_ptname    = []
a_ztname    = []

a_tcolors    = []
a_pcolors    = []
a_zcolors    = []
a_distcolors = []

for i in range(tndset):
  tt=0
  pmatch = []
  zmatch = []
  if nidx[i]>0:
    tmatch    = DT_match[i]
    distmatch = GCD_match[i]
    tlat      = drv_lat[i]
    tlon      = drv_lon[i]
    tVert_str = idx_file_vertstr[i]
    if tVert_str=="Pressure":
      pmatch  = Vert_match[i]
    else:
      zmatch  = Vert_match[i]
    tname     = dset_name[i]

    	# Arrays for histograms
    a_tname.append(tname)
    latD.append(tlat)
    lonD.append(tlon)
    a_tmatch.append(tmatch)
    a_distmatch.append(distmatch)
    a_tcolors.append(dcolors[i])
    a_distcolors.append(dcolors[i])

    if np.size(pmatch)>0:
      a_pmatch.append(pmatch)
      a_pcolors.append(dcolors[i])
      a_ptname.append(tname)
      py_name += "_"+str(tname)
      del pmatch
    if np.size(zmatch)>0:
      zzmatch = zmatch*1000.0		# convert km to m
      a_zmatch.append(zzmatch)
      a_zcolors.append(dcolors[i])
      a_ztname.append(tname)
      zy_name += "_"+str(tname)
      del zmatch,zzmatch
  
    del tmatch,distmatch,tname,tlat,tlon,tVert_str

np.warnings.filterwarnings('ignore', category=np.VisibleDeprecationWarning)

# TIME DIFFERENCE
match_str = "Time"
units = "min"
hist_diffs(nDuniq_listALL, latD, lonD, a_tmatch, x_name, a_tname, alphaval, match_str, a_tcolors, units)
outname = output_path+"HIST.Match_TimeDiff.Global."+str(dateIN)+".x_"+str(x_name)+".y"+str(y_name)+avgthin_str+".png"
plt.savefig(outname)
del outname,a_tmatch,units

# DISTANCE
match_str = "Colloc. Distances"
units = "km"
hist_diffs(nDuniq_listALL, latD, lonD, a_distmatch, x_name, a_tname, alphaval, match_str, a_distcolors, units)
outname = output_path+"HIST.Match_Distance.Global."+str(dateIN)+".x_"+str(x_name)+".y"+str(y_name)+avgthin_str+".png"
plt.savefig(outname)
del outname,a_distmatch,units
		
# VERTICAL
# ... PRESSURE DIFFERENCE
if np.size(a_pmatch)>0:
  match_str = "Pressure"
  units = "hPa"
  hist_diffs(nDuniq_listALL, latD, lonD, a_pmatch, x_name, a_ptname, alphaval, match_str, a_pcolors, units)
  outname = output_path+"HIST.Match_PressureDiff.Global."+str(dateIN)+".x_"+str(x_name)+".y"+str(py_name)+avgthin_str+".png"
  plt.savefig(outname)
  del outname,a_pmatch,units
# ... HEIGHT DIFFERENCE
if np.size(a_zmatch)>0:
  match_str = "Height"
  units = "m"
  hist_diffs(nDuniq_listALL, latD, lonD, a_zmatch, x_name, a_ztname, alphaval, match_str, a_zcolors, units)
  outname = output_path+"HIST.Match_HeightDiff.Global."+str(dateIN)+".x_"+str(x_name)+".y"+str(zy_name)+avgthin_str+".png"
  plt.savefig(outname)
  del outname,a_zmatch,units

del latD,lonD,py_name,zy_name,a_tname,a_ptname,a_ztname,a_tcolors,a_pcolors,a_zcolors,a_distcolors

plt.close("all")
	
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# ... Regional
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
print("HISTOGRAMS - Regional")

py_name = ""
zy_name = ""

NHlatD        = []
NHlonD        = []
NHa_tmatch    = []
NHa_pmatch    = []
NHa_zmatch    = []
NHa_distmatch = []
TRlatD        = []
TRlonD        = []
TRa_tmatch    = []
TRa_pmatch    = []
TRa_zmatch    = []
TRa_distmatch = []
SHlatD        = []
SHlonD        = []
SHa_tmatch    = []
SHa_pmatch    = []
SHa_zmatch    = []
SHa_distmatch = []

a_tname     = []
a_ptname    = []
a_ztname    = []

a_tcolors    = []
a_pcolors    = []
a_zcolors    = []
a_distcolors = []

vertstr = [[] for ir in range(nregions)]

for i in range(tndset):
  tt=0
  pmatch = []
  zmatch = []
  if nidx[i]>0:
    tmatch    = DT_match[i]
    distmatch = GCD_match[i]
    tlat      = drv_lat[i]
    tlon      = drv_lon[i]
    tVert_str = idx_file_vertstr[i]
    if tVert_str=="Pressure":
      pmatch  = Vert_match[i]
    else:
      zmatch  = Vert_match[i]
    tname     = dset_name[i]

    vertstr.append(tVert_str)

	# Split arrays into regions
    idxNH = np.where((tlat>=regions_latmin[0])*(tlat<regions_latmax[0]))
    idxTR = np.where((tlat>=regions_latmin[1])*(tlat<regions_latmax[1]))
    idxSH = np.where((tlat>=regions_latmin[2])*(tlat<regions_latmax[2]))

    	# Arrays for histograms
    NHlatD.append(tlat[idxNH])
    NHlonD.append(tlon[idxNH])
    NHa_tmatch.append(tmatch[idxNH])
    NHa_distmatch.append(distmatch[idxNH])
    TRlatD.append(tlat[idxTR])
    TRlonD.append(tlon[idxTR])
    TRa_tmatch.append(tmatch[idxTR])
    TRa_distmatch.append(distmatch[idxTR])
    SHlatD.append(tlat[idxSH])
    SHlonD.append(tlon[idxSH])
    SHa_tmatch.append(tmatch[idxSH])
    SHa_distmatch.append(distmatch[idxSH])

    a_tname.append(tname)
    a_tcolors.append(dcolors[i])
    a_distcolors.append(dcolors[i])

    if np.size(pmatch)>0:
      NHa_pmatch.append(pmatch[idxNH])
      TRa_pmatch.append(pmatch[idxTR])
      SHa_pmatch.append(pmatch[idxSH])
      a_pcolors.append(dcolors[i])
      a_ptname.append(tname)
      py_name += "_"+str(tname)
      del pmatch
    if np.size(zmatch)>0:
      zzmatch = zmatch*1000.0		# convert km to m
      NHa_zmatch.append(zzmatch[idxNH])
      TRa_zmatch.append(zzmatch[idxTR])
      SHa_zmatch.append(zzmatch[idxSH])
      a_zcolors.append(dcolors[i])
      a_ztname.append(tname)
      zy_name += "_"+str(tname)
      del zmatch,zzmatch
  
    del tmatch,distmatch,tname,tlat,tlon,tVert_str
    del idxNH,idxTR,idxSH

for ir in range(nregions):
  if regions[ir]=='NH':
    tnDuniq = nDuniq_listNH
    latD=NHlatD; lonD=NHlonD; a_tmatch=NHa_tmatch; a_distmatch=NHa_distmatch
    if vertstr[j]=="Pressure": a_pmatch=NHa_pmatch
    if vertstr[j]=="Height": a_zmatch=NHa_zmatch
  elif regions[ir]=='TR':
    tnDuniq = nDuniq_listTR
    latD=TRlatD; lonD=TRlonD; a_tmatch=TRa_tmatch; a_distmatch=TRa_distmatch
    if vertstr[j]=="Pressure": a_pmatch=TRa_pmatch
    if vertstr[j]=="Height": a_zmatch=TRa_zmatch
  elif regions[ir]=='SH':
    tnDuniq = nDuniq_listSH
    latD=SHlatD; lonD=SHlonD; a_tmatch=SHa_tmatch; a_distmatch=SHa_distmatch
    if vertstr[j]=="Pressure": a_pmatch=SHa_pmatch
    if vertstr[j]=="Height": a_zmatch=SHa_zmatch

  # TIME DIFFERENCE
  match_str = str(regions[ir])+" Time"
  units = "min"
  hist_diffs(tnDuniq, latD, lonD, a_tmatch, x_name, a_tname, alphaval, match_str, a_tcolors, units)
  outname = output_path+"HIST.Match_TimeDiff."+regions[ir]+"."+str(dateIN)+".x_"+str(x_name)+".y"+str(y_name)+avgthin_str+".png"
  plt.savefig(outname)
  del outname,a_tmatch,units

  # DISTANCE
  match_str = str(regions[ir])+" Colloc. Distances"
  units = "km"
  hist_diffs(tnDuniq, latD, lonD, a_distmatch, x_name, a_tname, alphaval, match_str, a_distcolors, units)
  outname = output_path+"HIST.Match_Distance."+regions[ir]+"."+str(dateIN)+".x_"+str(x_name)+".y"+str(y_name)+avgthin_str+".png"
  plt.savefig(outname)
  del outname,a_distmatch,units
		
  # VERTICAL
  # ... PRESSURE DIFFERENCE
  if vertstr[j]=="Pressure":
    match_str = str(regions[ir])+" Pressure"
    units = "hPa"
    hist_diffs(tnDuniq, latD, lonD, a_pmatch, x_name, a_ptname, alphaval, match_str, a_pcolors, units)
    outname = output_path+"HIST.Match_PressureDiff."+regions[ir]+"."+str(dateIN)+".x_"+str(x_name)+".y"+str(py_name)+avgthin_str+".png"
    plt.savefig(outname)
    del outname,a_pmatch,units
  # ... HEIGHT DIFFERENCE
  if vertstr[j]=="Height":
    match_str = str(regions[ir])+" Height"
    units = "m"
    hist_diffs(tnDuniq, latD, lonD, a_zmatch, x_name, a_ztname, alphaval, match_str, a_zcolors, units)
    outname = output_path+"HIST.Match_HeightDiff."+regions[ir]+"."+str(dateIN)+".x_"+str(x_name)+".y"+str(zy_name)+avgthin_str+".png"
    plt.savefig(outname)
    del outname,a_zmatch,units
  del tnDuniq

del y_name,zy_name,py_name, a_tname,a_ptname,a_ztname,a_tcolors,a_pcolors,a_zcolors,a_distcolors
del NHlatD,NHlonD,NHa_tmatch,NHa_pmatch,NHa_zmatch,NHa_distmatch
del TRlatD,TRlonD,TRa_tmatch,TRa_pmatch,TRa_zmatch,TRa_distmatch
del SHlatD,SHlonD,SHa_tmatch,SHa_pmatch,SHa_zmatch,SHa_distmatch

plt.close("all")

#==========================================================
#==========================================================
# Scatterplots and Timeseries of Matched Obs
#==========================================================
#==========================================================
print("SCATTERPLOTS and TIMESERIES")
	
marksize = 50
alphaval = 0.25

x_name = driver_name
y_name_hlos = ""
y_name_wspd = ""
y_name_pres = ""
y_name_hgt  = ""
if driver_name == 'Aeolus':
  x_name += "_"+aeolus_wind_type+"_"+aeolus_dset_type_str

yyyymmddD 	= []
yyyymmddhhD 	= []
latD  		= []
lonD  		= []
PyyyymmddD 	= []
PyyyymmddhhD 	= []
PlatD  		= []
PlonD  		= []
Pxvert 		= []
HyyyymmddD 	= []
HyyyymmddhhD 	= []
HlatD  		= []
HlonD  		= []
Hxvert 		= []

hlosD = []
hlos  = []
wspdD = []
wspd  = []
presD = []
pres  = []
hgtD  = []
hgt   = []

a_tname_hlos = []
a_tname_wspd = []
a_tname_pres = []
a_tname_hgt  = []
a_color_hlos = []
a_color_wspd = []
a_color_pres = []
a_color_hgt  = []
a_mark_hlos  = []
a_mark_wspd  = []
a_mark_pres  = []
a_mark_hgt   = []

Vert_str = []

zunits = []
zvar   = []
zstr   = []
txvert = []

nDpres = 0
nDhgts = 0

i=0
for j in range(tndset):
  tt=0
  if nidx[j]>0:
    tt=1
    tDyyyymmdd = Dyyyymmdd[j]
    tDyyyymmddhh = Dyyyymmddhh[j]
    tname = dset_name[j]
    txspd = drv_spd[j]
    txlat = drv_lat[j]
    txlon = drv_lon[j]
    tyspd = dset_spd[j]
    tVert_str = idx_file_vertstr[j]
    if tVert_str=="Pressure":
      txp = drv_prs[j]
      typ = dset_prs[j]
      ttxvert = txp
    else:
      if np.isnan(drv_hgt[j]).all():
        txz = prs_to_hgt(drv_prs[j])
      else:
        txz = drv_hgt[j]
      ttxvert = txz
      tyz = dset_hgt[j]

  if tt==1:
    if tname == 'Aeolus':
      tname += "_"+aeolus_wind_type+"_"+aeolus_dset_type_str

	# lists
    yyyymmddD.append(tDyyyymmdd)
    yyyymmddhhD.append(tDyyyymmddhh)
    latD.append(txlat)
    lonD.append(txlon)
    zstr.append(tVert_str)
    
    txvert.append(ttxvert)

	# append wind speeds
    if driver_name.find('Aeolus') != -1 or tname.find('Aeolus') != -1:	# find() != -1 --> string contains substring!
    	# HLOS wind
      hlosD.append(txspd)
      hlos.append(tyspd)
      a_tname_hlos.append(tname)
      a_color_hlos.append(dcolors[i])
      a_mark_hlos.append(imark[i])
        # y_name string for plot filename
      y_name_hlos += "_"+tname
    if driver_name.find('Aeolus') == -1 and tname.find('Aeolus') == -1:	# find() == -1 --> string does NOT contain substring!
  	# non-HLOS wind
      wspdD.append(txspd)
      wspd.append(tyspd)
      a_tname_wspd.append(tname)
      a_color_wspd.append(dcolors[i])
      a_mark_wspd.append(imark[i])
  	# y_name string for plot filename
      y_name_wspd += "_"+tname
	
	# append pressures or heights
    if tVert_str=="Pressure":
      zvar.append(txp)
      zunits.append("hPa")
      presD.append(txp)
      pres.append(typ)
      PyyyymmddD.append(tDyyyymmdd)
      PyyyymmddhhD.append(tDyyyymmddhh)
      PlatD.append(txlat)
      PlonD.append(txlon)
      Pxvert.append(ttxvert)
      a_tname_pres.append(tname)
      a_color_pres.append(dcolors[i])
      a_mark_pres.append(imark[i])
      nDpres += int(nDuniq_list[i])
  	# y_name string for plot filename
      y_name_pres += "_"+tname
      del txp,typ		
    if tVert_str!="Pressure":
      zvar.append(txz)
      zunits.append("km")
      hgtD.append(txz)
      hgt.append(tyz)
      HyyyymmddD.append(tDyyyymmdd)
      HyyyymmddhhD.append(tDyyyymmddhh)
      HlatD.append(txlat)
      HlonD.append(txlon)
      Hxvert.append(ttxvert)
      a_tname_hgt.append(tname)
      a_color_hgt.append(dcolors[i])
      a_mark_hgt.append(imark[i])
      nDhgts += int(nDuniq_list[i])
  	# y_name string for plot filename
      y_name_hgt += "_"+tname
      del txz,tyz		

    del txspd,tyspd,tname,ttxvert
  i+=1

shape_hlos = np.shape(hlosD)
shape_wspd = np.shape(wspdD)
shape_pres = np.shape(presD)
shape_hgt  = np.shape(hgtD )

adate	 = np.asarray(yyyymmddD   )
alatD    = np.asarray(latD        )
alonD    = np.asarray(lonD        )
azstr    = np.asarray(zstr)

aa_tname_hlos = np.asarray(a_tname_hlos)
amcolors_hlos = np.asarray(a_color_hlos)
ammarks_hlos  = np.asarray(a_mark_hlos )
aa_tname_wspd = np.asarray(a_tname_wspd)
amcolors_wspd = np.asarray(a_color_wspd)
ammarks_wspd  = np.asarray(a_mark_wspd )

axvert = np.asarray(txvert)
del txvert

lonmin = 0.
lonmax = 360.

#```````````````````````````````````````````
# HLOS WIND
if shape_hlos[0]>0:

  match_str = "HLOS"
  units = "m/s"
 
  aa_tname = aa_tname_hlos
  amcolors = amcolors_hlos
  ammarks  = ammarks_hlos

  aDs  = np.asarray(hlosD)
  aass = np.asarray(hlos )

	#`````````````````````````` 
  	# TIME SERIES of differences
		# Global
  regionstr = "Global"
  latmin = -90.
  latmax = 90.
  time_series(nDuniq_listALL,date_uniq,adate,alatD,alonD,axvert,aDs,aass,x_name,aa_tname,match_str,amcolors,units,regionstr,latmin,latmax)
  outname = output_path+"TIME_SERIES.HLOS_Velocity_Diff.DEP-DRV."+regionstr+"."+str(dateIN)+".x_"+str(x_name)+".y"+str(y_name_hlos)+avgthin_str+".png"
  plt.savefig(outname)
  del outname
  		# Regional
  for ir in range(nregions):
    regionstr = str(regions[ir])
    if regionstr=="NH": tnDuniq = nDuniq_listNH
    if regionstr=="TR": tnDuniq = nDuniq_listTR
    if regionstr=="SH": tnDuniq = nDuniq_listSH
    latmin    = regions_latmin[ir]
    latmax    = regions_latmax[ir]
    time_series(tnDuniq,date_uniq,adate,alatD,alonD,axvert,aDs,aass,x_name,aa_tname,match_str,amcolors,units,regionstr,latmin,latmax)
    outname = output_path+"TIME_SERIES.HLOS_Velocity_Diff.DEP-DRV."+regionstr+"."+str(dateIN)+".x_"+str(x_name)+".y"+str(y_name_hlos)+avgthin_str+".png"
    plt.savefig(outname)
    del outname
    del tnDuniq
  
	#`````````````````````````` 
  	# DIFF/SD Vertical Profiles

  zprs = []
  zhgt = []
  for iz in zstr:
    if iz=='Pressure': zprs = 1
    if iz=='Height': zhgt = 1
	
		# Global
  regionstr = "Global"
  latmin = -90.
  latmax = 90.
  if np.size(zprs)>0 or np.size(zhgt)>0:
 		# ... Pressure
    z_distr(fpres,fhgts,alatD,alonD,aDs,aass,zvar,x_name,aa_tname,match_str,azstr,amcolors,units,zunits,regionstr,latmin,latmax)
    outname = output_path+"VERT_DISTR.HLOS_Velocity_Diff.DEP-DRV."+regionstr+"."+str(dateIN)+".x_"+str(x_name)+".y_"+str(y_name_hlos)+avgthin_str+".png"
    plt.savefig(outname)
    del outname

  		# Regional
  for ir in range(nregions):
    regionstr = str(regions[ir])
    latmin    = regions_latmin[ir]
    latmax    = regions_latmax[ir]
    if np.size(zprs)>0 or np.size(zhgt)>0:
 		# ... Pressure
      z_distr(fpres,fhgts,alatD,alonD,aDs,aass,zvar,x_name,aa_tname,match_str,azstr,amcolors,units,zunits,regionstr,latmin,latmax)
      outname = output_path+"VERT_DISTR.HLOS_Velocity_Diff.DEP-DRV."+str(regions[ir])+"."+str(dateIN)+".x_"+str(x_name)+".y_"+str(y_name_hlos)+avgthin_str+".png"
      plt.savefig(outname)
      del outname

  del zprs,zhgt
 
  	#`````````````````````````` 
  	# DIFF/SD vs Wind Speed
	
  		# create bins of DRIVER wind speed
  binstride = 10
  bins = np.arange(-100,100,binstride)
	
		# Global
  regionstr = "Global"
  latmin = -90.
  latmax = 90.
  stats_windbins(bins,alatD,aDs,aass,x_name,aa_tname,match_str,amcolors,units,units,regionstr,latmin,latmax)
  outname = output_path+"STATSvsWSPD.HLOS_Velocity_Diff.DEP-DRV."+regionstr+"."+str(dateIN)+".x_"+str(x_name)+".y"+str(y_name_hlos)+avgthin_str+".png"
  plt.savefig(outname)
  del outname
 
  		# Regional
  for ir in range(nregions):
    regionstr = str(regions[ir])
    latmin    = regions_latmin[ir]
    latmax    = regions_latmax[ir]
    stats_windbins(bins,alatD,aDs,aass,x_name,aa_tname,match_str,amcolors,units,units,regionstr,latmin,latmax)
    outname = output_path+"STATSvsWSPD.HLOS_Velocity_Diff.DEP-DRV."+str(regions[ir])+"."+str(dateIN)+".x_"+str(x_name)+".y"+str(y_name_hlos)+avgthin_str+".png"
    plt.savefig(outname)
    del outname

  	#`````````````````````````` 
 	# SCATTERPLOT
		# sort datasets by size. Plot largest first ... to smallest last.
  sizesH = np.nan * np.ones(shape_hlos[0],dtype=int)
  for i in range(shape_hlos[0]):
    sizesH[i] = np.size(aDs[i])

  sort_size = np.sort(sizesH)[::-1]		# [::-1] = sort in descending order. For ascending, comment out [::-1]
  idx_sort = []
  for i in range(shape_hlos[0]):		# loop to go thru sort_size
    for j in range(shape_hlos[0]):		# loop to go thru size(aDs)
      if np.size(aDs[j])==sort_size[i]:
        idx_sort.append(int(j))

  Ddate   = adate[idx_sort]
  slatD   = alatD[idx_sort]
  slonD   = alonD[idx_sort]
  svert   = axvert[idx_sort]
  Ds      = aDs[idx_sort]
  ss      = aass[idx_sort]
  a_tname = aa_tname[idx_sort]
  mcolors = amcolors[idx_sort]
  mmarks  = ammarks[idx_sort]
  del sort_size,idx_sort,sizesH
  del hlosD,hlos,a_tname_hlos,a_color_hlos,a_mark_hlos
  
		# Global
  regionstr = "Global"
  latmin = -90.
  latmax = 90.
  scatter_matches(nDuniq_listALL,slatD,slonD,svert,Ds,ss,x_name,a_tname,marksize,alphaval,match_str,mcolors,units,mmarks,regionstr,latmin,latmax)
  outname = output_path+"SCATTER.HLOS_Velocity."+regionstr+"."+str(dateIN)+".x_"+str(x_name)+".y"+str(y_name_hlos)+avgthin_str+".png"
  plt.savefig(outname)	
  del outname

		# Regional
  for ir in range(nregions):
    regionstr = str(regions[ir])
    if regionstr=="NH": tnDuniq = nDuniq_listNH
    if regionstr=="TR": tnDuniq = nDuniq_listTR
    if regionstr=="SH": tnDuniq = nDuniq_listSH
    latmin    = regions_latmin[ir]
    latmax    = regions_latmax[ir]
    scatter_matches(tnDuniq,slatD,slonD,svert,Ds,ss,x_name,a_tname,marksize,alphaval,match_str,mcolors,units,mmarks,regionstr,latmin,latmax)
    outname = output_path+"SCATTER.HLOS_Velocity."+regionstr+"."+str(dateIN)+".x_"+str(x_name)+".y"+str(y_name_hlos)+avgthin_str+".png"
    plt.savefig(outname)	
    del outname
    del tnDuniq
	
  del aDs,aass	
  del units,Ddate,slatD,slonD,svert,Ds,ss,a_tname,mcolors,mmarks

  plt.close("all")
  
#```````````````````````````````````````````
# WIND SPEED (not HLOS)
if shape_wspd[0]>0:

  match_str = "Wind Speed"
  units = "m/s"

  aa_tname = aa_tname_wspd
  amcolors = amcolors_wspd
  ammarks  = ammarks_wspd

  aDs = np.asarray(wspdD)
  aass = np.asarray(wspd )

	#`````````````````````````` 
	# TIME SERIES
		# Global
  regionstr = "Global"
  latmin = -90.
  latmax = 90.
  time_series(nDuniq_listALL,date_uniq,adate,alatD,alonD,axvert,aDs,aass,x_name,aa_tname,match_str,amcolors,units,regionstr,latmin,latmax)
  outname = output_path+"TIME_SERIES.WindSpeed_Diff.DEP-DRV."+regionstr+"."+str(dateIN)+".x_"+str(x_name)+".y"+str(y_name_wspd)+avgthin_str+".png"
  plt.savefig(outname)
  del outname
	  	# Regional
  for ir in range(nregions):
    regionstr = str(regions[ir])
    if regionstr=="NH": tnDuniq = nDuniq_listNH
    if regionstr=="TR": tnDuniq = nDuniq_listTR
    if regionstr=="SH": tnDuniq = nDuniq_listSH
    latmin    = regions_latmin[ir]
    latmax    = regions_latmax[ir]
    time_series(tnDuniq,date_uniq,adate,alatD,alonD,axvert,aDs,aass,x_name,aa_tname,match_str,amcolors,units,regionstr,latmin,latmax)
    outname = output_path+"TIME_SERIES.WindSpeed_Diff.DEP-DRV."+regionstr+"."+str(dateIN)+".x_"+str(x_name)+".y"+str(y_name_wspd)+avgthin_str+".png"
    plt.savefig(outname)
    del outname
    del tnDuniq
 
	#`````````````````````````` 
  	# DIFF/SD Vertical Profiles

  zprs = []
  zhgt = []
  for iz in zstr:
    if iz=='Pressure': zprs = 1
    if iz=='Height': zhgt = 1
	
		# Global
  regionstr = "Global"
  latmin = -90.
  latmax = 90.
  if np.size(zprs)>0 or np.size(zhgt)>0:
 		# ... Pressure
    z_distr(fpres,fhgts,alatD,alonD,aDs,aass,zvar,x_name,aa_tname,match_str,azstr,amcolors,units,zunits,regionstr,latmin,latmax)
    outname = output_path+"VERT_DISTR.WindSpeed_Diff.DEP-DRV."+regionstr+"."+str(dateIN)+".x_"+str(x_name)+".y_"+str(y_name_wspd)+avgthin_str+".png"
    plt.savefig(outname)
    del outname
  
  		# Regional
  for ir in range(nregions):
    regionstr = str(regions[ir])
    latmin    = regions_latmin[ir]
    latmax    = regions_latmax[ir]
    if np.size(zprs)>0 or np.size(zhgt)>0:
    		# ... Pressure
      z_distr(fpres,fhgts,alatD,alonD,aDs,aass,zvar,x_name,aa_tname,match_str,azstr,amcolors,units,zunits,regionstr,latmin,latmax)
      outname = output_path+"VERT_DISTR.WindSpeed_Diff.DEP-DRV."+regionstr+"."+str(dateIN)+".x_"+str(x_name)+".y_"+str(y_name_wspd)+avgthin_str+".png"
      plt.savefig(outname)
      del outname
 
  del zprs,zhgt
 
  plt.close("all")

	#`````````````````````````` 
  	# DIFF/SD vs Wind Speed

		# create bins of DRIVER wind speed
  binstride = 10
  bins = np.arange(1,100,binstride)
	
		# Global
  regionstr = "Global"
  latmin = -90.
  latmax = 90.
  
  stats_windbins(bins,alatD,aDs,aass,x_name,aa_tname,match_str,amcolors,units,units,regionstr,latmin,latmax)
  outname = output_path+"STATSvsWSPD.WindSpeed_Diff.DEP-DRV."+regionstr+"."+str(dateIN)+".x_"+str(x_name)+".y"+str(y_name_wspd)+avgthin_str+".png"
  plt.savefig(outname)
  del outname

  		# Regional
  for ir in range(nregions):
    regionstr = str(regions[ir])
    latmin    = regions_latmin[ir]
    latmax    = regions_latmax[ir]
    
    stats_windbins(bins,alatD,aDs,aass,x_name,aa_tname,match_str,amcolors,units,units,regionstr,latmin,latmax)
    outname = output_path+"STATSvsWSPD.WindSpeed_Diff.DEP-DRV."+str(regions[ir])+"."+str(dateIN)+".x_"+str(x_name)+".y"+str(y_name_wspd)+avgthin_str+".png"
    plt.savefig(outname)
    del outname

  plt.close("all")

	#`````````````````````````` 
	# SCATTERPLOT
		# sort datasets by size. Plot largest first ... to smallest last.
  sizesW = np.nan * np.ones(shape_wspd[0],dtype=int)
  for i in range(shape_wspd[0]):
    sizesW[i] = np.size(aDs[i])

  sort_size = np.sort(sizesW)[::-1]		# [::-1] = sort in descending order. For ascending, comment out [::-1]
  idx_sort = []
  for i in range(shape_wspd[0]):		# loop to go thru sort_size
    for j in range(shape_wspd[0]):		# loop to go thru size(aDs)
      if np.size(aDs[j])==sort_size[i]:
        idx_sort.append(int(j))

  Ddate   = adate[idx_sort]
  slatD   = alatD[idx_sort]
  slonD   = alonD[idx_sort]
  svert   = axvert[idx_sort]
  Ds      = aDs[idx_sort]
  ss      = aass[idx_sort]
  a_tname = aa_tname[idx_sort]
  mcolors = amcolors[idx_sort]
  mmarks  = ammarks[idx_sort]
  del sort_size,idx_sort,sizesW
  del wspdD,wspd,a_tname_wspd,a_color_wspd,a_mark_wspd

  		# Global
  regionstr = "Global"
  latmin = -90.
  latmax = 90.
  scatter_matches(nDuniq_listALL,slatD,slonD,svert,Ds,ss,x_name,a_tname,marksize,alphaval,match_str,mcolors,units,mmarks,regionstr,latmin,latmax)
  outname = output_path+"SCATTER.WindSpeed."+regionstr+"."+str(dateIN)+".x_"+str(x_name)+".y"+str(y_name_wspd)+avgthin_str+".png"
  plt.savefig(outname)	
  del outname
  		# Regional
  for ir in range(nregions):
    regionstr = str(regions[ir])
    if regionstr=="NH": tnDuniq = nDuniq_listNH
    if regionstr=="TR": tnDuniq = nDuniq_listTR
    if regionstr=="SH": tnDuniq = nDuniq_listSH
    latmin    = regions_latmin[ir]
    latmax    = regions_latmax[ir]
    scatter_matches(tnDuniq,slatD,slonD,svert,Ds,ss,x_name,a_tname,marksize,alphaval,match_str,mcolors,units,mmarks,regionstr,latmin,latmax)
    outname = output_path+"SCATTER.WindSpeed."+regionstr+"."+str(dateIN)+".x_"+str(x_name)+".y"+str(y_name_wspd)+avgthin_str+".png"
    plt.savefig(outname)	
    del outname
    del tnDuniq 
 
  del aDs,aass
  del units,Ddate,slatD,slonD,Ds,ss,a_tname,mcolors,mmarks
  
  plt.close("all")

#```````````````````````````````````````````
# PRESSURE
if shape_pres[0]>0:

  match_str = "Pressure"
  units = "hPa"

  aDs = np.asarray(presD)
  aass = np.asarray(pres )
  aPdate = np.asarray(PyyyymmddD)
  aPlatD = np.asarray(PlatD)
  aPlonD = np.asarray(PlonD)
  aPxvert = np.asarray(Pxvert)
  Paa_tname = np.asarray(a_tname_pres)
  Pamcolors = np.asarray(a_color_pres)
  Pammarks  = np.asarray(a_mark_pres)

	#`````````````````````````` 
	# TIME SERIES
		# Global
  regionstr = "Global"
  latmin = -90.
  latmax = 90.
  time_series(nDuniq_listALL,date_uniq,aPdate,aPlatD,aPlonD,aPxvert,aDs,aass,x_name,Paa_tname,match_str,Pamcolors,units,regionstr,latmin,latmax)
  outname = output_path+"TIME_SERIES.Pressure_Diff.DEP-DRV."+regionstr+"."+str(dateIN)+".x_"+str(x_name)+".y"+str(y_name_pres)+avgthin_str+".png"
  plt.savefig(outname)
  del outname
  		# Regional
  for ir in range(nregions):
    regionstr = str(regions[ir])
    if regionstr=="NH": tnDuniq = nDuniq_listNH
    if regionstr=="TR": tnDuniq = nDuniq_listTR
    if regionstr=="SH": tnDuniq = nDuniq_listSH
    latmin    = regions_latmin[ir]
    latmax    = regions_latmax[ir]
    time_series(tnDuniq,date_uniq,aPdate,aPlatD,aPlonD,aPxvert,aDs,aass,x_name,Paa_tname,match_str,Pamcolors,units,regionstr,latmin,latmax)
    outname = output_path+"TIME_SERIES.Pressure_Diff.DEP-DRV."+regionstr+"."+str(dateIN)+".x_"+str(x_name)+".y"+str(y_name_pres)+avgthin_str+".png"
    plt.savefig(outname)
    del outname
    del tnDuniq 
 
  plt.close("all")

	#`````````````````````````` 
	# SCATTERPLOT
		# sort datasets by size. Plot largest first ... to smallest last.
  sizesP = np.nan * np.ones(shape_pres[0],dtype=int)
  for i in range(shape_pres[0]):
    sizesP[i] = np.size(aDs[i])

  sort_size = np.sort(sizesP)[::-1]	       # [::-1] = sort in descending order. For ascending, comment out [::-1]
  idx_sort = []
  for i in range(shape_pres[0]):	       # loop to go thru sort_size
    for j in range(shape_pres[0]):	       # loop to go thru size(aDs)
      if np.size(aDs[j])==sort_size[i]:
        idx_sort.append(int(j))

  Ddate   = aPdate[idx_sort]
  slatD   = aPlatD[idx_sort]
  slonD   = aPlonD[idx_sort]
  svert   = aPxvert[idx_sort]
  Ds	  = aDs[idx_sort]
  ss	  = aass[idx_sort]
  a_tname = Paa_tname[idx_sort]
  mcolors = Pamcolors[idx_sort]
  mmarks  = Pammarks[idx_sort]
  del sort_size,idx_sort,sizesP
  del presD,pres,a_tname_pres,a_color_pres,a_mark_pres

		# Global
  regionstr = "Global"
  latmin = -90.
  latmax = 90.
  scatter_matches(nDpres,slatD,slonD,svert,Ds,ss,x_name,a_tname,marksize,alphaval,match_str,mcolors,units,mmarks,regionstr,latmin,latmax)
  outname = output_path+"SCATTER.Pressure."+regionstr+"."+str(dateIN)+".x_"+str(x_name)+".y"+str(y_name_pres)+avgthin_str+".png"
  plt.savefig(outname)	
  del outname
  		# Regional
  for ir in range(nregions):
    regionstr = str(regions[ir])
    latmin    = regions_latmin[ir]
    latmax    = regions_latmax[ir]
    scatter_matches(nDpres,slatD,slonD,svert,Ds,ss,x_name,a_tname,marksize,alphaval,match_str,mcolors,units,mmarks,regionstr,latmin,latmax)
    outname = output_path+"SCATTER.Pressure."+regionstr+"."+str(dateIN)+".x_"+str(x_name)+".y"+str(y_name_pres)+avgthin_str+".png"
    plt.savefig(outname)	
    del outname
  	
  del aDs,aass
  del aPdate,aPlatD,aPlonD,Paa_tname,Pamcolors,Pammarks
  del units,Ddate,slatD,slonD,Ds,ss,a_tname,mcolors,mmarks
  
  plt.close("all")

#```````````````````````````````````````````
# HEIGHT
if shape_hgt[0]>0:

  match_str = "Height"
  units = "km"
  
  aDs = np.asarray(hgtD)
  aass = np.asarray(hgt )
  aHdate = np.asarray(HyyyymmddD)
  aHlatD = np.asarray(HlatD)
  aHlonD = np.asarray(HlonD)
  aHxvert = np.asarray(Hxvert)
  Haa_tname = np.asarray(a_tname_hgt)
  Hamcolors = np.asarray(a_color_hgt)
  Hammarks  = np.asarray(a_mark_hgt)

	#`````````````````````````` 
	# TIME SERIES
		# Global
  regionstr = "Global"
  latmin = -90.
  latmax = 90.
  time_series(nDuniq_listALL,date_uniq,aHdate,aHlatD,aHlonD,aHxvert,aDs,aass,x_name,Haa_tname,match_str,Hamcolors,units,regionstr,latmin,latmax)
  outname = output_path+"TIME_SERIES.Height_Diff.DEP-DRV."+regionstr+"."+str(dateIN)+".x_"+str(x_name)+".y"+str(y_name_hgt)+avgthin_str+".png"
  plt.savefig(outname)
  del outname
  		# Regional
  for ir in range(nregions):
    regionstr = str(regions[ir])
    if regionstr=="NH": tnDuniq = nDuniq_listNH
    if regionstr=="TR": tnDuniq = nDuniq_listTR
    if regionstr=="SH": tnDuniq = nDuniq_listSH
    latmin    = regions_latmin[ir]
    latmax    = regions_latmax[ir]
    time_series(tnDuniq,date_uniq,aHdate,aHlatD,aHlonD,aHxvert,aDs,aass,x_name,Haa_tname,match_str,Hamcolors,units,regionstr,latmin,latmax)
    outname = output_path+"TIME_SERIES.Height_Diff.DEP-DRV."+regionstr+"."+str(dateIN)+".x_"+str(x_name)+".y"+str(y_name_hgt)+avgthin_str+".png"
    plt.savefig(outname)
    del outname
    del tnDuniq
  
  plt.close("all")

	#`````````````````````````` 
  	# SCATTERPLOT
		# sort datasets by size. Plot largest first ... to smallest last.
  sizesH = np.nan * np.ones(shape_hgt[0],dtype=int)
  for i in range(shape_hgt[0]):
    sizesH[i] = np.size(aDs[i])

  sort_size = np.sort(sizesH)[::-1]	       # [::-1] = sort in descending order. For ascending, comment out [::-1]
  idx_sort = []
  for i in range(shape_hgt[0]): 	       # loop to go thru sort_size
    for j in range(shape_hgt[0]):	       # loop to go thru size(aDs)
      if np.size(aDs[j])==sort_size[i]:
        idx_sort.append(int(j))

  Ddate   = aHdate[idx_sort]
  slatD   = aHlatD[idx_sort]
  slonD   = aHlonD[idx_sort]
  svert   = aHxvert[idx_sort]
  Ds	  = aDs[idx_sort]
  ss	  = aass[idx_sort]
  a_tname = Haa_tname[idx_sort]
  mcolors = Hamcolors[idx_sort]
  mmarks  = Hammarks[idx_sort]
  del sort_size,idx_sort,sizesH
  del hgtD,hgt,a_tname_hgt,a_color_hgt,a_mark_hgt

  		# Global
  regionstr = "Global"
  latmin = -90.
  latmax = 90.
  scatter_matches(nDhgts,latD,lonD,svert,Ds,ss,x_name,a_tname,marksize,alphaval,match_str,mcolors,units,mmarks,regionstr,latmin,latmax)
  outname = output_path+"SCATTER.Height."+regionstr+"."+str(dateIN)+".x_"+str(x_name)+".y"+str(y_name_hgt)+avgthin_str+".png"
  plt.savefig(outname)	
  del outname
  		# Regional
  for ir in range(nregions):
    regionstr = str(regions[ir])
    latmin    = regions_latmin[ir]
    latmax    = regions_latmax[ir]
    scatter_matches(nDhgts,latD,lonD,svert,Ds,ss,x_name,a_tname,marksize,alphaval,match_str,mcolors,units,mmarks,regionstr,latmin,latmax)
    outname = output_path+"SCATTER.Height."+regionstr+"."+str(dateIN)+".x_"+str(x_name)+".y"+str(y_name_hgt)+avgthin_str+".png"
    plt.savefig(outname)	
    del outname
  	
  del aDs,aass
  del aHdate,aHlatD,aHlonD,Haa_tname,Hamcolors,Hammarks
  del units,Ddate,slatD,slonD,Ds,ss,a_tname,mcolors,mmarks

  plt.close("all")

del yyyymmddD,yyyymmddhhD,latD,lonD
del PyyyymmddD,PyyyymmddhhD,PlatD,PlonD,Pxvert
del HyyyymmddD,HyyyymmddhhD,HlatD,HlonD,Hxvert
del Vert_str,zunits,zvar,zstr#,txvert

#==========================================================
#==========================================================
# Maps: Locations of matched obs
#==========================================================
#==========================================================
print("MAPS")

alphaval = 0.25		# transparency factor (between 0 and 1, with 1=opaque)

		# select region to plot (CE map)
		#	0 = global
region_flag = 0

x_name = driver_name
y_name = ""	
if driver_name == 'Aeolus':
  x_name += "_"+aeolus_wind_type+"_"+aeolus_dset_type_str

	# y_name string for plot filename
for jj in range(np.size(names)):
  y_name += "_"+names[jj]
  if names[jj] == 'Aeolus':
    y_name += "_"+aeolus_wind_type+"_"+aeolus_dset_type_str

jDs = []
jss = []
jDx = []
jxx = []
jDy = []
jyy = []
ja_tname = []
jmcolors = []
jmmarks  = []
txvert   = []

for j in range(tndset):
  if nidx[j]>0:			# check if DEPENDENT dataset has any matches to DRIVER. If yes, proceed.
    Dts   = drv_spd[j]
    ts    = dset_spd[j]
    Dtx   = drv_lon[j]
    Dty   = drv_lat[j]
    tx    = dset_lon[j]
    ty    = dset_lat[j]
    tname = dset_name[j]
    
    tVert_str = idx_file_vertstr[j]
    if tVert_str=="Pressure":
      ttxvert = drv_prs[j]
    else:
      if np.isnan(drv_hgt[j]).all():
        ttxvert = prs_to_hgt(drv_prs[j])
      else:
        ttxvert = drv_hgt[j]
      
      # append data
    jDs.append(Dts)
    jss.append(ts)
    jDx.append(Dtx)
    jxx.append(tx)
    jDy.append(Dty)
    jyy.append(ty)
    ja_tname.append(tname)		      #!!! need to figure out what to do here (colors and marks too)
    jmcolors.append(dcolors[j])
    jmmarks.append(imark[j])
    txvert.append(ttxvert)

    del Dts,ts,Dtx,tx,Dty,ty,ttxvert
    
shape_D = np.shape(jss)
    
    	# sort datasets by size. Plot largest first ... to smallest last.
aDs	 = np.asarray(jDs     )
aass	 = np.asarray(jss     )
aDx	 = np.asarray(jDx     )
axx	 = np.asarray(jxx     )
aDy	 = np.asarray(jDy     )
ayy	 = np.asarray(jyy     )
atxvert  = np.asarray(txvert  )
aa_tname = np.asarray(ja_tname)
amcolors = np.asarray(jmcolors)
ammarks  = np.asarray(jmmarks )

sizesM = np.nan * np.ones(shape_D[0],dtype=int)
for i in range(shape_D[0]):
  sizesM[i] = np.size(aDs[i])

sort_size = np.sort(sizesM)[::-1]	      	# [::-1] = sort in descending order. For ascending, comment out [::-1]
idx_sort = []
for i in range(shape_D[0]):  	      		# loop to go thru sort_size
  for j in range(shape_D[0]):	      		# loop to go thru size(aDs)
    if np.size(aDs[j])==sort_size[i]:
      idx_sort.append(int(j))

Ds      = aDs[idx_sort]
ss      = aass[idx_sort]
Dx      = aDx[idx_sort]
xx      = axx[idx_sort]
Dy      = aDy[idx_sort]
yy      = ayy[idx_sort]
xvert   = atxvert[idx_sort]
a_tname = aa_tname[idx_sort]
mcolors = amcolors[idx_sort]
mmarks  = ammarks[idx_sort]

del sort_size,idx_sort,shape_D,sizesM,txvert,atxvert
del jDs,jss,jDx,jxx,jDy,jyy ,ja_tname,jmcolors,jmmarks
del aDs,aass,aDx,axx,aDy,ayy,aa_tname,amcolors,ammarks

	#```````````````````````````````````````````
	# Cylindrical Equidistant Map
print("Map: Cylindrical Equidistant")
	
marksize=15

map_locations_ce(nDuniq_listALL,Ds,ss,Dx,xx,Dy,yy,x_name,a_tname,region_flag,marksize,mmarks,alphaval,mcolors,dateIN)
outname = output_path+"MAP_CE.Match_Locations."+str(dateIN)+".x_"+str(x_name)+".y"+str(y_name)+avgthin_str+".png"
plt.savefig(outname)
del outname

	#```````````````````````````````````````````
	# Orthgraphic Maps
print("Maps: Orthographic")

	# North Pole
		# lat/lon representing center point of projection (ORTHO maps)
                #       lon = 0.0   --> Polar Stereographic Projection
                #       lat = 90.0  --> North Pole
                #       lat = 0.0   --> Equator
                #       lat = -90.0 --> South Pole
central_lon = 0.0
central_lat = 90.0

map_locations_ortho(nDuniq_listALL,Ds,ss,Dx,xx,Dy,yy,x_name,a_tname,marksize,mmarks,alphaval,mcolors,dateIN,central_lon,central_lat)
outname = output_path+"MAP_Ortho.Match_Locations_NorthPole."+str(dateIN)+".x_"+str(x_name)+".y"+str(y_name)+avgthin_str+".png"
plt.savefig(outname)
del outname

	# South Pole
                # lat/lon representing center point of projection (ORTHO maps)
                #       lon = 0.0   --> Polar Stereographic Projection
                #       lat = 90.0  --> North Pole
                #       lat = 0.0   --> Equator
                #       lat = -90.0 --> South Pole
central_lon = 0.0
central_lat = -90.0

map_locations_ortho(nDuniq_listALL,Ds,ss,Dx,xx,Dy,yy,x_name,a_tname,marksize,mmarks,alphaval,mcolors,dateIN,central_lon,central_lat)
outname = output_path+"MAP_Ortho.Match_Locations_SouthPole."+str(dateIN)+".x_"+str(x_name)+".y"+str(y_name)+avgthin_str+".png"
plt.savefig(outname)
del outname

	#```````````````````````````````````````````
	# Rotating Globe
print("... Rotating Globe")
	# Locations of collocated obs on rotating map
center_lat = 0.0

outname = output_path+"MAP_Rotate.Match_Locations.centerlat"+str(center_lat)+"."+str(dateIN)+".x_"+str(x_name)+".y"+str(y_name)+avgthin_str
map_locations_ortho_rotate(output_path,outname,nDuniq_listALL,Ds,ss,Dx,xx,Dy,yy,x_name,a_tname,marksize,mmarks,alphaval,mcolors,dateIN,center_lat)
del outname

	#```````````````````````````````````````````

del y_name,Ds,ss,Dx,xx,Dy,yy,a_tname,mcolors,mmarks,xvert
del y_name_hlos,y_name_wspd,y_name_pres,y_name_hgt,shape_hlos,shape_wspd,shape_pres,shape_hgt

plt.close("all")

#==========================================================
#==========================================================
# Density Scatterplots, Zonal/Meridional Means, and Vertical Cross-sections
#==========================================================
#==========================================================
print("DENSITY SCATTERPLOTS, ZONAL/MERIDIONAL MEANS, and VERT CROSS-SECTIONS")

x_name  = driver_name
zzlabel = idx_file_vertstr

# DEPENDENT DATASETS - WIND SPEED
Vert_str = []
for j in range(tndset):
  tt=0
  if nidx[j]>0:
    tdate = Dyyyymmdd[j]
    tname = dset_name[j]
    txspd = drv_spd[j]
    tyspd = dset_spd[j]
    txlat = drv_lat[j]
    txlon = drv_lon[j]
    tylat = dset_lat[j]
    tylon = dset_lon[j]
    tVert_str = idx_file_vertstr[j]
    if tVert_str=="Pressure":
      txzp = drv_prs[j]
      tyzp = dset_prs[j]
    else:
      if np.isnan(drv_hgt[j]).all():
        txzp = prs_to_hgt(drv_prs[j])
      else:
        txzp = drv_hgt[j]
      tyzp = dset_hgt[j]
    
    ssunits = "m/s"
    if x_name.find('Aeolus')!=-1 or tname.find('Aeolus')!=-1:
      sslabelfile = "Wind_HLOS"
      ssvar = "HLOS Wind Velocity"
    else:
      sslabelfile = "Wind"
      ssvar = "Wind Speed"

    if tname == 'Aeolus':
      tname += "_"+aeolus_wind_type+"_"+aeolus_dset_type_str

    # for 3D maps
    if tVert_str=="Pressure":
      txvert = txzp
      tyvert = tyzp
      vertvarstr = tVert_str+" (hPa)"
      fvert = fpres
    else:
      if np.isnan(txzp).all():
        txvert = prs_to_hgt(txzp)
        tyvert = prs_to_hgt(tyzp)
      else:
        txvert = txzp
        tyvert = tyzp
      fvert = fhgts
      vertvarstr = tVert_str+" (km)"

    # WIND SPEED strings
    sslabel = ssvar+" ("+ssunits+")"
    zzlabel = tVert_str
        
   #````````````````````````````
   # Zonal Means (VERT VS LAT)

	# CONVERT 1D arrays to 2D for plotting
    xaxis    = flats
    yaxis    = fvert
    xaxisstr = "Latitude"
    yaxisstr = "Vertical"

		# Dv2d = DRIVER on 2d grid
		# vv2d = DEPENDENT on 2d grid
		# diff = Mean of individual differences (DEPENDENT - DRIVER) per grid cell on 2d grid
		# SDdiff = SD of individual differences (DEPENDENT - DRIVER) per grid cell on 2d grid
    Dv2d,vv2d,diff,SDdiff = var_to_2d(xaxisstr,yaxisstr,xaxis,yaxis,txlat,txvert,txspd,tyspd)

	# PLOTS
    plottype = "2d"
    xaxisstr = "Latitude"
    yaxisstr = vertvarstr
   
    level    = " "
    levunits = " "
    regionstr = " "

    		#````````````````````````````
		# Wind
		# ... DRIVER (Dv2d)
    contour2d(Dv2d,ssvar,xaxis,yaxis,xaxisstr,yaxisstr,x_name+" matched with "+tname,ssunits,plottype,level,levunits,regionstr)
    outname = output_path+"VERTvsLAT_ZonalMean."+str(sslabelfile)+".DRV."+str(dateIN)+".x_"+str(x_name)+".y_"+str(tname)+avgthin_str+".png"
    plt.savefig(outname)
    del outname
    if ssvar.find("HLOS")!=-1:
        		# plot absolute value (ABS) of HLOS wind velocity
      contour2d(abs(Dv2d),ssvar,xaxis,yaxis,xaxisstr,yaxisstr,x_name+" matched with "+tname,ssunits,plottype,level,levunits,regionstr)
      outname = output_path+"VERTvsLAT_ZonalMean."+str(sslabelfile)+"_ABS.DRV."+str(dateIN)+".x_"+str(x_name)+".y_"+str(tname)+avgthin_str+".png"
      plt.savefig(outname)
      del outname    
 		# ... DEPENDENT (vv2d)
    contour2d(vv2d,ssvar,xaxis,yaxis,xaxisstr,yaxisstr,tname,ssunits,plottype,level,levunits,regionstr)
    outname = output_path+"VERTvsLAT_ZonalMean."+str(ssvar)+".DEP."+str(dateIN)+".x_"+str(x_name)+".y_"+str(tname)+avgthin_str+".png"
    plt.savefig(outname)
    del outname
    if ssvar.find("HLOS")!=-1:
        		# plot absolute value (ABS) of HLOS wind velocity
      contour2d(abs(vv2d),ssvar,xaxis,yaxis,xaxisstr,yaxisstr,tname,ssunits,plottype,level,levunits,regionstr)
      outname = output_path+"VERTvsLAT_ZonalMean."+str(sslabelfile)+"_ABS.DEP."+str(dateIN)+".x_"+str(x_name)+".y_"+str(tname)+avgthin_str+".png"
      plt.savefig(outname)
      del outname 
 		# ... Difference (diff): DEPENDENT - DRIVER  
    contour2d(diff,ssvar+" Diff",xaxis,yaxis,xaxisstr,yaxisstr,tname+"-"+x_name,ssunits,plottype,level,levunits,regionstr)
    outname = output_path+"VERTvsLAT_ZonalMean."+str(ssvar)+"_Diff.DEP-DRV."+str(dateIN)+".x_"+str(x_name)+".y_"+str(tname)+avgthin_str+".png"
    plt.savefig(outname)
    del outname
		# ... SD of Difference (SDdiff): DEPENDENT - DRIVER
    contour2d(SDdiff,ssvar+" SD",xaxis,yaxis,xaxisstr,yaxisstr,tname+"-"+x_name,ssunits,plottype,level,levunits,regionstr)
    outname = output_path+"VERTvsLAT_ZonalMean."+str(ssvar)+"_SDdiff.DEP-DRV."+str(dateIN)+".x_"+str(x_name)+".y_"+str(tname)+avgthin_str+".png"
    plt.savefig(outname)
    del outname

    del Dv2d,vv2d,diff,SDdiff

    		#````````````````````````````
    		# Number Density
    drv_count,dep_count = var_to_2d_counts(xaxisstr,yaxisstr,xaxis,yaxis,txlat,txvert,txlat,txvert)
    		# ... DRIVER (drv_count)
    contour2d(drv_count,"Count",xaxis,yaxis,xaxisstr,yaxisstr,x_name+" matched with "+tname,ssunits,plottype,level,levunits,regionstr)
    outname = output_path+"VERTvsLAT_ZonalMean.Nobs_Density.DRV."+str(dateIN)+".x_"+str(x_name)+".y_"+str(tname)+avgthin_str+".png"
    plt.savefig(outname)
    del outname
    		# ... DEPENDENT (dep_count)
    contour2d(dep_count,"Count",xaxis,yaxis,xaxisstr,yaxisstr,tname,ssunits,plottype,level,levunits,regionstr)
    outname = output_path+"VERTvsLAT_ZonalMean.Nobs_Density.DEP."+str(dateIN)+".x_"+str(x_name)+".y_"+str(tname)+avgthin_str+".png"
    plt.savefig(outname)
    del outname
    
    del drv_count,dep_count,xaxis,yaxis,xaxisstr,yaxisstr
    
    #````````````````````````````
    # Meridional Means (VERT VS LON)

	# CONVERT 1D arrays to 2D for plotting
    xaxis    = flons180
    yaxis    = fvert
    xaxisstr = "Longitude"
    yaxisstr = "Vertical"

		# convert [0,360] longitudes to [-180,180] for plotting here
    txlon180 = [txlon[i]-360.0 if txlon[i]>180.0 else txlon[i] for i in range(np.size(txlon))]
    tylon180 = [tylon[i]-360.0 if tylon[i]>180.0 else tylon[i] for i in range(np.size(tylon))]

		# Dv2d = DRIVER on 2d grid
		# vv2d = DEPENDENT on 2d grid
		# diff = Mean of individual differences (DEPENDENT - DRIVER) per grid cell on 2d grid
		# SDdiff = SD of individual differences (DEPENDENT - DRIVER) per grid cell on 2d grid
    Dv2d,vv2d,diff,SDdiff = var_to_2d(xaxisstr,yaxisstr,xaxis,yaxis,txlon180,txvert,txspd,tyspd)

	# PLOTS
    plottype = "2d"
    xaxisstr = "Longitude"
    yaxisstr = vertvarstr
   
    level    = " "
    levunits = " "
    regionstr = " "

    		#````````````````````````````
		# Wind
		# ... DRIVER (Dv2d)
    contour2d(Dv2d,ssvar,xaxis,yaxis,xaxisstr,yaxisstr,x_name+" matched with "+tname,ssunits,plottype,level,levunits,regionstr)
    outname = output_path+"VERTvsLON_MeridMean."+str(sslabelfile)+".DRV."+str(dateIN)+".x_"+str(x_name)+".y_"+str(tname)+avgthin_str+".png"
    plt.savefig(outname)
    del outname
    if ssvar.find("HLOS")!=-1:
        # plot absolute value (ABS) of HLOS wind velocity
      contour2d(abs(Dv2d),ssvar,xaxis,yaxis,xaxisstr,yaxisstr,x_name+" matched with "+tname,ssunits,plottype,level,levunits,regionstr)
      outname = output_path+"VERTvsLON_MeridMean."+str(sslabelfile)+"_ABS.DRV."+str(dateIN)+".x_"+str(x_name)+".y_"+str(tname)+avgthin_str+".png"
      plt.savefig(outname)
      del outname
 		# ... DEPENDENT (vv2d)
    contour2d(vv2d,ssvar,xaxis,yaxis,xaxisstr,yaxisstr,tname,ssunits,plottype,level,levunits,regionstr)
    outname = output_path+"VERTvsLON_MeridMean."+str(sslabelfile)+".DEP."+str(dateIN)+".x_"+str(x_name)+".y_"+str(tname)+avgthin_str+".png"
    plt.savefig(outname)
    del outname
    if ssvar.find("HLOS")!=-1:
        # plot absolute value (ABS) of HLOS wind velocity
      contour2d(abs(vv2d),ssvar,xaxis,yaxis,xaxisstr,yaxisstr,tname,ssunits,plottype,level,levunits,regionstr)
      outname = output_path+"VERTvsLON_MeridMean."+str(sslabelfile)+"_ABS.DEP."+str(dateIN)+".x_"+str(x_name)+".y_"+str(tname)+avgthin_str+".png"
      plt.savefig(outname)
      del outname
 		# ... Difference (diff): DEPENDENT - DRIVER  
    contour2d(diff,ssvar+" Diff",xaxis,yaxis,xaxisstr,yaxisstr,tname+"-"+x_name,ssunits,plottype,level,levunits,regionstr)
    outname = output_path+"VERTvsLON_MeridMean."+str(sslabelfile)+"_Diff.DEP-DRV."+str(dateIN)+".x_"+str(x_name)+".y_"+str(tname)+avgthin_str+".png"
    plt.savefig(outname)
    del outname
		# ... SD of Difference (SDdiff): DEPENDENT - DRIVER
    contour2d(SDdiff,ssvar+" SD",xaxis,yaxis,xaxisstr,yaxisstr,tname+"-"+x_name,ssunits,plottype,level,levunits,regionstr)
    outname = output_path+"VERTvsLON_MeridMean."+str(sslabelfile)+"_SDdiff.DEP-DRV."+str(dateIN)+".x_"+str(x_name)+".y_"+str(tname)+avgthin_str+".png"
    plt.savefig(outname)
    del outname

    del Dv2d,vv2d,diff,SDdiff

    		#````````````````````````````
    		# Number Density
    drv_count,dep_count = var_to_2d_counts(xaxisstr,yaxisstr,xaxis,yaxis,txlon180,txvert,txlon180,txvert)
    		# ... DRIVER (drv_count)
    contour2d(drv_count,"Count",xaxis,yaxis,xaxisstr,yaxisstr,x_name+" matched with "+tname,ssunits,plottype,level,levunits,regionstr)
    outname = output_path+"VERTvsLON_MeridMean.Nobs_Density.DRV."+str(dateIN)+".x_"+str(x_name)+".y_"+str(tname)+avgthin_str+".png"
    plt.savefig(outname)
    del outname
    		# ... DEPENDENT (dep_count)
    contour2d(dep_count,"Count",xaxis,yaxis,xaxisstr,yaxisstr,tname,ssunits,plottype,level,levunits,regionstr)
    outname = output_path+"VERTvsLON_MeridMean.Nobs_Density.DEP."+str(dateIN)+".x_"+str(x_name)+".y_"+str(tname)+avgthin_str+".png"
    plt.savefig(outname)
    del outname
    
    del drv_count,dep_count,xaxis,yaxis,xaxisstr,yaxisstr
    del txlon180,tylon180

    #````````````````````````````
    # Vertical Cross-sections at specified longitudes (VERT VS LAT plot)

		# convert 1D arrays to 2D for plotting
    xaxis    = flats
    yaxis    = fvert
    xaxisstr = "Latitude"
    yaxisstr = "Vertical"

    plottype = "2d"
    xaxisstr = "Latitude"
    yaxisstr = vertvarstr
   
    level    = " "
    levunits = " "
    regionstr = " "

    for ilon in range(np.size(vertcrossLONSmin)):
      xidxLON = np.where((txlon>=vertcrossLONSmin[ilon])*(txlon<vertcrossLONSmax[ilon]))
      yidxLON = np.where((tylon>=vertcrossLONSmin[ilon])*(tylon<vertcrossLONSmax[ilon]))

                # Dv2d = DRIVER on 2d grid
                # vv2d = DEPENDENT on 2d grid
                # diff = Mean of individual differences (DEPENDENT - DRIVER) per grid cell on 2d grid
                # SDdiff = SD of individual differences (DEPENDENT - DRIVER) per grid cell on 2d grid
      Dv2d,vv2d,diff,SDdiff  = var_to_2d(xaxisstr,yaxisstr,xaxis,yaxis,txlat[xidxLON],txvert[xidxLON],txspd[xidxLON],tyspd[xidxLON])

      del xidxLON,yidxLON

    		#````````````````````````````
		# Wind
		# ... DRIVER (Dv2d)
      contour2d(Dv2d,ssvar,xaxis,yaxis,xaxisstr,yaxisstr,x_name+" matched with "+tname,ssunits,plottype,level,levunits,regionstr)
      outname = output_path+"VERTvsLAT_LonRange"+str(vertcrossLONSmin[ilon])+"-"+str(vertcrossLONSmax[ilon])+"."+str(sslabelfile)+".DRV."+str(dateIN)+".x_"+str(x_name)+".y_"+str(tname)+avgthin_str+".png"
      plt.savefig(outname)
      del outname
      if ssvar.find("HLOS")!=-1:
        		# plot absolute value (ABS) of HLOS wind velocity
        contour2d(abs(Dv2d),ssvar,xaxis,yaxis,xaxisstr,yaxisstr,x_name+" matched with "+tname,ssunits,plottype,level,levunits,regionstr)
        outname = output_path+"VERTvsLAT_LonRange"+str(vertcrossLONSmin[ilon])+"-"+str(vertcrossLONSmax[ilon])+"."+str(sslabelfile)+"_ABS.DRV."+str(dateIN)+".x_"+str(x_name)+".y_"+str(tname)+avgthin_str+".png"
        plt.savefig(outname)
        del outname
 		# ... DEPENDENT (vv2d)
      contour2d(vv2d,ssvar,xaxis,yaxis,xaxisstr,yaxisstr,tname,ssunits,plottype,level,levunits,regionstr)
      outname = output_path+"VERTvsLAT_LonRange"+str(vertcrossLONSmin[ilon])+"-"+str(vertcrossLONSmax[ilon])+"."+str(sslabelfile)+".DEP."+str(dateIN)+".x_"+str(x_name)+".y_"+str(tname)+avgthin_str+".png"
      plt.savefig(outname)
      del outname
      if ssvar.find("HLOS")!=-1:
        		# plot absolute value (ABS) of HLOS wind velocity
        contour2d(abs(vv2d),ssvar,xaxis,yaxis,xaxisstr,yaxisstr,tname,ssunits,plottype,level,levunits,regionstr)
        outname = output_path+"VERTvsLAT_LonRange"+str(vertcrossLONSmin[ilon])+"-"+str(vertcrossLONSmax[ilon])+"."+str(sslabelfile)+"_ABS.DEP."+str(dateIN)+".x_"+str(x_name)+".y_"+str(tname)+avgthin_str+".png"
        plt.savefig(outname)
        del outname
 		# ... Difference (diff): DEPENDENT - DRIVER  
      contour2d(diff,ssvar+" Diff",xaxis,yaxis,xaxisstr,yaxisstr,tname+"-"+x_name,ssunits,plottype,level,levunits,regionstr)
      outname = output_path+"VERTvsLAT_LonRange"+str(vertcrossLONSmin[ilon])+"-"+str(vertcrossLONSmax[ilon])+"."+str(sslabelfile)+"_Diff.DEP-DRV."+str(dateIN)+".x_"+str(x_name)+".y_"+str(tname)+avgthin_str+".png"
      plt.savefig(outname)
      del outname
      		# ... SD of Difference (diff): DEPENDENT - DRIVER  
      contour2d(SDdiff,ssvar+" SD",xaxis,yaxis,xaxisstr,yaxisstr,tname+"-"+x_name,ssunits,plottype,level,levunits,regionstr)
      outname = output_path+"VERTvsLAT_LonRange"+str(vertcrossLONSmin[ilon])+"-"+str(vertcrossLONSmax[ilon])+"."+str(sslabelfile)+"_SDdiff.DEP-DRV."+str(dateIN)+".x_"+str(x_name)+".y_"+str(tname)+avgthin_str+".png"
      plt.savefig(outname)
      del outname
    
      del Dv2d,vv2d,diff,SDdiff
    
    		#````````````````````````````
    		# Number Density
      drv_count,dep_count = var_to_2d_counts(xaxisstr,yaxisstr,xaxis,yaxis,txlat,txvert,txlat,txvert)
    		# ... DRIVER (drv_count)
      contour2d(drv_count,"Count",xaxis,yaxis,xaxisstr,yaxisstr,x_name+" matched with "+tname,ssunits,plottype,level,levunits,regionstr)
      outname = output_path+"VERTvsLAT_LonRange"+str(vertcrossLONSmin[ilon])+"-"+str(vertcrossLONSmax[ilon])+".Nobs_Density.DRV."+str(dateIN)+".x_"+str(x_name)+".y_"+str(tname)+avgthin_str+".png"
      plt.savefig(outname)
      del outname
    		# ... DEPENDENT (dep_count)
      contour2d(dep_count,"Count",xaxis,yaxis,xaxisstr,yaxisstr,tname,ssunits,plottype,level,levunits,regionstr)
      outname = output_path+"VERTvsLAT_LonRange"+str(vertcrossLONSmin[ilon])+"-"+str(vertcrossLONSmax[ilon])+".Nobs_Density.DEP."+str(dateIN)+".x_"+str(x_name)+".y_"+str(tname)+avgthin_str+".png"
      plt.savefig(outname)
      del outname
    
      del drv_count,dep_count

      plt.close("all")

    del xaxis,yaxis,xaxisstr,yaxisstr

    #````````````````````````````
    # VERT vs TIME

		# convert 1D arrays to 2D for plotting
    xaxis    = date_uniq
    yaxis    = fvert
    xaxisstr = "Time"
    yaxisstr = vertvarstr

    	#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    	# Global
    	#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    plottype = "2d"
    level    = " "
    levunits = " "
    regionstr = "Global"
     
    		#````````````````````````````
    		# Number Density
    drv_count,dep_count = var_to_2d_counts(xaxisstr,yaxisstr,xaxis,yaxis,tdate,txvert,tdate,txvert)
    		# ... DRIVER (drv_count)
    contour2d(drv_count,"Count",xaxis,yaxis,xaxisstr,yaxisstr,x_name+" matched with "+tname,ssunits,plottype,level,levunits,regionstr)
    outname = output_path+"VERTvsTIME.Nobs_Density.DRV."+regionstr+"."+str(dateIN)+".x_"+str(x_name)+".y_"+str(tname)+avgthin_str+".png"
    plt.savefig(outname)
    del outname
    		# ... DEPENDENT (dep_count)
    contour2d(dep_count,"Count",xaxis,yaxis,xaxisstr,yaxisstr,tname,ssunits,plottype,level,levunits,regionstr)
    outname = output_path+"VERTvsTIME.Nobs_Density.DEP."+regionstr+"."+str(dateIN)+".x_"+str(x_name)+".y_"+str(tname)+avgthin_str+".png"
    plt.savefig(outname)
    del outname
    
    del drv_count,dep_count#,xaxis,yaxis,xaxisstr,yaxisstr
    del regionstr

    plt.close("all")
     
    	#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	# Regional (NH,TR,SH)
	#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    for ir in range(nregions):
      RidxD = np.where((txlat>=regions_latmin[ir])*(txlat<regions_latmax[ir]))

      str2d = " "

      plottype = "2d"
      level    = " "
      levunits = " "
      regionstr = str(regions[ir])

    		#````````````````````````````
    		# Number Density
      drv_count,dep_count = var_to_2d_counts(xaxisstr,yaxisstr,xaxis,yaxis,tdate[RidxD],txvert[RidxD],tdate[RidxD],txvert[RidxD])
    		# ... DRIVER (drv_count)
      contour2d(drv_count,"Count",xaxis,yaxis,xaxisstr,yaxisstr,x_name+" matched with "+tname,ssunits,plottype,level,levunits,regionstr)
      outname = output_path+"VERTvsTIME.Nobs_Density.DRV."+regionstr+"."+str(dateIN)+".x_"+str(x_name)+".y_"+str(tname)+avgthin_str+".png"
      plt.savefig(outname)
      del outname
    		# ... DEPENDENT (dep_count)
      contour2d(dep_count,"Count",xaxis,yaxis,xaxisstr,yaxisstr,tname,ssunits,plottype,level,levunits,regionstr)
      outname = output_path+"VERTvsTIME.Nobs_Density.DEP."+regionstr+"."+str(dateIN)+".x_"+str(x_name)+".y_"+str(tname)+avgthin_str+".png"
      plt.savefig(outname)
      del outname
    
      del drv_count,dep_count
      del regionstr
      del RidxD
  
      plt.close("all")

    del xaxis,yaxis,xaxisstr,yaxisstr

    #````````````````````````````
    # LAT vs LON: 2D maps
	
    var2dstr = "Wind Speed"
    var2dstr = "Wind"
    if x_name.find('Aeolus')!=-1 or tname.find('Aeolus')!=-1:
      var2dstr = "HLOS Wind Velocity"
      var2dstrfile = "Wind_HLOS"
    units = "m/s"

    if tVert_str=="Pressure":
      levunits = "hPa"
      levs     = levsP		# levels to plot
      dlev     = plevstride
    elif tVert_str=="Height":
      levunits = "km"
      levs     = levsZ		# levels to plot
      dlev     = hlevstride

    daterange = dateSTART+"-"+dateEND

    xaxis    = flons
    yaxis    = flats
    xaxisstr = "Longitude"
    yaxisstr = "Latitude"

		#````````````````````````````
    		# More Number Density Maps
    drv_count,dep_count = var_to_2d_counts(xaxisstr,yaxisstr,xaxis,yaxis,txlon,txlat,txlon,txlat)

    			# copy count arrays to plot rotating globe later on in script
    drv_count_rotate = drv_count
    dep_count_rotate = dep_count

		# MAP: CE
    plottype = "map"
    regionstr = " "
   
		# ... DRIVER
    contour2d(drv_count,"Count",xaxis,yaxis,xaxisstr,yaxisstr,x_name+" matched with "+tname,ssunits,plottype,level,levunits,regionstr)
    outname = output_path+"MAP_2Dcontour.Nobs_Density.DRV."+str(dateIN)+".x_"+str(x_name)+".y_"+str(tname)+avgthin_str+".png"
    plt.savefig(outname)
    del outname
        	# ... DEPENDENT
    contour2d(dep_count,"Count",xaxis,yaxis,xaxisstr,yaxisstr,tname,ssunits,plottype,level,levunits,regionstr)
    outname = output_path+"MAP_2Dcontour.Nobs_Density.DEP."+str(dateIN)+".x_"+str(x_name)+".y_"+str(tname)+avgthin_str+".png"
    plt.savefig(outname)
    del outname

    		# MAP: Orthographic
			# North Pole
    center_lon = 0.0
    center_lat = 90.0
    			# ... DRIVER
    contour2d_orthomap(daterange,drv_count,"Count",xaxis,yaxis,xaxisstr,yaxisstr,x_name+" matched with "+tname,ssunits,center_lon,center_lat,level,levunits)
    outname = output_path+"MAP_Ortho_NorthPole.Nobs_Density.DRV."+str(dateIN)+".x_"+str(x_name)+".y_"+str(tname)+avgthin_str+".png"
    plt.savefig(outname)
    del outname
                	# ... DEPENDENT
    contour2d_orthomap(daterange,dep_count,"Count",xaxis,yaxis,xaxisstr,yaxisstr,tname,ssunits,center_lon,center_lat,level,levunits)
    outname = output_path+"MAP_Ortho_NorthPole.Nobs_Density.DEP."+str(dateIN)+".x_"+str(x_name)+".y_"+str(tname)+avgthin_str+".png"
    plt.savefig(outname)
    del outname

                        # South Pole
    center_lon = 0.0
    center_lat = -90.0
                        # ... DRIVER
    contour2d_orthomap(daterange,drv_count,"Count",xaxis,yaxis,xaxisstr,yaxisstr,x_name+" matched with "+tname,ssunits,center_lon,center_lat,level,levunits)
    outname = output_path+"MAP_Ortho_SouthPole.Nobs_Density.DRV."+str(dateIN)+".x_"+str(x_name)+".y_"+str(tname)+avgthin_str+".png"
    plt.savefig(outname)
    del outname
                        # ... DEPENDENT
    contour2d_orthomap(daterange,dep_count,"Count",xaxis,yaxis,xaxisstr,yaxisstr,tname,ssunits,center_lon,center_lat,level,levunits)
    outname = output_path+"MAP_Ortho_SouthPole.Nobs_Density.DEP."+str(dateIN)+".x_"+str(x_name)+".y_"+str(tname)+avgthin_str+".png"
    plt.savefig(outname)
    del outname
    
    del drv_count,dep_count
    
    plt.close("all")

		#````````````````````````````
    		# Mapped Means at Specified Levels
		
    Dv2d_rotate = []
    vv2d_rotate = []
    diff_rotate = []
    SDdiff_rotate = []
    levs_rotate = []

    for ip in range(np.size(levs)):
	# make 2d variables
      hmin = levs[ip]-dlev
      hmax = levs[ip]+dlev
      Lidx = np.where((txvert>=hmin)*(txvert<hmax))

                # Dv2d = DRIVER on 2d grid
                # vv2d = DEPENDENT on 2d grid
                # diff = Mean of individual differences (DEPENDENT - DRIVER) per grid cell on 2d grid
                # SDdiff = SD of individual differences (DEPENDENT - DRIVER) per grid cell on 2d grid
      Dv2d,vv2d,diff,SDdiff = var_to_2d(xaxisstr,yaxisstr,xaxis,yaxis,txlon[Lidx],txlat[Lidx],txspd[Lidx],tyspd[Lidx])

      del Lidx,hmin,hmax

	# copy arrays for plotting on rotating globe later on in script
      Dv2d_rotate.append(Dv2d)
      vv2d_rotate.append(vv2d)
      diff_rotate.append(diff)
      SDdiff_rotate.append(SDdiff)
      levs_rotate.append(int(levs[ip]))

      plottype  = "map"
      level     = str(int(levs[ip]))
      regionstr = " "

	# ... DRIVER (x_grid)
      contour2d(Dv2d,var2dstr,xaxis,yaxis,xaxisstr,yaxisstr,x_name+" matched with "+tname,ssunits,plottype,level,levunits,regionstr)
      outname = output_path+"MAP_2Dcontour."+str(var2dstrfile)+".DRV."+str(level)+str(levunits)+"."+str(dateIN)+".x_"+str(x_name)+".y_"+str(tname)+avgthin_str+".png"
      plt.savefig(outname)
      del outname
      if var2dstr.find("HLOS")!=-1:
        	# plot absolute value (ABS) of HLOS wind velocity
        contour2d(abs(Dv2d),var2dstr,xaxis,yaxis,xaxisstr,yaxisstr,x_name+" matched with "+tname,ssunits,plottype,level,levunits,regionstr)
        outname = output_path+"MAP_2Dcontour."+str(var2dstrfile)+"_ABS.DRV."+str(level)+str(levunits)+"."+str(dateIN)+".x_"+str(x_name)+".y_"+str(tname)+avgthin_str+".png"
        plt.savefig(outname)
        del outname
	# ... DEPENDENT (y_grid)
      contour2d(vv2d,var2dstr,xaxis,yaxis,xaxisstr,yaxisstr,tname,ssunits,plottype,level,levunits,regionstr)
      outname = output_path+"MAP_2Dcontour."+str(var2dstrfile)+".DEP."+str(level)+str(levunits)+"."+str(dateIN)+".x_"+str(x_name)+".y_"+str(tname)+avgthin_str+".png"
      plt.savefig(outname)
      del outname
      if var2dstr.find("HLOS")!=-1:
        	# plot absolute value (ABS) of HLOS wind velocity
        contour2d(abs(vv2d),var2dstr,xaxis,yaxis,xaxisstr,yaxisstr,tname,ssunits,plottype,level,levunits,regionstr)
        outname = output_path+"MAP_2Dcontour."+str(var2dstrfile)+"_ABS.DEP."+str(level)+str(levunits)+"."+str(dateIN)+".x_"+str(x_name)+".y_"+str(tname)+avgthin_str+".png"
        plt.savefig(outname)
        del outname
	# ... Difference (diff): DEPENDENT - DRIVER
      contour2d(diff,var2dstr+" Diff",xaxis,yaxis,xaxisstr,yaxisstr,tname+"-"+x_name,ssunits,plottype,level,levunits,regionstr)
      outname = output_path+"MAP_2Dcontour."+str(var2dstrfile)+"_Diff.DEP-DRV."+str(level)+str(levunits)+"."+str(dateIN)+".x_"+str(x_name)+".y_"+str(tname)+avgthin_str+".png"
      plt.savefig(outname)
      del outname
        # ... SD of Difference (SDdiff): DEPENDENT - DRIVER
      contour2d(SDdiff,var2dstr+" SD",xaxis,yaxis,xaxisstr,yaxisstr,tname+"-"+x_name,ssunits,plottype,level,levunits,regionstr)
      outname = output_path+"MAP_2Dcontour."+str(var2dstrfile)+"_SDdiff.DEP-DRV."+str(level)+str(levunits)+"."+str(dateIN)+".x_"+str(x_name)+".y_"+str(tname)+avgthin_str+".png"
      plt.savefig(outname)
      del outname

      plt.close("all")

	# MAP: Orthographic
	# ... North Pole
      center_lon = 0.0
      center_lat = 90.0
    		# ... DRIVER
      contour2d_orthomap(daterange,Dv2d,var2dstr,xaxis,yaxis,xaxisstr,yaxisstr,x_name+" matched with "+tname,ssunits,center_lon,center_lat,level,levunits)
      outname = output_path+"MAP_Ortho_NorthPole."+str(var2dstrfile)+".DRV."+str(level)+str(levunits)+"."+str(dateIN)+".x_"+str(x_name)+".y_"+str(tname)+avgthin_str+".png"
      plt.savefig(outname)
      del outname
      if var2dstr.find("HLOS")!=-1:
        		# plot absolute value (ABS) of HLOS wind velocity
        contour2d_orthomap(daterange,abs(Dv2d),var2dstr,xaxis,yaxis,xaxisstr,yaxisstr,x_name+" matched with "+tname,ssunits,center_lon,center_lat,level,levunits)
        outname = output_path+"MAP_Ortho_NorthPole."+str(var2dstrfile)+"_ABS.DRV."+str(level)+str(levunits)+"."+str(dateIN)+".x_"+str(x_name)+".y_"+str(tname)+avgthin_str+".png"
        plt.savefig(outname)
        del outname
      		# ... DEPENDENT
      contour2d_orthomap(daterange,vv2d,var2dstr,xaxis,yaxis,xaxisstr,yaxisstr,tname,ssunits,center_lon,center_lat,level,levunits)
      outname = output_path+"MAP_Ortho_NorthPole."+str(var2dstrfile)+".DEP."+str(level)+str(levunits)+"."+str(dateIN)+".x_"+str(x_name)+".y_"+str(tname)+avgthin_str+".png"
      plt.savefig(outname)
      del outname
      if var2dstr.find("HLOS")!=-1:
        		# plot absolute value (ABS) of HLOS wind velocity
        contour2d_orthomap(daterange,abs(vv2d),var2dstr,xaxis,yaxis,xaxisstr,yaxisstr,tname,ssunits,center_lon,center_lat,level,levunits)
        outname = output_path+"MAP_Ortho_NorthPole."+str(var2dstrfile)+"_ABS.DEP."+str(level)+str(levunits)+"."+str(dateIN)+".x_"+str(x_name)+".y_"+str(tname)+avgthin_str+".png"
        plt.savefig(outname)
        del outname
      		# ... Difference: DEPENDENT - DRIVER
      contour2d_orthomap(daterange,diff,var2dstr+" Diff",xaxis,yaxis,xaxisstr,yaxisstr,tname+"-"+x_name,ssunits,center_lon,center_lat,level,levunits)
      outname = output_path+"MAP_Ortho_NorthPole."+str(var2dstrfile)+"_Diff.DEP-DRV."+str(level)+str(levunits)+"."+str(dateIN)+".x_"+str(x_name)+".y_"+str(tname)+avgthin_str+".png"
      plt.savefig(outname)
      del outname
		# ... SD of Difference: DEPENDENT - DRIVER
      contour2d_orthomap(daterange,SDdiff,var2dstr+" SD",xaxis,yaxis,xaxisstr,yaxisstr,tname+"-"+x_name,ssunits,center_lon,center_lat,level,levunits)
      outname = output_path+"MAP_Ortho_NorthPole."+str(var2dstrfile)+"_SDdiff.DEP-DRV."+str(level)+str(levunits)+"."+str(dateIN)+".x_"+str(x_name)+".y_"+str(tname)+avgthin_str+".png"
      plt.savefig(outname)
      del outname

      plt.close("all")

      	# ... South Pole
      center_lon = 0.0
      center_lat = -90.0
    		# ... DRIVER
      contour2d_orthomap(daterange,Dv2d,var2dstr,xaxis,yaxis,xaxisstr,yaxisstr,x_name+" matched with "+tname,ssunits,center_lon,center_lat,level,levunits)
      outname = output_path+"MAP_Ortho_SouthPole."+str(var2dstrfile)+".DRV."+str(level)+str(levunits)+"."+str(dateIN)+".x_"+str(x_name)+".y_"+str(tname)+avgthin_str+".png"
      plt.savefig(outname)
      del outname
      if var2dstr.find("HLOS")!=-1:
        		# plot absolute value (ABS) of HLOS wind velocity
        contour2d_orthomap(daterange,abs(Dv2d),var2dstr,xaxis,yaxis,xaxisstr,yaxisstr,x_name+" matched with "+tname,ssunits,center_lon,center_lat,level,levunits)
        outname = output_path+"MAP_Ortho_SouthPole."+str(var2dstrfile)+"_ABS.DRV."+str(level)+str(levunits)+"."+str(dateIN)+".x_"+str(x_name)+".y_"+str(tname)+avgthin_str+".png"
        plt.savefig(outname)
        del outname
      		# ... DEPENDENT
      contour2d_orthomap(daterange,vv2d,var2dstr,xaxis,yaxis,xaxisstr,yaxisstr,tname,ssunits,center_lon,center_lat,level,levunits)
      outname = output_path+"MAP_Ortho_SouthPole."+str(var2dstrfile)+".DEP."+str(level)+str(levunits)+"."+str(dateIN)+".x_"+str(x_name)+".y_"+str(tname)+avgthin_str+".png"
      plt.savefig(outname)
      del outname
      if var2dstr.find("HLOS")!=-1:
        		# plot absolute value (ABS) of HLOS wind velocity
        contour2d_orthomap(daterange,abs(vv2d),var2dstr,xaxis,yaxis,xaxisstr,yaxisstr,tname,ssunits,center_lon,center_lat,level,levunits)
        outname = output_path+"MAP_Ortho_SouthPole."+str(var2dstrfile)+"_ABS.DEP."+str(level)+str(levunits)+"."+str(dateIN)+".x_"+str(x_name)+".y_"+str(tname)+avgthin_str+".png"
        plt.savefig(outname)
        del outname
      		# ... Difference: DEPENDENT - DRIVER
      contour2d_orthomap(daterange,diff,var2dstr+" Diff",xaxis,yaxis,xaxisstr,yaxisstr,tname+"-"+x_name,ssunits,center_lon,center_lat,level,levunits)
      outname = output_path+"MAP_Ortho_SouthPole."+str(var2dstrfile)+"_Diff.DEP-DRV."+str(level)+str(levunits)+"."+str(dateIN)+".x_"+str(x_name)+".y_"+str(tname)+avgthin_str+".png"
      plt.savefig(outname)
      del outname
		# ... SD of Difference: DEPENDENT - DRIVER
      contour2d_orthomap(daterange,SDdiff,var2dstr+" SD",xaxis,yaxis,xaxisstr,yaxisstr,tname+"-"+x_name,ssunits,center_lon,center_lat,level,levunits)
      outname = output_path+"MAP_Ortho_SouthPole."+str(var2dstrfile)+"_SDdiff.DEP-DRV."+str(level)+str(levunits)+"."+str(dateIN)+".x_"+str(x_name)+".y_"+str(tname)+avgthin_str+".png"
      plt.savefig(outname)
      del outname     

      plt.close("all")

      del level
      del Dv2d,vv2d,diff,SDdiff
    del xaxis,yaxis,xaxisstr,yaxisstr,levs,levunits

    		#````````````````````````````
		# 2D maps of scattered points
	
    xaxis    = flons
    yaxis    = flats

    if tVert_str=="Pressure":
      levunits  = "hPa"
      levs     = levsP		# levels to plot
    elif tVert_str=="Height":
      levunits = "km"
      levs     = levsZ						# levels to plot

	# plot all observations pairs within +/-levstride of each vertical level in 'levs'
    xtmp_var  = txspd
    xtmp_vert = txvert
    xtmp_lat  = txlat
    xtmp_lon  = txlon
    ytmp_var  = tyspd
    ytmp_vert = tyvert
    ytmp_lat  = tylat
    ytmp_lon  = tylon
    xmlat = [[] for n in range(np.size(levs))]
    xmlon = [[] for n in range(np.size(levs))]
    xmvar = [[] for n in range(np.size(levs))]
    ymlat = [[] for n in range(np.size(levs))]
    ymlon = [[] for n in range(np.size(levs))]
    ymvar = [[] for n in range(np.size(levs))]
    diffvar = [[] for n in range(np.size(levs))]
    for ip in range(np.size(levs)):
      if tVert_str=="Pressure":
        levdiff = plevstride                  # half the difference between each gridded z value
      elif tVert_str=="Height":
        levdiff = hlevstride                  # half the difference between each gridded z value
      levmin  = levs[ip]-levdiff
      levmax  = levs[ip]+levdiff
      del levdiff
      levarr = xtmp_vert
      Lidx = np.where((levarr>=levmin)*(levarr<levmax))
      xmlat[ip] = np.append(xmlat[ip],xtmp_lat[Lidx],axis=0)
      xmlon[ip] = np.append(xmlon[ip],xtmp_lon[Lidx],axis=0)
      xmvar[ip] = np.append(xmvar[ip],xtmp_var[Lidx],axis=0)
      ymlat[ip] = np.append(ymlat[ip],ytmp_lat[Lidx],axis=0)
      ymlon[ip] = np.append(ymlon[ip],ytmp_lon[Lidx],axis=0)
      ymvar[ip] = np.append(ymvar[ip],ytmp_var[Lidx],axis=0)
      txx = np.asarray(xtmp_var[Lidx])
      tyy = np.asarray(ytmp_var[Lidx])
      diff = tyy - txx
      diffvar[ip] = np.append(diffvar[ip],diff,axis=0)
      del Lidx,levarr
      del txx,tyy,diff
    del xtmp_var,xtmp_vert,xtmp_lat,xtmp_lon
    del ytmp_var,ytmp_vert,ytmp_lat,ytmp_lon

    opt=0
    marksize=15

    for ip in range(np.size(levs)):
      if np.size(xmvar[ip])>0:
        level = str(int(levs[ip]))
        # ... DRIVER
        map_points2d_ce(var2dstr,xmlon[ip],xmlat[ip],xmvar[ip],xaxis,yaxis,x_name+" matched with "+tname,marksize,alphaval,daterange,units,levs[ip],levunits,opt)
        outname = output_path+"MAP_2Dpoints."+str(var2dstrfile)+".DRV."+level+str(levunits)+"."+str(dateIN)+".x_"+str(x_name)+".y_"+str(tname)+avgthin_str+".png"
        plt.savefig(outname)
        del outname
        if var2dstr.find("HLOS")!=-1:
        	# plot absolute value (ABS) of HLOS wind velocity
          map_points2d_ce(var2dstr,xmlon[ip],xmlat[ip],abs(xmvar[ip]),xaxis,yaxis,x_name+" matched with "+tname,marksize,alphaval,daterange,units,levs[ip],levunits,opt)
          outname = output_path+"MAP_2Dpoints."+str(var2dstrfile)+"_ABS.DRV."+level+str(levunits)+"."+str(dateIN)+".x_"+str(x_name)+".y_"+str(tname)+avgthin_str+".png"
          plt.savefig(outname)
          del outname
         # ... DEPENDENT
        opt=0
        map_points2d_ce(var2dstr,ymlon[ip],ymlat[ip],ymvar[ip],xaxis,yaxis,tname,marksize,alphaval,daterange,units,levs[ip],levunits,opt)
        outname = output_path+"MAP_2Dpoints."+str(var2dstrfile)+".DEP."+level+str(levunits)+"."+str(dateIN)+".x_"+str(x_name)+".y_"+str(tname)+avgthin_str+".png"
        plt.savefig(outname)
        del outname
        if var2dstr.find("HLOS")!=-1:
  	      # plot absolute value (ABS) of HLOS wind velocity
          map_points2d_ce(var2dstr,ymlon[ip],ymlat[ip],abs(ymvar[ip]),xaxis,yaxis,tname,marksize,alphaval,daterange,units,levs[ip],levunits,opt)
          outname = output_path+"MAP_2Dpoints."+str(var2dstrfile)+"_ABS.DEP."+level+str(levunits)+"."+str(dateIN)+".x_"+str(x_name)+".y_"+str(tname)+avgthin_str+".png"
          plt.savefig(outname)
          del outname
        # ... Difference: DEPENDENT - DRIVER
        opt=1
        var2dstr += " Diff"
        map_points2d_ce(var2dstr,xmlon[ip],xmlat[ip],diffvar[ip],xaxis,yaxis,tname+"-"+x_name,marksize,alphaval,daterange,units,levs[ip],levunits,opt)
        outname = output_path+"MAP_2Dpoints."+str(var2dstrfile)+"_Diff.DEP-DRV."+level+str(levunits)+"."+str(dateIN)+".x_"+str(x_name)+".y_"+str(tname)+avgthin_str+".png"
        plt.savefig(outname)
        del outname
        del level

    del xaxis,yaxis,levunits,levs
    del xmlon,xmlat,xmvar,ymlon,ymlat,ymvar,diffvar

    plt.close("all")

    #````````````````````````````
    # DENSITY SCATTERPLOTS
	
	#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	# Global
	#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	
    regionstr = "Global"
	
		#WIND SPEED
    units = "m/s"
    #density_scatter(nDuniq_list[j],txlat,txlon,txspd, tyspd, x_name, tname, units,regionstr)
    density_scatter(txspd, tyspd, x_name, tname, units,regionstr)
    outname = output_path+"DENSITY_SCATTER.Wind."+regionstr+"."+str(dateIN)+".x_"+str(x_name)+".y_"+str(tname)+avgthin_str+".png"
    plt.savefig(outname)
    del outname

    		#PRESSURE
    if tVert_str=="Pressure":
      units = "hPa"
      #density_scatter(nDuniq_list[j],txlat,txlon,txzp, tyzp, x_name, tname, units,regionstr)
      density_scatter(txzp, tyzp, x_name, tname, units,regionstr)
      outname = output_path+"DENSITY_SCATTER.Pressure."+regionstr+"."+str(dateIN)+".x_"+str(x_name)+".y_"+str(tname)+avgthin_str+".png"
      plt.savefig(outname)
      del outname
    else:
    		#HEIGHT
      units = "km"
      #density_scatter(nDuniq_list[j],txlat,txlon,txzp, tyzp, x_name, tname, units,regionstr)
      density_scatter(txzp, tyzp, x_name, tname, units,regionstr)
      outname = output_path+"DENSITY_SCATTER.Height."+regionstr+"."+str(dateIN)+".x_"+str(x_name)+".y_"+str(tname)+avgthin_str+".png"
      plt.savefig(outname)
      del outname

	#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	# Regional (NH,TR,SH)
	#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    for ir in range(nregions):
      Ridx = np.where((txlat>=regions_latmin[ir])*(txlat<regions_latmax[ir])*(tylat>=regions_latmin[ir])*(tylat<regions_latmax[ir]))
      
      regionstr = regions[ir]
      if regionstr=="NH": tnDuniq = nDuniq_listNH#[j]
      if regionstr=="TR": tnDuniq = nDuniq_listTR#[j]
      if regionstr=="SH": tnDuniq = nDuniq_listSH#[j]
	
		#WIND SPEED
      units = "m/s"
      #density_scatter(tnDuniq,txlat[Ridx],txlon[Ridx],txspd[Ridx], tyspd[Ridx], x_name, tname, units,regionstr)
      density_scatter(txspd[Ridx], tyspd[Ridx], x_name, tname, units,regionstr)
      outname = output_path+"DENSITY_SCATTER.Wind."+regionstr+"."+str(dateIN)+".x_"+str(x_name)+".y_"+str(tname)+avgthin_str+".png"
      plt.savefig(outname)
      del outname
    		#PRESSURE
      if tVert_str=="Pressure":
        units = "hPa"
        #density_scatter(tnDuniq,txlat[Ridx],txlon[Ridx],txzp[Ridx], tyzp[Ridx], x_name, tname, units,regionstr)
        density_scatter(txzp[Ridx], tyzp[Ridx], x_name, tname, units,regionstr)
        outname = output_path+"DENSITY_SCATTER.Pressure."+regionstr+"."+str(dateIN)+".x_"+str(x_name)+".y_"+str(tname)+avgthin_str+".png"
        plt.savefig(outname)
        del outname
      elif tVert_str=="Height":
    		#HEIGHT
        units = "km"
        #density_scatter(tnDuniq,txlat[Ridx],txlon[Ridx],txzp[Ridx], tyzp[Ridx], x_name, tname, units,regionstr)
        density_scatter(txzp[Ridx], tyzp[Ridx], x_name, tname, units,regionstr)
        outname = output_path+"DENSITY_SCATTER.Height."+regionstr+"."+str(dateIN)+".x_"+str(x_name)+".y_"+str(tname)+avgthin_str+".png"
        plt.savefig(outname)
        del outname

      del Ridx,regionstr,tnDuniq

    #````````````````````````````
    # Plot 3D MAP of DEPENDENT DATASET wind speed profiles
   
	# All Obs 
		# Wind Velocity
		# ... DRIVER    
    tspd         = txspd
    dset_plotted = "DRV"
    tmp_varname  = sslabelfile
    dset_plotted_name = "DRIVER "+x_name+" matched with "+tname
    if np.size(tspd)>0:
      map_3d_profile(nDuniq_list[j],tspd,txlon,txlat,txvert,x_name,tname,dset_plotted,dset_plotted_name,dateIN,sslabel,zzlabel)
      outname = output_path+"MAP_3D."+str(tmp_varname)+"."+str(dset_plotted)+"."+str(dateIN)+".x_"+str(x_name)+".y_"+str(tname)+avgthin_str+".png"
      plt.savefig(outname)
      del outname
      if tmp_varname.find("HLOS")!=-1:
			# plot absolute value (ABS) of HLOS wind velocity
        map_3d_profile(nDuniq_list[j],abs(tspd),txlon,txlat,txvert,x_name,tname,dset_plotted,dset_plotted_name,dateIN,sslabel,zzlabel)
        outname = output_path+"MAP_3D."+str(tmp_varname)+"_ABS."+str(dset_plotted)+"."+str(dateIN)+".x_"+str(x_name)+".y_"+str(tname)+avgthin_str+".png"
        plt.savefig(outname)
        del outname
    del tspd,dset_plotted,tmp_varname
		# ... DEPENDENT
    tspd         = tyspd
    dset_plotted = "DEP"
    tmp_varname  = sslabelfile
    dset_plotted_name = tname+" matched with DRIVER "+x_name
    if np.size(tspd)>0:
      map_3d_profile(nDuniq_list[j],tspd,tylon,tylat,tyvert,x_name,tname,dset_plotted,dset_plotted_name,dateIN,sslabel,zzlabel)
      outname = output_path+"MAP_3D."+str(tmp_varname)+"."+str(dset_plotted)+"."+str(dateIN)+".x_"+str(x_name)+".y_"+str(tname)+avgthin_str+".png"
      plt.savefig(outname)
      del outname
      if tmp_varname.find("HLOS")!=-1:
        		# plot absolute value (ABS) of HLOS wind velocity
        map_3d_profile(nDuniq_list[j],abs(tspd),tylon,tylat,tyvert,x_name,tname,dset_plotted,dset_plotted_name,dateIN,sslabel,zzlabel)
        outname = output_path+"MAP_3D."+str(tmp_varname)+"_ABS."+str(dset_plotted)+"."+str(dateIN)+".x_"+str(x_name)+".y_"+str(tname)+avgthin_str+".png"
        plt.savefig(outname)
        del outname
    del tspd,dset_plotted,tmp_varname
		# Wind Difference (DEPENDENT - DRIVER)
    tspd         = tyspd - txspd
    dset_plotted = "DEP-DRV"
    tmp_varname  = sslabelfile+"_Diff"
    dset_plotted_name = tname+" minus "+x_name
    if np.size(tspd)>0:
      map_3d_profile(nDuniq_list[j],tspd,txlon,txlat,txvert,x_name,tname,dset_plotted,dset_plotted_name,dateIN,"Diff in "+sslabel,zzlabel)
      outname = output_path+"MAP_3D."+str(tmp_varname)+"."+str(dset_plotted)+"."+str(dateIN)+".x_"+str(x_name)+".y_"+str(tname)+avgthin_str+".png"
      plt.savefig(outname)
      del outname
    del tspd,dset_plotted,tmp_varname

	#````````````````````````````````
	# Strong Winds (> |25| m/s)
    speedthresh  = 25.0
                # Wind Velocity
                # ... DRIVER
    idxD	 = np.where(abs(txspd) > speedthresh)
    tspd         = txspd[idxD]
    dset_plotted = "DRV"
    tmp_varname  = sslabelfile
    dset_plotted_name = "DRIVER "+x_name+" matched with "+tname

    if np.size(tspd)>0:
      map_3d_profile(nDuniq_list[j],tspd,txlon[idxD],txlat[idxD],txvert[idxD],x_name,tname,dset_plotted,dset_plotted_name,dateIN,"Strong "+sslabel,zzlabel)
      outname = output_path+"MAP_3D."+str(tmp_varname)+".Wind_gt_"+str(int(speedthresh))+"."+str(dset_plotted)+"."+str(dateIN)+".x_"+str(x_name)+".y_"+str(tname)+avgthin_str+".png"
      plt.savefig(outname)
      del outname
      if tmp_varname.find("HLOS")!=-1:
		        # plot absolute value (ABS) of HLOS wind velocity
        map_3d_profile(nDuniq_list[j],abs(tspd),txlon[idxD],txlat[idxD],txvert[idxD],x_name,tname,dset_plotted,dset_plotted_name,dateIN,"Strong "+sslabel,zzlabel)
        outname = output_path+"MAP_3D."+str(tmp_varname)+"_ABS.Wind_gt_"+str(int(speedthresh))+"."+str(dset_plotted)+"."+str(dateIN)+".x_"+str(x_name)+".y_"+str(tname)+avgthin_str+".png"
        plt.savefig(outname)
        del outname

    del tspd,dset_plotted,tmp_varname,idxD

                # ... DEPENDENT
    idxDep	 = np.where(abs(tyspd) > speedthresh)
    tspd         = tyspd[idxDep]
    dset_plotted = "DEP"
    tmp_varname  = sslabelfile
    dset_plotted_name = tname+" matched with DRIVER "+x_name

    if np.size(tspd)>0:
      map_3d_profile(nDuniq_list[j],tspd,tylon[idxDep],tylat[idxDep],tyvert[idxDep],x_name,tname,dset_plotted,dset_plotted_name,dateIN,"Strong "+sslabel,zzlabel)
      outname = output_path+"MAP_3D."+str(tmp_varname)+".Wind_gt_"+str(int(speedthresh))+"."+str(dset_plotted)+"."+str(dateIN)+".x_"+str(x_name)+".y_"+str(tname)+avgthin_str+".png"
      plt.savefig(outname)
      del outname
      if tmp_varname.find("HLOS")!=-1:
        		# plot absolute value (ABS) of HLOS wind velocity
        map_3d_profile(nDuniq_list[j],abs(tspd),tylon[idxDep],tylat[idxDep],tyvert[idxDep],x_name,tname,dset_plotted,dset_plotted_name,dateIN,"Strong "+sslabel,zzlabel)
        outname = output_path+"MAP_3D."+str(tmp_varname)+"_ABS.Wind_gt_"+str(int(speedthresh))+"."+str(dset_plotted)+"."+str(dateIN)+".x_"+str(x_name)+".y_"+str(tname)+avgthin_str+".png"
        plt.savefig(outname)
        del outname
    del tspd,dset_plotted,tmp_varname,idxDep

                # Large Wind Differences (DEPENDENT - DRIVER): Diffs > |5| m/s
    tspd         = tyspd - txspd
    diffthresh   = 5.0
    idxDiff 	 = np.where(abs(tspd) > diffthresh)
    dset_plotted = "DEP-DRV"
    tmp_varname  = sslabelfile+"_Diff"
    dset_plotted_name = tname+" minus "+x_name

    if np.size(tspd)>0:
      map_3d_profile(nDuniq_list[j],tspd[idxDiff],txlon[idxDiff],txlat[idxDiff],txvert[idxDiff],x_name,tname,dset_plotted,dset_plotted_name,dateIN,"Large Diff in "+sslabel,zzlabel)
      outname = output_path+"MAP_3D."+str(tmp_varname)+".DIFF_gt_"+str(int(diffthresh))+"."+str(dset_plotted)+"."+str(dateIN)+".x_"+str(x_name)+".y_"+str(tname)+avgthin_str+".png"
      plt.savefig(outname)
      del outname
    del tspd,dset_plotted,tmp_varname,idxDiff

    del speedthresh,diffthresh
    del txspd,tyspd,txlat,txlon,tylat,tylon,tVert_str,sslabel,zzlabel
    del tdate
    
    plt.close("all")

    #````````````````````````````  
    # Rotating Globe Maps
    # 		These plots take the longest to process, which is why they are the last to be plotted

    daterange = dateSTART+"-"+dateEND
    level    = " "
    levunits = " "
    xaxis    = flons
    yaxis    = flats
    xaxisstr = "Longitude"
    yaxisstr = "Latitude"

      # Set latitude around which projection rotates
    center_lat = 0	      # center_lat=0: Projection rotates around Equator

    outname = output_path+"MAP_Rotate.Nobs_Density.DRV.centerlat"+str(center_lat)+"."+str(dateIN)+".x_"+str(x_name)+".y_"+str(tname)+avgthin_str
    contour2d_orthomap_rotate(daterange,output_path,outname,drv_count_rotate,"Collocation Count",xaxis,yaxis,xaxisstr,yaxisstr,x_name+" matched with "+tname,ssunits,center_lat,level,levunits)
    del outname    

    outname = output_path+"MAP_Rotate.Nobs_Density.DEP.centerlat"+str(center_lat)+"."+str(dateIN)+".x_"+str(x_name)+".y_"+str(tname)+avgthin_str
    contour2d_orthomap_rotate(daterange,output_path,outname,dep_count_rotate,"Collocation Count",xaxis,yaxis,xaxisstr,yaxisstr,tname,ssunits,center_lat,level,levunits)
    del outname

    del drv_count_rotate,dep_count_rotate,Dv2d_rotate,vv2d_rotate,diff_rotate,SDdiff_rotate,levs_rotate
    del tname

    plt.close("all")

#---------------------------------------------------------------
#---------------------------------------------------------------

plt.close("all")

#******************************************************************************************************************
#******************************************************************************************************************
# END PLOTS
#******************************************************************************************************************
#******************************************************************************************************************
    
print("========== END MAIN PROGRAM ==========")
print("======================================")
###################################################################################################################
###################################################################################################################
# END MAIN PROGRAM
###################################################################################################################
###################################################################################################################
