#!/bin/bash
##############################################################
# Main script to collocate wind datasets with Aeolus winds, and save to user-specified output directory
#
# 	In this script, the user chooses and assigns the following criteria for collocation: paths, dates, datasets, collocation limits, QC, HPC partition
#
# 	To run this script on the command line: ./NameOfThisScript.bash
#
# CONTRIBUTORS: Katherine E. Lukens 		NOAA/NESDIS/STAR, CISESS at U. of Maryland
#		Kevin Garrett			NOAA/NWS/OSTI
#		Kayo Ide			U. of Maryland
#		David Santek			CIMSS at U. of Wisconsin-Madison
#		Brett Hoover			NOAA/NWS/NCEP/EMC, Lynker Technologies
#               David Huber                     NOAA/NWS/NCEP/EMC, Redline Performance Solutions, LLC
#		Ross N. Hoffman			NOAA/NESDIS/STAR, CISESS at U. of Maryland
#		Hui Liu				NOAA/NESDIS/STAR, CISESS at U. of Maryland
#
# NOTES: 
#	1. This script runs all jobs in the background. If using S4 supercomputer, it is recommended to run everything on the
#	   's4' partition (6-hour runtime limit). The partition is assigned in this script.
#
# HISTORY:
#	2021-11-18	K.E. Lukens	Created
#	2022-06-08	K.E. Lukens	Expanded ingest capability and added comments throughout
#	2023-04-14	K.E. Lukens	Finalized for upload to SAWC archive
#
##############################################################

#set -x

#==============================================
# BEGIN USER INPUT
#==============================================

#----------------------------------------------
# Set date range over which to run collocation
#	Collocation index files will be generated automatically for the dates specified below.
#	For computational efficiency, it is recommended to loop through one month at a time. 
#	Choose one year/month combination and an array of days.
#
# yyyy 	= year
# mm	= month
# ddarr = day array

	#`````````````````````````````````````
	# Year
yyyy=2022

	#`````````````````````````````````````
	# Month
mm=(01)

	#`````````````````````````````````````
	# Day array - choose ONE, comment out the rest

	# Months with 31 days: Jan, Mar, May, Jul, Aug, Oct, Dec
ddarr=(02 03 04)  #(01 02 03 04 05 06 07 08 09 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31)
	# Months with 30 days: Apr, Jun, Sep, Nov
#ddarr=(01 02 03 04 05 06 07 08 09 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30)
	# Feb - LEAP YEAR
#ddarr=(01 02 03 04 05 06 07 08 09 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29)
	# Feb - NON-LEAP YEAR
#ddarr=(01 02 03 04 05 06 07 08 09 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28)

	# Custom day array
#ddarr=(01 02)

#----------------------------------------------
# Set full paths
#       !!! NOTE: Always end path names with a slash "/"

	# Set home path 'dir_home': This is the current working directory (i.e., where this script is located)
dir_home="/home/jlocke/SAWC_v1.2/SAWC_Collocation_Application_v1.2.0/"
echo 'WORKING DIRECTORY = '${dir_home}

	# Set output path 'dir_out': This is where the output index files are saved
dir_out="/home/jlocke/SAWC_v1.2/index_file/"
echo 'OUTPUT DIRECTORY = '${dir_out}

	# Set path 'archive_parent': This is the directory where the /wind-datasets and /collocation directories are located
	#	Example: If the full path to /aircraft is /Full/Path/To/wind-datasets/aircraft ,
	#	         then 
	#	 	 archive_parent = /Full/Path/To/
	#
archive_parent=/home/jlocke/

#----------------------------------------------
# Set datasets to collocate
#	Users must choose a driver dataset.
#	Users must choose at least one dependent dataset.
#	Users should only choose wind datasets from the following list, as written:
#		Aeolus		... Refers to NetCDF Aeolus winds (both Rayleigh-clear and Mie-cloudy) converted from ESA's Earth Explorer (EE) files.
#		Aircraft	... Refers to NetCDF Aircraft data converted from NCEP aircft prepBUFR files.
#		AMV_4th_Int	... Refers to NetCDF AMVs for 4th International Comparison
#		AMV_NCEP	... Refers to NetCDF AMVs converted from NCEP satwnd prepBUFR files.
#		Loon		... Refers to NetCDF Loon winds.
#		Sonde		... Refers to NetCDF sonde (radio-, drop-, or rawinsonde) winds converted from NCEP adpupa prepBUFR files.

	#set dataset names to use. += appends the names onto the end of the variable 'dataset_names'
	#	NOTES:
	#		1. The first name should be the driver dataset. All subsequent names should be dependent datasets.
	#		2. Make sure to add a comma at the end of each name <-- this is important for the collocation program.
	#		3. BEWARE! While the order of the DEPENDENT dataset names does not matter, the program assumes that 
	#		   the 'qc_flags' order corresponds to the order of 'dataset_names'.
	
dataset_names=("Aeolus" "Aircraft" "AMV_NCEP" "Sonde")
	
AMV_center="none"		# Not to be changed.

#----------------------------------------------
# Define collocation criteria
#       Each criterion should have the same number of elements, and that number should equal the number of DEPENDENT datasets

dst_max=(100.0 100.0 100.0 150.0)          # collocation distance maximum in km ... horizontal distance
prs_max=(0.04 0.04 0.04 0.04)              # colloation pressure difference maximum in log10(hPa) ... for Datasets with given pressures
tim_max=(60.0 60.0 60.0 90.0)              # collocation time difference maximum in minutes
hgt_max=(1.0 1.0 1.0 1.0)                  # collocation height difference maximum in km ... vertical distance ... for Datasets given heights

#----------------------------------------------
# Set quality control (QC) flags for each dataset
#	QC flags should be set for each dataset named above and written as a single string, 
#	with each flag separated by a comma "," <-- this is important for the collocation program.
#
#	Set to 0 if no QC is to be applied
#	Set to 1 if QC is to be applied

qc_flags=(0 0 0 0)

#----------------------------------------------
# Set quality indicator (QI) in percent (%) for AMV dataset(s)
#	QI flags should be set for each AMV dataset named above and written as a single string, 
#	with each flag separated by a comma "," <-- this is important for the collocation program.
#
#	Set to value >= 0 and <= 100
#
#	NOTES: 
#		1. This variable must be set even if AMVs are not used.
#		2. If more than 1 AMV dataset is used, list each QI and separate with comma "," delimiter
#			For example: amv_qi_flags='80,60'
#		   Keep in mind that order matters! First (second) QI listed refers to first (second) AMV dataset listed in 'dataset_names', etc.

#amv_qi_flags=\"0\"
#amv_qi_flags=\"60\"
amv_qi_flags=\"80\"

#----------------------------------------------
# Set quality indicator (QI) choice(s)
#	'qi_choices' is a string
#
#	NO_FC (default value) --> QI without forecast
#	YES_FC                --> QI with forecast
#
#       NOTES:
#               1. This variable must be set even if AMVs are not used.
#               2. If more than 1 AMV dataset is used, list each QI and separate with comma "," delimiter
#                       For example: qi_choices='NO_FC,YES_FC'
#                  Keep in mind that order matters! First (second) QI listed refers to first (second) AMV dataset listed in 'dataset_names', etc.

qi_choices=\"NO_FC\"
#qi_choices=\"YES_FC\"

#----------------------------------------------
# Number of collocations allowed per DRIVER observation

n_max=50

#----------------------------------------------
# Number of processors to use during parallelization

nproc=50

#----------------------------------------------
# Set HPC account, partition, and runtime limit for jobs
#	All variables are strings (use double quotes)

account="star"

qos="normal"

partition="s4"

#timelimit="01:00:00"
timelimit="03:00:00"
#timelimit="06:00:00"

#----------------------------------------------
# For Aeolus dataset only: Choose if using reprocessed or non-reprocessed (original) dataset
#	'bline' and 'dtype' must be specificed.

	#orig = original = Aeolus dataset not reprocessed
bline=orig
dtype=original

	#reprocessed Aeolus dataset: use baseline number as bline (e.g., B10)
#bline=B10
#bline=B11
#bline=B14
#dtype=reprocessed/2${bline}

#----------------------------------------------
# Load Modules
#	This list is taylored for use on the S4 supercomputer.
#	If using another HPC, customize following its own protocol.

module purge

cd

module load license_intel/S4
module load intel/18.0.4  #emc-hpc-stack/2020-q3

module load hdf/
module load hdf5/
module load netcdf4/
module load udunits2
module load ncview
module load nco
module load ncl


#==============================================
# END USER INPUT
#==============================================

#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#!!!!! USER SHOULD NOT CHANGE ANYTHING BELOW THIS LINE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#----------------------------------------------
# Create output directory if it doesn't already exist

if [[ ! -d ${dir_out}/ ]] ; then
  mkdir -p ${dir_out}/
fi

#----------------------------------------------
# Find sizes of arrays

ndd=${#ddarr[@]}

ndsets=${#dataset_names[@]}
echo 'number of datasets = '$ndsets

	#find if Aeolus is a dataset
itmp=0
while [[ itmp -lt $ndsets ]]
do
  if echo ${dataset_names[$itmp]} | grep -q "Aeolus"; then
    TFaeolus=True
    break
  else
    TFaeolus=False
  fi
  let itmp=itmp+1
done

#----------------------------------------------
# Set working paths
#	Always end path names with a slash "/"

dir_coll=${dir_home}src/           #location of collocation source code
echo 'COLLOCATION CODE DIRECTORY = '${dir_coll}

#----------------------------------------------
# Set name of main collocation code script (to pass as argument to run job script)

colloc_script=MAIN.match_driver_dependents.py

#----------------------------------------------
# Set name of run job script (to pass as argument to run job script)

run_job_script=run_collocation_code.job

#==============================================
#==============================================
# LOOP thru days of selected month
#==============================================

idy=0
while [[ idy -lt $ndd ]]
do
  dd=${ddarr[$idy]}		#day

  echo "========== DATE: $yyyy $mm $dd =========="
  yyyymmdd=$yyyy$mm$dd

#----------------------------------------------
# Activate .yml file to run python codes

  cd ${dir_coll}

  #source activate bhoover-Universal_AMV_Matching &		# activates conda environment needed to run collocation code, and run in background (&)
  source activate SAWC_conda_env &		# activates conda environment needed to run collocation code, and run in background (&)

  cd ${dir_home}

#----------------------------------------------
# Run Collocation (matching) Program using SLURM (sbatch command)
#	The following initiates conlocations with Aeolus Rayleigh-clear and Mie-cloudy winds 
#	for all four 6-hour periods in a single day (+/- 3 hours around each center analysis 
#	time (00, 06, 12, and 18 UTC)). Each combination is run as a separate job, and is 
#	initiated at the same time; thus, 8 jobs will run at once (4 for each Aeolus wind type).
#
#	Run command format for each job:
#		sbatch $output_logfile_name $job_name $input_directory $conloc_script $arg1 /
#		$arg2 $arg3 $arg4 $arg5 $arg6 $arg7 $arg8 $arg9 $arg10 $arg11 $arg12 $arg13 /
#               $arg14 $arg15 $arg16 $partition $timelimit $run_job_script
#
#       Arguments (arg) to customize conlocation:
#               1. Aeolus (driver) wind type: RayClear, MieCloud
#               2. QC applied to each dataset? 0=no, 1=yes. Separated by comma "," delimiter
#               3. If using AMVs, this is the lower limit of AMV QI in %, for QC
#               4. If using AMVs, this is the choice of QI type to use, for QC
#               5. Date in YYYYMMDDHH
#               6. Aeolus dataset type: for non-reprocessed data, use orig (original); for reprocessed files, use baseline number (example: B10)
#               7. Output directory
#		8. Colnocation distance maximum
#	        9. Colnocation pressure difference maximum
#	       10. Colnocation time difference maximum
#	       11. Colnocation height difference maximum
#	       12. Names of all datasets to use for colnocation, separated by comma "," delimiter
#	       13. Archive parent path: path where home archive directory is located
#              14. Maximum number of matches allowed per data point
#              15. Number of processors to use during parallelization
#	       16. Abbreviation of wind-producing center

  #arg1: automatically set in loop below
  #arg2=${qc_flags}
  arg3=${amv_qi_flags}
  arg4=${qi_choices}
  #arg5: automatically set in loop below
  arg6=${bline}
  arg7=${dir_out}
  #arg8=${dst_max}
  #arg9=${prs_max}
  #arg10=${tim_max}
  #arg11=${hgt_max}
  #arg12=${dataset_names}
  arg13=${archive_parent}
  arg14=${n_max}
  arg15=${nproc}
  arg16=${AMV_center}


if [[ $TFaeolus == True ]]; then		# check dataset_names (all of them) if 'Aeolus' is present
	# 'Aeolus' is a listed dataset. Submit one job per 6-h analysis cycle (00, 06, 12, 18 UTC) per Aeolus wind type (Rayleigh-clear or Mie-cloudy) = 8 total jobs.

  	# LOOP through datasets
  idset=1
  while [[ idset -lt $ndsets ]]
  do
    driv=${dataset_names[0]}		#driver dataset
    dset=${dataset_names[$idset]}	#dependent dataset

	#set additional job arguments
    arg2=\"${qc_flags[0]},${qc_flags[$idset]}\"
 
    idsetmax=$idset
    let idsetmax=idsetmax-1
    arg8=\"${dst_max[$idsetmax]}\" 
    arg9=\"${prs_max[$idsetmax]}\" 
    arg10=\"${tim_max[$idsetmax]}\" 
    arg11=\"${hgt_max[$idsetmax]}\" 

    arg12=\"${driv},${dset}\"

    #`````````````````````````````````````````
    # Collocate with Aeolus Rayleigh-clear

    arg1=RayClear

		#job names
    Rjob00=COLray_${yyyy}${mm}${dd}00_${driv}${dset}
    Rjob06=COLray_${yyyy}${mm}${dd}06_${driv}${dset}
    Rjob12=COLray_${yyyy}${mm}${dd}12_${driv}${dset}
    Rjob18=COLray_${yyyy}${mm}${dd}18_${driv}${dset}

		#job logfile names
    Rout00=LOG_COLray_${yyyy}${mm}${dd}00_${driv}${dset}
    Rout06=LOG_COLray_${yyyy}${mm}${dd}06_${driv}${dset}
    Rout12=LOG_COLray_${yyyy}${mm}${dd}12_${driv}${dset}
    Rout18=LOG_COLray_${yyyy}${mm}${dd}18_${driv}${dset}

    rm $Rout00 $Rout06 $Rout12 $Rout18

		# hour = 00 UTC
    arg5=${yyyymmdd}00
    sbatch --output=${Rout00} --job-name=${Rjob00} --account=${account} --partition=${partition} --qos=${qos} --time=${timelimit} --export=INDIR=${dir_coll},SCRIPT=${colloc_script},ARG1=${arg1},ARG2=${arg2},ARG3=${arg3},ARG4=${arg4},ARG5=${arg5},ARG6=${arg6},ARG7=${arg7},ARG8=${arg8},ARG9=${arg9},ARG10=${arg10},ARG11=${arg11},ARG12=${arg12},ARG13=${arg13},ARG14=${arg14},ARG15=${arg15},ARG16=${arg16} ${run_job_script}
		# hour = 06 UTC
    arg5=${yyyymmdd}06
    sbatch --output=${Rout06} --job-name=${Rjob06} --account=${account} --partition=${partition} --qos=${qos} --time=${timelimit} --export=INDIR=${dir_coll},SCRIPT=${colloc_script},ARG1=${arg1},ARG2=${arg2},ARG3=${arg3},ARG4=${arg4},ARG5=${arg5},ARG6=${arg6},ARG7=${arg7},ARG8=${arg8},ARG9=${arg9},ARG10=${arg10},ARG11=${arg11},ARG12=${arg12},ARG13=${arg13},ARG14=${arg14},ARG15=${arg15},ARG16=${arg16} ${run_job_script}
		# hour = 12 UTC
    arg5=${yyyymmdd}12
    sbatch --output=${Rout12} --job-name=${Rjob12} --account=${account} --partition=${partition} --qos=${qos} --time=${timelimit} --export=INDIR=${dir_coll},SCRIPT=${colloc_script},ARG1=${arg1},ARG2=${arg2},ARG3=${arg3},ARG4=${arg4},ARG5=${arg5},ARG6=${arg6},ARG7=${arg7},ARG8=${arg8},ARG9=${arg9},ARG10=${arg10},ARG11=${arg11},ARG12=${arg12},ARG13=${arg13},ARG14=${arg14},ARG15=${arg15},ARG16=${arg16} ${run_job_script}
		# hour = 18 UTC
    arg5=${yyyymmdd}18
    sbatch --output=${Rout18} --job-name=${Rjob18} --account=${account} --partition=${partition} --qos=${qos} --time=${timelimit} --export=INDIR=${dir_coll},SCRIPT=${colloc_script},ARG1=${arg1},ARG2=${arg2},ARG3=${arg3},ARG4=${arg4},ARG5=${arg5},ARG6=${arg6},ARG7=${arg7},ARG8=${arg8},ARG9=${arg9},ARG10=${arg10},ARG11=${arg11},ARG12=${arg12},ARG13=${arg13},ARG14=${arg14},ARG15=${arg15},ARG16=${arg16} ${run_job_script}

    #`````````````````````````````````````````
    # Collocate with Aeolus Mie-cloud

    arg1=MieCloud
  
  		#job names
    Mjob00=COLmie_${yyyy}${mm}${dd}00_${driv}${dset}
    Mjob06=COLmie_${yyyy}${mm}${dd}06_${driv}${dset}
    Mjob12=COLmie_${yyyy}${mm}${dd}12_${driv}${dset}
    Mjob18=COLmie_${yyyy}${mm}${dd}18_${driv}${dset}

		#job logfile names
    Mout00=LOG_COLmie_${yyyy}${mm}${dd}00_${driv}${dset}
    Mout06=LOG_COLmie_${yyyy}${mm}${dd}06_${driv}${dset}
    Mout12=LOG_COLmie_${yyyy}${mm}${dd}12_${driv}${dset}
    Mout18=LOG_COLmie_${yyyy}${mm}${dd}18_${driv}${dset}

    rm $Mout00 $Mout06 $Mout12 $Mout18

                # hour = 00 UTC
    arg5=${yyyymmdd}00
    sbatch --output=${Mout00} --job-name=${Mjob00} --account=${account} --partition=${partition} --qos=${qos} --time=${timelimit} --export=INDIR=${dir_coll},SCRIPT=${colloc_script},ARG1=${arg1},ARG2=${arg2},ARG3=${arg3},ARG4=${arg4},ARG5=${arg5},ARG6=${arg6},ARG7=${arg7},ARG8=${arg8},ARG9=${arg9},ARG10=${arg10},ARG11=${arg11},ARG12=${arg12},ARG13=${arg13},ARG14=${arg14},ARG15=${arg15},ARG16=${arg16} ${run_job_script}
                # hour = 06 UTC
    arg5=${yyyymmdd}06
    sbatch --output=${Mout06} --job-name=${Mjob06} --account=${account} --partition=${partition} --qos=${qos} --time=${timelimit} --export=INDIR=${dir_coll},SCRIPT=${colloc_script},ARG1=${arg1},ARG2=${arg2},ARG3=${arg3},ARG4=${arg4},ARG5=${arg5},ARG6=${arg6},ARG7=${arg7},ARG8=${arg8},ARG9=${arg9},ARG10=${arg10},ARG11=${arg11},ARG12=${arg12},ARG13=${arg13},ARG14=${arg14},ARG15=${arg15},ARG16=${arg16} ${run_job_script}
                # hour = 12 UTC
    arg5=${yyyymmdd}12
    sbatch --output=${Mout12} --job-name=${Mjob12} --account=${account} --partition=${partition} --qos=${qos} --time=${timelimit} --export=INDIR=${dir_coll},SCRIPT=${colloc_script},ARG1=${arg1},ARG2=${arg2},ARG3=${arg3},ARG4=${arg4},ARG5=${arg5},ARG6=${arg6},ARG7=${arg7},ARG8=${arg8},ARG9=${arg9},ARG10=${arg10},ARG11=${arg11},ARG12=${arg12},ARG13=${arg13},ARG14=${arg14},ARG15=${arg15},ARG16=${arg16} ${run_job_script}
                # hour = 18 UTC
    arg5=${yyyymmdd}18
    sbatch --output=${Mout18} --job-name=${Mjob18} --account=${account} --partition=${partition} --qos=${qos} --time=${timelimit} --export=INDIR=${dir_coll},SCRIPT=${colloc_script},ARG1=${arg1},ARG2=${arg2},ARG3=${arg3},ARG4=${arg4},ARG5=${arg5},ARG6=${arg6},ARG7=${arg7},ARG8=${arg8},ARG9=${arg9},ARG10=${arg10},ARG11=${arg11},ARG12=${arg12},ARG13=${arg13},ARG14=${arg14},ARG15=${arg15},ARG16=${arg16} ${run_job_script}

    #`````````````````````````````````````````

    let idset=idset+1
  done

elif [[ $TFaeolus == False ]]; then               # check dataset_names (all of them) if 'Aeolus' is present
	# 'Aeolus' is NOT a listed dataset. Submit one job per 6-h analysis cycle (00, 06, 12, 18 UTC) = 4 total jobs.

	# LOOP through datasets
  idset=1
  while [[ idset -lt $ndsets ]]
  do
    driv=${dataset_names[0]}		#driver dataset
    dset=${dataset_names[$idset]}	#dependent dataset

	#set additional job arguments
    arg1="NA"				#not Aeolus

    arg2=\"${qc_flags[0]},${qc_flags[$idset]}\"
 
    idsetmax=$idset
    let idsetmax=idsetmax-1
    arg8=\"${dst_max[$idsetmax]}\" 
    arg9=\"${prs_max[$idsetmax]}\" 
    arg10=\"${tim_max[$idsetmax]}\" 
    arg11=\"${hgt_max[$idsetmax]}\" 

    arg12=\"${driv},${dset}\"

		#job names
    job00=COL_${yyyy}${mm}${dd}00_${driv}${dset}
    job06=COL_${yyyy}${mm}${dd}06_${driv}${dset}
    job12=COL_${yyyy}${mm}${dd}12_${driv}${dset}
    job18=COL_${yyyy}${mm}${dd}18_${driv}${dset}

		#job logfile name
    out00=LOG_COL_${yyyy}${mm}${dd}00_${driv}${dset}
    out06=LOG_COL_${yyyy}${mm}${dd}06_${driv}${dset}
    out12=LOG_COL_${yyyy}${mm}${dd}12_${driv}${dset}
    out18=LOG_COL_${yyyy}${mm}${dd}18_${driv}${dset}

    rm $out00 $out06 $out12 $out18
    
    #`````````````````````````````````````````

  		# hour = 00 UTC
    arg5=${yyyymmdd}00
    sbatch --output=${out00} --job-name=${job00} --account=${account} --partition=${partition} --qos=${qos} --time=${timelimit} --export=INDIR=${dir_coll},SCRIPT=${colloc_script},ARG1=${arg1},ARG2=${arg2},ARG3=${arg3},ARG4=${arg4},ARG5=${arg5},ARG6=${arg6},ARG7=${arg7},ARG8=${arg8},ARG9=${arg9},ARG10=${arg10},ARG11=${arg11},ARG12=${arg12},ARG13=${arg13},ARG14=${arg14},ARG15=${arg15},ARG16=${arg16} ${run_job_script}
                # hour = 06 UTC
    arg5=${yyyymmdd}06
    sbatch --output=${out06} --job-name=${job06} --account=${account} --partition=${partition} --qos=${qos} --time=${timelimit} --export=INDIR=${dir_coll},SCRIPT=${colloc_script},ARG1=${arg1},ARG2=${arg2},ARG3=${arg3},ARG4=${arg4},ARG5=${arg5},ARG6=${arg6},ARG7=${arg7},ARG8=${arg8},ARG9=${arg9},ARG10=${arg10},ARG11=${arg11},ARG12=${arg12},ARG13=${arg13},ARG14=${arg14},ARG15=${arg15},ARG16=${arg16} ${run_job_script}
                # hour = 12 UTC
    arg5=${yyyymmdd}12
    sbatch --output=${out12} --job-name=${job12} --account=${account} --partition=${partition} --qos=${qos} --time=${timelimit} --export=INDIR=${dir_coll},SCRIPT=${colloc_script},ARG1=${arg1},ARG2=${arg2},ARG3=${arg3},ARG4=${arg4},ARG5=${arg5},ARG6=${arg6},ARG7=${arg7},ARG8=${arg8},ARG9=${arg9},ARG10=${arg10},ARG11=${arg11},ARG12=${arg12},ARG13=${arg13},ARG14=${arg14},ARG15=${arg15},ARG16=${arg16} ${run_job_script}
                # hour = 18 UTC
    arg5=${yyyymmdd}18
    sbatch --output=${out18} --job-name=${job18} --account=${account} --partition=${partition} --qos=${qos} --time=${timelimit} --export=INDIR=${dir_coll},SCRIPT=${colloc_script},ARG1=${arg1},ARG2=${arg2},ARG3=${arg3},ARG4=${arg4},ARG5=${arg5},ARG6=${arg6},ARG7=${arg7},ARG8=${arg8},ARG9=${arg9},ARG10=${arg10},ARG11=${arg11},ARG12=${arg12},ARG13=${arg13},ARG14=${arg14},ARG15=${arg15},ARG16=${arg16} ${run_job_script}

    #`````````````````````````````````````````

    let idset=idset+1
  done

fi				#END IF: Aeolus in 'dataset_names'

#----------------------------------------------

  let idy=idy+1
done

#==============================================
# END Date LOOP
#==============================================
#==============================================

# END
##############################################################
