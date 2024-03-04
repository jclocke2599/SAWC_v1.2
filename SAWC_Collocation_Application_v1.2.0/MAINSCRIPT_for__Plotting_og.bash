#!/bin/bash
##############################################################
# Main script to collocate wind datasets with Aeolus winds, and save to user-specified output directory
#
# 	In this script, the user chooses and assigns the following criteria for collocation: paths, dates, datasets, collocation limits, QC, HPC partition
#
# 	To run this script on the command line: ./NameOfThisScript.bash
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
# NOTES: 
#	1. This script runs all jobs in the background. If using S4 supercomputer, it is recommended to run everything on the
#	   's4' partition (6-hour runtime limit). The partition is assigned in this script.
#
# HISTORY:
#	2021-11-18	K.E. Lukens	Created
#	2022-06-08	K.E. Lukens	Expanded ingest capability and added comments throughout
#	2023-04-14	K.E. Lukens	Finalized for upload to SAWC archive
#	2023-04-21	K.E. Lukens	Bug fix: Now uses correct input path for archived index files
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
	# Start Date
yyyyS="2022"
mmS="01"
ddS="02"
hhS="00"

	#`````````````````````````````````````
	# End Date
yyyyE="2022"
mmE="01"
ddE="04"
hhE="18"

#----------------------------------------------
# Set full paths
#       !!! NOTE: Always end path names with a slash "/"

	# Set path 'archive_parent': This is the directory where the /wind-datasets and /collocation directories are located
        #       Example: If the full path to /aircraft is /Full/Path/To/wind-datasets/aircraft ,
        #                then
        #                archive_parent = /Full/Path/To/
        #
archive_parent=/home/jlocke/

	# Set home path 'dir_home': This is the current working directory (i.e., where this script is located)
dir_home="/home/jlocke/SAWC_v1.2/SAWC_Collocation_Application_v1.2.0/"
echo 'WORKING DIRECTORY = '${dir_home}

	# Set input path 'dir_in': This is where the collocation index files are located
	#	Note: 
	#		If index files are located within $archive_parent:
	#			dir_in=${archive_parent}/collocation/index_files/
	#		Otherwise, dir_in should be set to the path where the index files are located
	#
dir_in="/home/jlocke/SAWC_v1.2/index_file/"
echo 'INPUT DIRECTORY = '${dir_in}

	# Set index file suffixes
		# Set datasets to compare and plot
		#	Current Options: Aeolus, Aircraft, AMV_NCEP, Loon, Sonde
		#		If Aeolus is driver, user must add Aeolus wind type: e.g., Aeolus_RayClear_ReprocB11
		#	dsets = sring array containing dependent dataset names
		#	dset_qc = string array indicating whether QC was applied to each dataset during collocation. 
		#		Options: 'QC' if qc was applied; 'NoQC' if qc was not applied
dep_dset=("Aircraft" "AMV_NCEP" "Sonde")

dep_dset_qc=("NoQC" "NoQC" "NoQC")

		# This is the portion of the filename after the DATE. Include all punctuation.
		# This should be the same for ALL index files.
drv_dset='Aeolus_MieCloud_Orig'
#drv_dset='Aeolus_RayClear_ReprocB11'

drv_dset_qc='NoQC'

			# !!! USER SHOULD NOT MODIFY THIS PARAMETER
filename_segment='.drv_'${drv_dset}'__'${drv_dset_qc}'.dset1_'	

	# Set output path 'dir_out': This is where the plots are saved
dir_out="/home/jlocke/SAWC_v1.2/plots/"
echo 'OUTPUT DIRECTORY = '${dir_out}

#----------------------------------------------
# Select choice to super-ob (average) multiple collocations per DRIVER observation, or plot all matches
#      -1 = plot all matches (i.e., do nothing)
#	0 = super-ob (average all matches)

#superob_choice=-1
superob_choice=0

#----------------------------------------------
# Set HPC account, partition, and runtime limit for jobs
#	All variables are strings (use double quotes)

account="star"

qos="normal"

partition="s4"

timelimit="03:00:00"

#----------------------------------------------
# Load modules
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
#!!!!! USERS SHOULD NOT CHANGE ANYTHING BELOW THIS LINE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#----------------------------------------------
# Create output directory (where figures go) if it doesn't already exist

if [[ ! -d ${dir_out}/ ]] ; then
  mkdir -p ${dir_out}/
fi

#----------------------------------------------
# Find sizes of arrays

ndd=${#ddarr[@]}

	# number of dependent datasets to compare to driver
ndsets=${#dep_dset[@]}
echo 'ndsets = '$ndsets

#----------------------------------------------
# Find if Aeolus is a dataset
TFaeolus=False
		# check DRIVER
if echo "$drv_dset" | grep -q "Aeolus"; then #${drv_dset} == "Aeolus" ]; then
  TFaeolus=True
fi
		# check DEPENDENTS
itmp=0
while [[ itmp -lt $ndsets ]]
do
  if echo ${dep_dset[$itmp]} | grep -q "Aeolus"; then
    TFaeolus=True
  fi
  let itmp=itmp+1
done

#----------------------------------------------
# Find if AMV is a dataset
TFamv=False
                # check DRIVER
if echo "$drv_dset" | grep -q "AMV"; then #${drv_dset} == "Aeolus" ]; then
  TFamv=True
fi
if echo "$drv_dset" | grep -q "amv"; then #${drv_dset} == "Aeolus" ]; then
  TFamv=True
fi
                # check DEPENDENTS
itmp=0
while [[ itmp -lt $ndsets ]]
do
  if echo ${dep_dset[$itmp]} | grep -q "AMV"; then
    TFamv=True
  fi
  if echo ${dep_dset[$itmp]} | grep -q "amv"; then
    TFamv=True
  fi
  let itmp=itmp+1
done

#----------------------------------------------
# Set suffix for each input index file to be called for analysis and plotting

dset_in=\"
files_in=\"			#begin string of suffixes with '\"'

itmp=0
while [[ itmp -lt $ndsets ]]
do
  dset_in+=${dep_dset[$itmp]}
  tfile=${filename_segment}${dep_dset[$itmp]}"__"${dep_dset_qc[$itmp]}".nc4"
  files_in+=${tfile}

  let ndsetscheck=$ndsets-1
  if [ $itmp -lt $ndsetscheck ]; then
    dset_in+=","
    files_in+=","
  fi

  let itmp=itmp+1
done
dset_in+=\"
files_in+=\"			#add '\"' to end of string of suffixes
echo "dset_in = "$dset_in
#echo "files_in = "$files_in

#----------------------------------------------
# Set working paths
#	Always end path names with a slash "/"

dir_plot=${dir_home}src/           #location of collocation source code
echo 'PLOTTING CODE DIRECTORY = '${dir_plot}

#----------------------------------------------
# Set name of main collocation code script (to pass as argument to run job script)

plotting_script=MAIN.plotting_code.py
plotting_script_AMV=MAIN.plotting_code.AMVtypes.py

#----------------------------------------------
# Set name of run job script (to pass as argument to run job script)

run_job_script=run_plotting_code.job

#==============================================
#==============================================
# LOOP thru days of selected month
#==============================================

echo "-- Start Date: $yyyyS $mmS $ddS $hhS"
echo "-- End Date:   $yyyyE $mmE $ddE $hhE"
dateSTART=$yyyyS$mmS$ddS$hhS
dateEND=$yyyyE$mmE$ddE$hhE

cd ${dir_home}

#----------------------------------------------
# Run Collocation (matching) algorithm using SLURM (sbatch command)
#	The following initiates collocations with Aeolus Rayleigh-clear and Mie-cloudy winds 
#	for all four 6-hour periods in a single day (+/- 3 hours around each center analysis 
#	time (00, 06, 12, and 18 UTC)). Each combination is run as a separate job, and is 
#	initiated at the same time; thus, 8 jobs will run at once (4 for each Aeolus wind type).
#
#	Run command format for each job:
#		sbatch $output_logfile_name $job_name $input_directory $colloc_script $arg1 /
#		$arg2 $arg3 $arg4 $arg5 $arg6 $arg7 $arg8 $partition $timelimit $run_job_script
#
#       Arguments (arg) to customize collocation:
#               1. Aeolus (driver) wind type: RayClear, MieCloud
#               2. Start Date in YYYYMMDDHH
#               3. End Date in YYYYMMDDHH
#               4. Input directory for index files
#		5. Input dependent dataset names
#		6. Input index file suffix
#               7. Output directory
#		8. Archive parent path: path where home archive directory is located
#		9. Choice to super-ob, thin, or plot all matches

#arg1: automatically set in loop below
arg2=${dateSTART}
arg3=${dateEND}
arg4=${dir_in}
arg5=${dset_in}
arg6=${files_in}
arg7=${dir_out}
arg8=${archive_parent}
arg9=${superob_choice}

if [[ $TFaeolus == True ]]; then		# check dataset_names (all of them) if 'Aeolus' is present
  if [[ "$filename_segment" == *"_Ray"* ]]; then	# check if DRIVER obs (Aeolus) to be plotted are categorized as Rayleigh-clear (Ray)
	# 'Aeolus' is a listed dataset. Submit one job per 6-h analysis cycle (00, 06, 12, 18 UTC) per Aeolus wind type (Rayleigh-clear or Mie-cloudy) = 8 total jobs.

	#`````````````````````````````````````````
        # Collocate with Aeolus Rayleigh-clear

    arg1=RayClear

		#job logfile names
    Rout=LOG_PLOTray_${dateSTART}_${dateEND}
    jRout=PLOTray_${dateSTART}_${dateEND}

    rm $Rout

    sbatch --output=${Rout} --job-name=${jRout} --account=${account} --partition=${partition} --qos=${qos} --time=${timelimit} --export=INDIR=${dir_plot},SCRIPT=${plotting_script},ARG1=${arg1},ARG2=${arg2},ARG3=${arg3},ARG4=${arg4},ARG5=${arg5},ARG6=${arg6},ARG7=${arg7},ARG8=${arg8},ARG9=${arg9} ${run_job_script}

    if [[ $TFamv == True ]]; then
      RoutAMV=LOG_PLOTray_AMV_${dateSTART}_${dateEND}
      jRoutAMV=PLOTray_AMV_${dateSTART}_${dateEND}

      rm $RoutAMV

      sbatch --output=${RoutAMV} --job-name=${jRoutAMV} --account=${account} --partition=${partition} --qos=${qos} --time=${timelimit} --export=INDIR=${dir_plot},SCRIPT=${plotting_script_AMV},ARG1=${arg1},ARG2=${arg2},ARG3=${arg3},ARG4=${arg4},ARG5=${arg5},ARG6=${arg6},ARG7=${arg7},ARG8=${arg8},ARG9=${arg9} ${run_job_script}
    fi    

  elif [[ "$filename_segment" == *"_Mie"* ]]; then    # check if DRIVER obs (Aeolus) to be plotted are categorized as Mie-cloudy (Mie)

	#`````````````````````````````````````````
	# Collocate with Aeolus Mie-cloud

    arg1=MieCloud

		#job logfile names
    Mout=LOG_PLOTmie_${dateSTART}_${dateEND}
    jMout=PLOTmie_${dateSTART}_${dateEND}

    rm $Mout

    sbatch --output=${Mout} --job-name=${jMout} --account=${account} --partition=${partition} --qos=${qos} --time=${timelimit} --export=INDIR=${dir_plot},SCRIPT=${plotting_script},ARG1=${arg1},ARG2=${arg2},ARG3=${arg3},ARG4=${arg4},ARG5=${arg5},ARG6=${arg6},ARG7=${arg7},ARG8=${arg8},ARG9=${arg9} ${run_job_script}

    if [[ $TFamv == True ]]; then
      MoutAMV=LOG_PLOTmie_AMV_${dateSTART}_${dateEND}
      jMoutAMV=PLOTmie_AMV_${dateSTART}_${dateEND}

      rm $MoutAMV

      sbatch --output=${MoutAMV} --job-name=${jMoutAMV} --account=${account} --partition=${partition} --qos=${qos} --time=${timelimit} --export=INDIR=${dir_plot},SCRIPT=${plotting_script_AMV},ARG1=${arg1},ARG2=${arg2},ARG3=${arg3},ARG4=${arg4},ARG5=${arg5},ARG6=${arg6},ARG7=${arg7},ARG8=${arg8},ARG9=${arg9} ${run_job_script}
    fi

  fi

elif [[ $TFaeolus == False ]]; then                # check dataset_names (all of them) if 'Aeolus' is present
	# 'Aeolus' is NOT a listed dataset. Submit one job per 6-h analysis cycle (00, 06, 12, 18 UTC) = 4 total jobs.

  arg1="NA"	# not Aeolus

  out=LOG_PLOT_${dateSTART}_${dateEND}
  jout=PLOT_${dateSTART}_${dateEND}

  rm $out

  sbatch --output=${out} --job-name=${jout} --account=${account} --partition=${partition} --qos=${qos} --time=${timelimit} --export=INDIR=${dir_plot},SCRIPT=${plotting_script},ARG1=${arg1},ARG2=${arg2},ARG3=${arg3},ARG4=${arg4},ARG5=${arg5},ARG6=${arg6},ARG7=${arg7},ARG8=${arg8},ARG9=${arg9} ${run_job_script}

  if [[ $TFamv == True ]]; then
    outAMV=LOG_PLOT_${dateSTART}_${dateEND}
    joutAMV=PLOT_${dateSTART}_${dateEND}

    rm $outAMV

    sbatch --output=${outAMV} --job-name=${joutAMV} --account=${account} --partition=${partition} --qos=${qos} --time=${timelimit} --export=INDIR=${dir_plot},SCRIPT=${plotting_script_AMV},ARG1=${arg1},ARG2=${arg2},ARG3=${arg3},ARG4=${arg4},ARG5=${arg5},ARG6=${arg6},ARG7=${arg7},ARG8=${arg8},ARG9=${arg9} ${run_job_script}
  fi

fi

#----------------------------------------------

#==============================================
# END Date LOOP
#==============================================
#==============================================

# END
##############################################################
