#!/bin/bash
#SBATCH --qos=regular
#SBATCH --time=05:00:00
#SBATCH --nodes=1
#SBATCH --job-name=CMBS4_todcompress
#SBATCH --licenses=SCRATCH
#SBATCH --constraint=cpu
#SBATCH --account=mp107

# set parent directory corresponding to one frequency band
PARENT_DIRECTORY="/global/cfs/cdirs/cmbs4/dc/dc0/staging/noise_sim/outputs_rk/LAT0_CHLAT/f090"

# set range of observations to compress by this script (inclusive)
START_INDEX=1
END_INDEX=30

echo "Listing all observations in $PARENT_DIRECTORY"

# list all observations in parent directory and save to a variable
SUBDIR_LIST=$(find "$PARENT_DIRECTORY" -mindepth 1 -maxdepth 1 -type d | sort)
# extract observations names for printing to console
SUBDIR_NAMES=$(echo "$SUBDIR_LIST" | xargs -I{} basename {})
echo "Observations found: "
echo "$SUBDIR_NAMES"

echo "Proceeding to compress observations indexed in range: $START_INDEX-$END_INDEX"
# select subset of observations (subdirectories) based on range and save to a variable
SELECTED_SUBDIRS=$(echo "$SUBDIR_LIST" | sed -n "${START_INDEX},${END_INDEX}p")
# extract selected observations names for printing to console
SELECTED_SUBDIR_NAMES=$(echo "$SELECTED_SUBDIRS" | xargs -I{} basename {})
echo "Selected observations: "
echo "$SELECTED_SUBDIR_NAMES"

# loop through selected subdirectories and process each one
for subdir in $SELECTED_SUBDIRS; do
  echo "Processing observation: $(basename $subdir)"
  # search for files with the expected starting keywords : 'RISING' or 'SETTING'
  FILE_LIST=$(find "$subdir" -type f \( -name "RISING*" -o -name "SETTING*" \) -printf "%p ")
  # extract file names for printing to console
  echo "Files to compress: "
  for filename in $FILE_LIST; do
    echo $(basename $filename)
  done

  # call compression script on the list of files
  date
  echo "Calling flac compression script ..."
  srun -n 128 python compress_hdf5.py --verify --obs $FILE_LIST > "log_$(basename $subdir).txt" 2>&1

  # if python script runs without error, delete backup files
  if [ $? -eq 0 ]; then
    echo "FLAC compression script ran successfully. Deleting backup files..."
    if find "$subdir" -type f -name "*uncompressed.h5" -delete; then
      echo "Backup files deleted successfully."
      date
    else
      # If backup files deletion fails for some reason we stop the everything to avoid any risk of running out of disk memory.
      echo "Error deleting backup files. Exiting loop over observations."
      date
      break
    fi
  else
    echo "FLAC compression script encountered an error. Not deleting any files."
    date
  fi
done

echo "Observation batch $START_INDEX-$END_INDEX processing in $(basename $PARENT_DIRECTORY) band done. Please verify log files."