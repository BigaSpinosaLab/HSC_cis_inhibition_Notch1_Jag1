#!/bin/bash

#SBATCH --job-name=trimming 
#SBATCH --partition=long
#SBATCH --nodes=1  
#SBATCH --output=logs/trim_cutad_fqmerged.out
#SBATCH --error=logs/trim_cutad_fqmerged.err
#SBATCH --array=1-14%2

#=========================
# User defined parameters: relevant paths
#=========================
# SPECIFY Root directory in the cluster (usually /projects/cancer)
ROOTDIR="/projects/cancer"

# SPECIFY your project working directory. 
# NOTE: It is assumed that Raw data is within 'raw data' folder
#       and Trimmed data is within 'trimmed_data' folder

WKD=$ROOTDIR'/RNAseq_Hematopoietic_AGM_Notch1_RThambyrajah'

#=========================
# General configuration
#=========================
START=$(date +%s)
# Enable Singularity image to look into the general path (equivalent to -B)
export SINGULARITY_BIND=$ROOTDIR 
# Path to images folder in cluster
IMAGES_PATH=$ROOTDIR"/images"
# Path to databases folder in cluster
DB_PATH=$ROOTDIR"/db_files"

################################################################################
##       Reads trimming with cutadapt
################################################################################

# Link to cutadapt manual
# https://cutadapt.readthedocs.io/en/stable/

###########################
## 1. Other relevant paths
###########################

# Folder where raw data is stored
RAWDATA=$WKD'/raw_data_merged'

# Folder where trimmed data will be stored
TRIMDATA=$WKD'/trimmed_merged_data'

#################################################
## 2. Singularity image and Tool Parametrization
#################################################

# Specify image/s name to be used (tool-related)
CUTADAPT='cutadapt_v4.2.sif'  #This image inludes CUTADAPT 4.2

# Specify any particular tool parameters

# Clontech SMARTer II A Oligonucleotide
ADAPTER='AAGCAGTGGTATCAACGCAGAGTAC'

# Maximum error rate (adapter trimming). By default is 10%
ERROR_RATE='0.05'
#No. of allowed errors: 1

################################################################################
## 3. Command file preparation: to execute batch array
################################################################################

cd $RAWDATA

for FILENAME in *.fastq.gz
do
    NAME=${FILENAME%.fastq.gz}
    SAMPLE=$(basename $NAME)
    
    # Construct the full execution command
    echo "singularity exec $IMAGES_PATH/$CUTADAPT cutadapt -g ^$ADAPTER \
                                                          -o $TRIMDATA/$SAMPLE.trimmed.fastq.gz $FILENAME  \
                                                          --cores=0 \
                                                          --no-indels \
                                                          -e $ERROR_RATE 1>> $TRIMDATA/Cutadapt.summary.report.txt"
                                                  
done > $WKD'/scripts/cmds/cutadapt_samples.cmd'

################################################################################
## 4. TrimGalore execution for all available samples
################################################################################

DATE=$(date +%m-%d-%Y--%T)
echo "Starting Trimming in array mode: $DATE"
echo ''

SEEDFILE=$WKD'/scripts/cmds/cutadapt_samples.cmd'
SEED=$(sed -n ${SLURM_ARRAY_TASK_ID}p $SEEDFILE)
eval $SEED

################################################################################
## 4. End
################################################################################

END=$(date +%s)
DIFF=$(( $END - $START ))
echo 'Trimming with cutadapt completed' 
echo "Processing Time: $DIFF seconds"
