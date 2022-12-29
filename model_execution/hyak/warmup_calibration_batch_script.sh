#!/bin/bash
#job name
#SBATCH --job-name=epi_model_calibration_HIV_TB_KZN_SA
#SBATCH --account=icrc
#SBATCH --partition=compute

## Resources
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=2:00:00
#SBATCH --mem=6GB

## You can get email notifications when your job starts and stops - useful fr long running jobs
#SBATCH --mail-user=cgreene3@uw.edu
#SBATCH --mail-type=ALL
##SBATCH --output=job.%J.out # tell it to store the output console text to a file called job.<assigned job number>.out
##SBATCH --error=job.%J.err # tell it to store the error messages from the program (if it doesn't write them to normal console output) to a file called job.<assigned job muber>.err

## set array
#######SBATCH --array=1-5000
#SBATCH --array=708,804,2015,2019,2023,2024,2027,2029,2030,2035,2045

module load apptainer

apptainer run --bind /gscratch/icrc/cgreene3/ /gscratch/icrc/cgreene3/warmup_calibration_runs_image.sif $SLURM_ARRAY_TASK_ID

exit 0

