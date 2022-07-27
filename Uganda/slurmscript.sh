#!/bin/bash
#SBATCH --job-name=RCCTraining
#SBATCH --account=stf
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=00:02:00
#SBATCH --mem=4GB
# You can get email notifications when your job starts and stops - useful for long running jobs
##SBATCH --mail-user=<yournetID>@uw.edu
##SBATCH --mail-type=ALL


module load cesg/python

python FizzBuzz.py > output.log

   

exit 0
