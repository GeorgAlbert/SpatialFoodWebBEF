#!/bin/bash
## baseline julia script


#SBATCH --job-name=spatialBEF
#SBATCH --chdir=/Julia/
#SBATCH --time=24:00:00
#SBATCH --array=1:5100
#SBATCH -c 1
#SBATCH --mem-per-cpu=8G

# specify output and error files (folde needs to exist!)
#SBATCH --output=/dat/3_outerr/sim_%A-%a.out
#SBATCH --error=/dat/3_outerr/sim_%A-%a.err


# make the module system available and load Julia
module load Julia

# gives some feedback and indicates when it starts (shows in output file)
now=$(date +"%T")
echo "will do my job and call the run script _ $now"

# start script (hands over parameters to be accessible in Julia --> create identifiable output
julia 3_ODEsimulate.jl $SLURM_ARRAY_JOB_ID $SLURM_ARRAY_TASK_ID 

# gives some feedback and indicates when it's done
now=$(date +"%T")
echo "did my job _ $now"
