#!/bin/bash
#SBATCH  --array=1-2134
#SBATCH  --account=def-aliameli
#SBATCH  --time=5:00:00
#SBATCH  --mem=64000M
#SBATCH  --cpus-per-task=10
#SBATCH  --job-name='HFLM-DS_pred'
#SBATCH  --output=hflmds_pred_logs/sim_%A_%a.out

# Load the standard environment built for Fir's AMD Zen 5 architecture
module load StdEnv/2023 r/4.4.0

# CRITICAL: Prevent OpenBLAS/MKL from spawning implicit threads.
# Allowing native BLAS multi-threading while running a PSOCK parallel backend 
# causes catastrophic core oversubscription and node freezing.
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1

# Ensure the logs directory exists for the array outputs
mkdir -p hflmds_pred_logs

Rscript ./pred_testing.R $SLURM_ARRAY_TASK_ID