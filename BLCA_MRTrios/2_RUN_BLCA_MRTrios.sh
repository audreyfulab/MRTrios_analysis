#!/bin/bash
#SBATCH --job-name=MRTrios_BLCA
#SBATCH --output=/wsu/home/hb/hb68/hb6890/fulab/MRTrios/logs/BLCA_MRTrios_part%a_%j.out
#SBATCH --error=/wsu/home/hb/hb68/hb6890/fulab/MRTrios/logs/BLCA_MRTrios_part%a_%j.err
#SBATCH --partition=prip
#SBATCH --qos=primary
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=150G
#SBATCH --time=7-00:00:00
#SBATCH --array=1-10

# ── Create directories ─────────────────────────────────────
mkdir -p /wsu/home/hb/hb68/hb6890/fulab/MRTrios/logs
mkdir -p /wsu/home/hb/hb68/hb6890/fulab/MRTrios/Output_BLCA

# ── Load modules ───────────────────────────────────────────
module purge
module load gnu9/9.1.0
module load r/4.5.0
module load gsl/2.8


# ── Verify R loaded correctly ──────────────────────────────
echo "R version: $(Rscript --version 2>&1)"
echo "R path   : $(which Rscript)"

# ── Print job info ─────────────────────────────────────────
echo "========================================"
echo "Job ID     : $SLURM_JOB_ID"
echo "Array ID   : $SLURM_ARRAY_TASK_ID"
echo "Node       : $SLURMD_NODENAME"
echo "Memory     : 150G"
echo "Start time : $(date)"
echo "========================================"

# ── Run R script ─────────────────────────────────────
Rscript /wsu/home/hb/hb68/hb6890/fulab/MRTrios/Code/Code_BLCA_MRTrios/BLCA_MRTrios.R \
$SLURM_ARRAY_TASK_ID 10

echo "========================================"
echo "Part $SLURM_ARRAY_TASK_ID finished: $(date)"
echo "========================================"
