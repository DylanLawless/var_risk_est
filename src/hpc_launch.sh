#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=62:00:00
#SBATCH --job-name=varriskest
#SBATCH --output=/project/home/lawless/var_risk_est/log/hpc_launch_id_%a_%A_%J.out
#SBATCH --error=/project/home/lawless/var_risk_est/log/hpc_launch_id_%a_%A_%J.err
#SBATCH --array=0-446
# #SBATCH --partition=all-cpu-nodes
#SBATCH --partition=dynamic-64cores-256g-1gpu-A100-40g

set -e
echo "START AT $(date)"
module load r/4.3.1

PanelAppRex_ID=$(sed -n "${SLURM_ARRAY_TASK_ID}p" par_ids.txt)
echo "Now running panel ID: $PanelAppRex_ID"

Rscript inheritance_prob_generalised_hpc.R $PanelAppRex_ID

echo "FINISH AT $(date)"
