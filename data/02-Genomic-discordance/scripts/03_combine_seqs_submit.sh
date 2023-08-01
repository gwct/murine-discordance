#!/bin/bash
#SBATCH --job-name=03_combine_seqs
#SBATCH --output=slurm-logs/%x-%j.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=greggwct@gmail.com
#SBATCH --partition=good_lab_cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=64g
#SBATCH --time=24:00:00

python 03_combine_seqs.py -w 10