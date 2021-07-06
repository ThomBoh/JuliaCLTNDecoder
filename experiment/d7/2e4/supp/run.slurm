#!/bin/bash
#
#SBATCH --job-name=MLdecoderdz7p2e-4
#SBATCH --output=job_%j.log   # Standard output and error log
#SBATCH --nodes=1          # Run all processes on a single node
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --time=150:00:00
#SBATCH --partition=cpu48

julia mlsurfnc.jl 7 $1 4000 1e-15 10000 $2
