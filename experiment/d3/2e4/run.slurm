#!/bin/bash
#
#SBATCH --job-name=MLdecoderdz3p2e-4
#SBATCH --output=job_%j.log   # Standard output and error log
#SBATCH --nodes=1          # Run all processes on a single node
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=96
#SBATCH --time=150:00:00
#SBATCH --partition=cpu96

julia mlsurfnc.jl 3 $1 5000 1e-15 10000 $2
