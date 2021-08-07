#!/bin/bash
#
#SBATCH --job-name=MLdecoderdz11p7e-5
#SBATCH --output=job_%j.log   # Standard output and error log
#SBATCH --nodes=1          # Run all processes on a single node
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --time=150:00:00
#SBATCH --partition=cpu48

julia mlsurfnc.jl 11 $1 5000 1e-15 10000 $2
