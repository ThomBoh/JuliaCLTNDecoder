#!/bin/bash
#
#SBATCH --job-name=MLdecoderdz9p2e-5
#SBATCH --output=job_%j.log   # Standard output and error log
#SBATCH --nodes=1          # Run all processes on a single node
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --time=150:00:00
#SBATCH --partition=cpu32

julia mlsurfnc.jl 9 $1 46000 1e-15 10000 $2
