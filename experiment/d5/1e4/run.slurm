#!/bin/bash
#
#SBATCH --job-name=MLdecoderdz5p1e-4
#SBATCH --output=job_%j.log   # Standard output and error log
#SBATCH --nodes=1          # Run all processes on a single node
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --time=150:00:00
#SBATCH --partition=cpu32

julia mlsurfnc.jl 5 $1 24000 1e-15 10000 $2
