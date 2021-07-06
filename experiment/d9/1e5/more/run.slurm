#!/bin/bash
#
#SBATCH --job-name=MLdecoderdz9p7e-5
#SBATCH --output=job_%j.log   # Standard output and error log
#SBATCH --nodes=1          # Run all processes on a single node
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=96
#SBATCH --time=150:00:00
#SBATCH --partition=cpu96

julia mlsurfnc.jl 9 $1 1000000 1e-15 10000 $2
