#!/bin/bash
#
#SBATCH --job-name=MLdecoderdz3p1e-5
#SBATCH --output=job_%j.log   # Standard output and error log
#SBATCH --nodes=1          # Run all processes on a single node
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --time=150:00:00
#SBATCH --partition=compute

julia mlsurfnc.jl 3 $1 40000 1e-15 10000 $2
