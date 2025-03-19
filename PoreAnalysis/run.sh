#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --time=8:00:00
#SBATCH --mem=30G
#SBATCH --job-name=SWpb
#SBATCH --account=nawimem
#SBATCH --error=%x-%j.err
#SBATCH --output=%x-%j.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=smahajan29@wisc.edu

~/software/PoreBlazer/src/poreblazer_gfortran.exe < input.dat > output.dat
