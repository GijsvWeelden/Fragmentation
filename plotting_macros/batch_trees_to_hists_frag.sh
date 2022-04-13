#!/bin/bash

# Use bash shell on SGE (only relevant on quark)
#$ -S /bin/bash
#
# Use current working directory
#$ -cwd
#
# use log dir for logiles
#$ -o log

#PBS -o $HOME/jewel/analysis/log/analyse_files_batch.o$PBS_O_JOBID
#PBS -e $HOME/jewel/analysis/log/analyse_files_batch.e$PBS_O_JOBID
#PBS -V

root -b ./trees_to_hists_frag.C
