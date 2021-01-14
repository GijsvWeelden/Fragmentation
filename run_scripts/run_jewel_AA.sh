#!/bin/bash
#$ -cwd
#$ -shell yes

export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${HOME}/jewel/lhapdf/lib
export LHAPATH=/cvmfs/alice.cern.ch/x86_64-2.6-gnu-4.1.2/Packages/AliRoot/v5-09-21-1/LHAPDF/PDFsets/

echo $LD_LIBRARY_PATH
cd $PBS_O_WORKDIR
pwd

# Remove old output (needed for dCache file system)
rm -f example.log example.hepmc
rm -f pdfs.dat splitint.dat xsecs.dat 

# /scratch/staff/leeuw179/JEWEL/jewel-2.0.2/jewel-2.0.2-vac params.pp.dat
time $HOME/jewel/jewel-2.2.0/jewel-2.2.0-simple params.AA.dat
