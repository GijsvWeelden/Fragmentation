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

# . /scratch/software/ns-sap-env.sh v5-05-56-AN
# . /scratch/software/ns-sap-env.sh vAN-20141001
# . /scratch/software/ns-sap-env.sh -c v5-05-Rev-21-test
# . /scratch/software/ns-sap-env.sh vAN-20150818
# source /cvmfs/alice.cern.ch/etc/login.sh
# eval $(alienv printenv VO_ALICE@AliPhysics::vAN-20180913-1)

lastfile=$firstfile+$nfile
echo "nfile = $nfile"
echo "lastfile = $lastfile"

cd $PBS_O_WORKDIR

if [ ! -d $outdir ]
then
  mkdir $outdir
fi

for ((ifile=$firstfile; ifile<$lastfile; ifile++))
do
  inpathfile=$indir/$ifile/$infile
  if [ -f $inpathfile ]
  then
  # ./analyze_hepmc_jet_frag $indir/$ifile/$infile $outdir/jet_frag_$ifile
  # ./analyze_hepmc_dihadrons $indir/$ifile/$infile $outdir/dihadrons_$ifile
  # ./analyze_hepmc_jet_spectra $indir/$ifile/$infile $outdir/jet_spectra_$ifile
  # ./analyze_hepmc_jet_shapes $indir/$ifile/$infile $outdir/jet_shapes_${ifile}_nobkgsub_full.root
  # ./analyze_hepmc_jet_shapes_constsub $indir/$ifile/$infile $outdir/jet_shapes_${ifile}_constsub_full.root 
  # ./analyze_hepmc_recoil_shapes_constsub $indir/$ifile/$infile $outdir/recoil_shapes_${ifile}_constsub_charged.root 

  # without bkg subtraction
  # ./analyze_hepmc_jet_shapes_constsub_eventwise --chargedjets --nobkg $indir/$ifile/$infile $outdir/jet_shapes_constsub_eventwise_${ifile} 
  ./analyze_hepmc_jet_shapes_constsub_eventwise_treeout --chargedjets --nobkg $indir/$ifile/$infile $outdir/jet_shapes_constsub_eventwise_tree_${ifile} 
  # ./analyze_hepmc_hjet --gammajet --fulljets --nobkg $inpathfile $outdir/gjet_hists_${ifile}
  # ./analyze_hepmc_jet_frag --fulljets --nobkg $indir/$ifile/$infile $outdir/jet_frag_$ifile
  # ./analyze_hepmc_jet_strangeness --fulljets --nobkg $indir/$ifile/$infile $outdir/jet_frag_strangeness_$ifile
  # ./analyze_hepmc_jet_Rdiff $indir/$ifile/$infile $outdir/jet_rdiff_$ifile

  # with bkg subtraction
  # ./analyze_hepmc_jet_shapes_constsub_eventwise --chargedjets $indir/$ifile/$infile $outdir/jet_shapes_constsub_eventwise_${ifile} 
  # ./analyze_hepmc_hjet --chargedjets $inpathfile $outdir/hjet_hists_${ifile}
  #  ./analyze_hepmc_hjet --gammajet --chargedjets $inpathfile $outdir/gjet_hists_${ifile}
  else
    echo "input file $inpathfile not found"
  fi
done
