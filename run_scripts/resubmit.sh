for ((ifile=1; ifile < 200; ifile++)) 
do   
  if [ ! -s $ifile/example.hepmc ] 
  then 
    echo resubmit $ifile 
    cd $ifile 
    rm -f example.hepmc example.log
    rm -f splitint.dat pdfs.dat xsecs.dat
    qsub run_jewel_AA.sh 
    cd - 
  fi
done
