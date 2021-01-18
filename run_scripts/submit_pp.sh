if [ $# -lt 2 ]
then
  echo "Specify start, end job nr"
  exit 
fi

start=$1
end=$2

for (( i=$start; i<=$end; i++ ))
do
  mkdir $i
  if [ $? -ne 0 ]
  then
    echo "Directory $i already exists. Skip"
    continue
  fi
  sed 's/\#njob\#/'$i'/' params.pp.dat > $i/params.pp.dat
  cp run_jewel_pp.sh $i
  cd $i
  qsub run_jewel_pp.sh
  cd -
done
