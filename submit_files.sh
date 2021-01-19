
if [ $# -lt 3 ]
then
  echo "Need input directory name: $0 <indir> <startfolder> <endfolder> [nfolder]"
  exit 255
fi

# indir=/home/staff/leeuw179/JEWEL/$1
indir=$1
if [ ! -d $indir ]
then
  echo "Input directory $indir not found"
  exit 1
fi
outdir=/dcache/alice/marcovl/jewel_analysis/`basename $1`
firstfile=$2
lastfile=$3
nfile=${4:-10}
#firstfile=5
#lastfile=20

if [ ! -d $outdir ]
then
  echo "Creating output directory $outdir"
  mkdir -p $outdir
fi

for ((ifile=$firstfile;ifile<=$lastfile;ifile+=$nfile))
do
  qsub -N analyse_files_batch -q long -o log -e log -v indir=$indir,outdir=$outdir,infile=example.hepmc,firstfile=$ifile,nfile=$nfile analyze_files_batch.sh
done

