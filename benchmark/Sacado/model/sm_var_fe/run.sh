#export OMP_NUM_THREADS=$1
#$TDS/ParAnalRun -i input.hcl -p input.par -n OPN -Sk 1
proj=$PWD
cfle=input
#cfle=input_acm
cd $TDS
tds_rev=$(hg id --num)
cd $proj
if [ -d data_r${tds_rev} ]
then
  rm -rf data_r${tds_rev}
fi
mkdir data_r${tds_rev}
cd data_r${tds_rev}
cp ../${cfle}.hcl ./
echo TDSolve revision r${tds_rev} 2>&1 |tee report
$TDS/TumorDynRun -i ${cfle}.hcl 2>&1 |tee -a report
cd ../
