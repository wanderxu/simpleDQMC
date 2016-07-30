#!/bin/bash

WORKDIR="$PWD"
datadir=$WORKDIR/../run/
source cal_para.sh
echo $WORKDIR
cd $WORKDIR

rm L${L}b${beta}_Ekin_vs_U.dat
rm L${L}b${beta}_Docc_vs_U.dat
rm L${L}b${beta}_Spipi_vs_U.dat

for rhub in ${rhub_array}; do
    cd $datadir
    subdir=U${rhub}L${L}b${beta}
    if [ -d $subdir ]; then
        cd $subdir
        echo " Analysing data for parameters U = $rhub, L=$L, beta=$beta "
        awk '{for(i=1;i<=NF;i++) {sum[i] += $i; sumsq[i] += ($i)^2}} \
              END {for (i=1;i<=NF;i++)  { printf( "%12.6f %12.6f \n", sum[i]/NR, sqrt((sumsq[i]-sum[i]^2/NR)/NR) )} }' \
              ener1.bin > tmp.dat
        awk '{if(NR==2) print rhubv, $0 }' rhubv=$rhub tmp.dat >> $WORKDIR/L${L}b${beta}_Ekin_vs_U.dat
        awk '{if(NR==3) print rhubv, $0 }' rhubv=$rhub tmp.dat >> $WORKDIR/L${L}b${beta}_Docc_vs_U.dat
        rm tmp.dat

        awk '{for(i=1;i<=NF;i++) {sum[i] += $i; sumsq[i] += ($i)^2}} \
              END {for (i=1;i<=NF;i++)  { printf( "%12.6f %12.6f \n", sum[i]/NR, sqrt((sumsq[i]-sum[i]^2/NR)/NR) )} }' \
              spin_corrlt.bin |awk '{if(NR==1) print rhubv, $0 }' rhubv=$rhub >> $WORKDIR/L${L}b${beta}_Spipi_vs_U.dat
    fi
done
