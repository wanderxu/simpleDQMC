#!/bin/bash

WORKDIR="$PWD"
source cal_para.sh
echo $WORKDIR
EXE=$WORKDIR/../src/ftdqmc
cd $WORKDIR
for rhub in ${rhub_array}; do
        cd $WORKDIR
        jobdir=U${rhub}L${L}b${beta}
        if [ ! -d $jobdir ]; then
            mkdir $jobdir
        fi

        cd $jobdir

        if [ -f confout ]; then
            cp confout confin
        else
            echo "0" > confin
        fi
cat>ftdqmc.in<<endin
$rhub $L $beta $dtau $nwrap $nsweep $nbin $ltau
 rhub, l, beta, dtau, nwrap, nsweep, nbin, ltau
endin
        echo " Running with parameters U = $rhub, L=$L, beta=$beta "
        $EXE >logs 
done
