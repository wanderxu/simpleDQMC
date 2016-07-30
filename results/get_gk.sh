#!/bin/bash
source cal_para.sh
WORKDIR="$PWD"
cd $WORKDIR
anaexe=${WORKDIR}/../analysis/cov_tau.out
datadir=${WORKDIR}/../run/

pi=3.141592653

for rhub in ${rhub_array}; do
    cd $datadir
    subdir=U${rhub}L${L}b${beta}
    if [ -f $subdir/gtau.bin ]; then
        cd $subdir

        cp gtau.bin intau
        $anaexe
        rm -f gkbhalf_map.dat

        for((iy=1;iy<=$L;iy++)); do
            ky=` echo "scale=9; (($iy-1)*2*$pi/$L ) - $pi" | bc `
            for((ix=1;ix<=$L;ix++)); do
                ik=` echo "(${iy}-1)*${L} + ${ix}" | bc `
                kx=` echo "scale=9; (($ix-1)*2*$pi/$L ) - $pi" | bc `
                awk '{ if( NR>2 && $1==betav/2 ) print kxv, kyv, $2 }' betav=$beta kxv=$kx kyv=$ky  g_k${ik} >> gkbhalf_map.dat
            done
            kx=$pi
            ik=` echo "(${iy}-1)*${L} + 1" | bc `
            awk '{ if( NR>2 && $1==betav/2 ) print kxv, kyv, $2 }' betav=$beta kxv=$kx kyv=$ky  g_k${ik} >> gkbhalf_map.dat
            echo "" >> gkbhalf_map.dat
        done

        ky=$pi
        iy=1
        for((ix=1;ix<=$L;ix++)); do
            ik=` echo "(${iy}-1)*${L} + ${ix}" | bc `
            kx=` echo "scale=9; (($ix-1)*2*$pi/$L ) - $pi" | bc `
            awk '{ if( NR>2 && $1==betav/2 ) print kxv, kyv, $2 }' betav=$beta kxv=$kx kyv=$ky  g_k${ik} >> gkbhalf_map.dat
        done
        kx=$pi
        ik=` echo "(${iy}-1)*${L} + 1" | bc `
        awk '{ if( NR>2 && $1==betav/2 ) print kxv, kyv, $2 }' betav=$beta kxv=$kx kyv=$ky  g_k${ik} >> gkbhalf_map.dat

cat>plot_gkbhalf_map.gnu<<endin
set terminal postscript eps enhanced color
set output 'gkbhalf_map.eps'
set xtics font 'Times-Roman, 24'
set ytics font 'Times-Roman, 24'
#set key font 'Times-Roman,24' at 5,-0.05
set xlabel "k_x" font 'Times-Roman-Italic,28' offset 0,-1
set ylabel "k_y" font 'Times-Roman-Italic,28' offset -2,0
set xrange[-3.1416:3.1416]
set yrange[-3.1416:3.1416]
set cbrange[0:0.5]
set pm3d map
set size ratio -1
splot "gkbhalf_map.dat" u 1:2:(abs(\$3)) notitle
endin
        gnuplot plot_gkbhalf_map.gnu
        mv gkbhalf_map.eps $WORKDIR/L${L}b${beta}U${rhub}_gkbhalf_map.eps
    fi
done
