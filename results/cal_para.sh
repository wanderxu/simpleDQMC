#!/bin/bash
rhub_array=$(awk 'BEGIN{for(i=1.0;i<=3.01;i+=1.0) printf("%6.2f",i)}')
#rhub_array=$(echo '5.00')
L=4
beta=4.0
dtau=0.05
nwrap=10
nsweep=200
nbin=1
ltau=T

echo " U = " $rhub_array
echo " L = " $L
echo " beta = " $beta
