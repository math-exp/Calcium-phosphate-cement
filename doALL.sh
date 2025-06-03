#!/usr/local/bin/bash

# select_phi(.)
Sphi=1
IPc=(0.0 0.05 0.1 0.15) # for N=128
# IPc=(0.01 0.02 0.03 0.04) # for N=200
# IPc=(0.0 0.02 0.03) # for N=200
# IPc=(0.0) # for N=128 TEST

for ((i=0;i<${#IPc[@]};i++))
do
    echo "#ifndef HEAD_SPHI_IPC" > head_Sphi_IPc.h
    echo "#define HEAD_SPHI_IPC" >> head_Sphi_IPc.h
    echo " " >> head_Sphi_IPc.h
    echo "#define select_phi ($Sphi)" >> head_Sphi_IPc.h
    echo "#define IPc (${IPc[i]})" >> head_Sphi_IPc.h
    echo "#endif" >> head_Sphi_IPc.h

    ./do.sh
done
