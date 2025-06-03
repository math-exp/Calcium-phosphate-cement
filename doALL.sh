#!/usr/local/bin/bash

# parameters in the paper
# Sphi=1
# IPc=(0.0 0.02 0.03), N=200

# select_phi(.)
Sphi=1

IPc=(0.0) # for N=128 TEST
# IPc=(0.0 0.02 0.03) # for N=200

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
