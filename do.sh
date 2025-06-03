#!/bin/bash

# remove old data
rm mkdatadir.sh
# cd data
# rm *.dat
# cd ..
# get initial curve
# cd init_curve
# make clean
# make
# ./init
# cd ..
# make data directory
cd evolution
make clean
gcc paste-mkdir.c
./a.out
cd ../
chmod +x mkdatadir.sh
./mkdatadir.sh
# get evolution of curves
cd evolution
make clean
make
./paste
# gnuplot
cd ../
chmod +x gnuplot.sh
./gnuplot.sh
# movie
# cp ./data/*/*.mp4 ./movie/
cp ./data/*/*.mpeg ./movie/
cp ./data/*/*.gif ./movie/
