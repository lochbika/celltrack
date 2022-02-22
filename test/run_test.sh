#!/bin/bash

set -ex

mkdir -p test_5_timestep_nometa
cd test_5_timestep_nometa
../../build/src/celltrack -i ../data/test_5ts_nc4_zip6.nc -thres 0.001 -var surfprec -nometa -v | tee output.txt
cd ..

mkdir -p test_5_timestep_nometa_buffer
cd test_5_timestep_nometa_buffer
../../build/src/celltrack -i ../data/test_5ts_nc4_zip6.nc -thres 0.001 -var surfprec -buffer 5 -nometa -v | tee output.txt
cd ..

exit 0
