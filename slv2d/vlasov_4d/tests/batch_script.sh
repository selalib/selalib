#!/bin/bash
set NPROCS = 16
for CASE_NAME in "landau_cos_prod" "landau_cos_sum" "landau_x" "landau_y" "tsi"
do
   for VA_NUMBER in 0 1 2 3
   do
      cp $CASE_NAME.nml test.nml
      echo "&algo_charge va=$VA_NUMBER /" >> test.nml
      #mpirun -np $NPROCS ../bin/vm4d_spectral_charge ./test.nml
      #mv thf.dat thf-landau-cos-prod-va=0.dat
      touch thf-$CASE_NAME-va=$VA_NUMBER.dat
   done
done
