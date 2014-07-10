#!/bin/bash
export NPROCS=8
for CASE_NAME in "landau_cos_prod" "landau_cos_sum" "landau_x" "landau_y" "tsi"
do
   for VA_NUMBER in 0 1 2 3
   do
      cp $CASE_NAME.nml test.nml
      echo "&algo_charge va=$VA_NUMBER /" >> test.nml
      mpirun -np $NPROCS ../bin/vm4d_spectral_charge ./test.nml
      mv thf.dat thf-spectral-$CASE_NAME-va=$VA_NUMBER.dat
   done
done
for VA_NUMBER in 0 2
do
   for METH in 1 2 4
   do 
      for TEST_CASE in "landau_x" "landau_cos_prod" "tsi"
      do
         cp $CASE_NAME.nml test.nml
         echo "&algo_charge va=$VA_NUMBER meth=$METH /" >> test.nml
         echo "&solver_poisson mud_case=0/" >> test.nml
         mpirun -np $NPROCS ../bin/vm4d_transpose ./test.nml
         mv thf.dat thf-csl-$CASE_NAME-va=$VA_NUMBER-meth=$METH.dat
      done
   done
done
