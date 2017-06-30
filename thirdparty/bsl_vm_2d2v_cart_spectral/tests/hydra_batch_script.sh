#!/bin/bash
export NPROCS=8
echo "#!/bin/bash" > sub.sh
for CASE_NAME in "landau_cos_prod" "landau_cos_sum" "landau_x" "landau_y" "tsi"
do
   for VA_NUMBER in 0 1 2
   do
      export PREFIX="thf-spectral-$CASE_NAME-va$VA_NUMBER"
      mkdir $PREFIX
      cp $CASE_NAME.nml $PREFIX/input.nml
      cp hydra.ll $PREFIX/batch.ll
      echo "&algo_charge va=$VA_NUMBER /" >> $PREFIX/input.nml
      echo "cd $PREFIX" >> $PREFIX/batch.ll
      echo "poe ../../bin/vm4d_spectral_charge ./input.nml -procs 64" >> $PREFIX/batch.ll
      echo "llsubmit '$PREFIX/batch.ll'" >> sub.sh
   done
done
for VA_NUMBER in 0 2
do
   for METH in 1 2 4
   do 
      for CASE_NAME in "landau_x" "landau_cos_prod" "tsi"
      do
         export PREFIX="thf-csl-$CASE_NAME-va$VA_NUMBER-meth$METH"
         mkdir $PREFIX
         cp $CASE_NAME.nml $PREFIX/input.nml
         cp hydra.ll $PREFIX/batch.ll
         echo "&algo_charge va=$VA_NUMBER meth=$METH /" >> $PREFIX/input.nml
         echo "&solver_poisson mud_case=0/" >> $PREFIX/input.nml
         echo "cd $PREFIX" >> $PREFIX/batch.ll
         echo "poe ../../bin/vm4d_transpose ./input.nml -procs 64" >> $PREFIX/batch.ll
         echo "llsubmit '$PREFIX/batch.ll'" >> sub.sh
      done
   done
done
