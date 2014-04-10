#!/bin/bash
echo "#!/bin/bash" > sub.sh
for CASE_NAME in "landau_cos_prod" 
do
   for VA_NUMBER in 0 1
   do
      export PREFIX="thf-spectral-$CASE_NAME-va$VA_NUMBER"
      mkdir -p $PREFIX
      cp $CASE_NAME.nml $PREFIX/input.nml
      cp ccamu-intel.oar $PREFIX/batch.oar
      echo "&algo_charge va=$VA_NUMBER /" >> $PREFIX/input.nml
      echo "cd $PREFIX" >> $PREFIX/batch.oar
      echo "mpirun -np 32 ../../bin/vm4d_spectral_charge ./test.nml" >> $PREFIX/batch.oar
      chmod +x $PREFIX/batch.oar
      echo "oarsub './$PREFIX/batch.oar'" >> sub.sh
   done
done
