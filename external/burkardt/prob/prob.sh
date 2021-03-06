#!/bin/bash
#
mkdir temp
cd temp
rm *
~/bin/$ARCH/f90split ../prob.F90
#
for FILE in `ls -1 *.F90`;
do
  gfortran -c -g $FILE >& compiler.txt
  if [ $? -ne 0 ]; then
    echo "Errors compiling " $FILE
    exit
  fi
  rm compiler.txt
done
rm *.F90
#
ar qc libprob.a *.o
rm *.o
#
mv libprob.a ~/lib/$ARCH
cd ..
rmdir temp
#
echo "Library installed as ~/lib/$ARCH/libprob.a"
