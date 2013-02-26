#!/bin/bash
 
HOST=`hostname -s`
echo "${HOST}"
ARCH=`uname -s`
echo "ARCH:${ARCH}"
git log --date-order --date=short | \
sed -e '/^commit.*$/d' | \
awk '/^Author/ {sub(/\\$/,""); getline t; print $0 t; next}; 1' | \
sed -e 's/^Author: //g' | \
sed -e 's/>Date:   \([0-9]*-[0-9]*-[0-9]*\)/>\t\1/g' | \
sed -e 's/^\(.*\) \(\)\t\(.*\)/\3    \1    \2/g' > ChangeLog
mkdir build
cd build; {
cmake ..
make NightlyUpdate
make NightlyConfigure
make NightlyBuild
make NightlyTest || true
if [ -f Testing/TAG ] ; then
   xsltproc ../cmake/ctest2junix.xsl Testing/`head -n 1 < Testing/TAG`/Test.xml > CTestResults.xml
fi
make NightlySubmit
}; cd -
rm -rf build
exit 0
