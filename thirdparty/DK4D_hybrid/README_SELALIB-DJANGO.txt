!************************************************
! HOW TO COMPILE SELALIB AND DJANGO ON POINCARE
!************************************************
source ~/bin/config_gnu_django.sh

1) SELALIB installation
--------------------------------
git clone git+ssh://gforge//gitroot/selalib/selalib.git SELALIB_proto_git
cd SELALIB_proto_git
mkdir build
cd build
ccmake ../src or cmake -DCMAKE_BUILD_TYPE=Release -DBUILD_PACKAGE=ON -DCMAKE_INSTALL_PREFIX=$SELALIB_DIR ../src
make
make install
--> The library will be installed in SELALIB_proto_git/usr. If you want to install it somewhere else specify it with the variable CMAKE_INSTALL_PREFIX

2) SPM installation
--------------------------------
git clone git+ssh://gdgirard@scm.gforge.inria.fr//gitroot/spmmanager/spmmanager.git SPM_git
cd SPM_git
mkdir build
cd build
cmake -DCMAKE_Fortran_COMPILER=gfortran -DCMAKE_BUILD_TYPE=Release -DSPM_DEBUG_TRACE=OFF -DSPM_SAVE_MATRIX=OFF -DCMAKE_INSTALL_PREFIX=$SPM_DIR ..
make && make install

2) DJANGO installation
--------------------------------
git clone git+ssh://gdgirard@scm.gforge.inria.fr//gitroot/jorek/jorek.git DJANGO_git
cd DJANGO_git
git checkout devel-ARA
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=$JOREK_DIR ..
make && make install

3) SELALIB-simulations installation
-----------------------------------------
git clone git+ssh://gdgirard@scm.gforge.inria.fr//gitroot/selalib/sll_simulations.git SELALIB-simulations_git
cd SELALIB-simulations_git
cd VP4D_hybrid
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Release ../src
make
