!************************************************
! HOW TO COMPILE USE CAID WITH PYTHON ON 
!  POINCARE MACHINE
!************************************************
source ~/bin/config_gnu_django.sh;module rm epd;module load python

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$PIGASUS_DIR/lib/python2.7/site-packages/pigasus:$PIGASUS_DIR/lib


!---> Important files
!-----------------------------
$HOME/caid/caid/cad_geometry.py
$HOME/caid/src/patchActions.py
$HOME/caid/src/geometryActions.py
