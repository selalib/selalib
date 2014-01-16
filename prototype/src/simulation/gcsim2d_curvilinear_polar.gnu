# PARAMETERS ARE SET FOR THE FOLLOWING GNUPLOT DIAGNOSTIC
# script to put in sim2d_gc_curvilinear.gnu
  
set term aqua 1 fsize 20

set logscale y
set format y "%g"

set key bottom right

set title "guiding center 2d polar"

p   'thdiagp.dat' u 1: 10 w l lw 3,\
     1.e-11*exp(0.3673*x) lw 2
  #'../post_proc/gcsim2d_polar_ref.dat' u 1:10 w l lw 2 title 'reference',\
