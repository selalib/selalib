# PARAMETERS ARE SET FOR THE FOLLOWING GNUPLOT DIAGNOSTIC
# script to put in sim2d_gc_curvilinear.gnu
  
set term aqua 1 fsize 20

set logscale y
set format y "%g"

set key bottom right

set title "guiding center 2d polar"

p   'thdiagp.dat' u 1: (sqrt($10)) w l lw 3,\
     '../post_proc/gcsim2d_polar_ref.dat' u 1:($10) w l lw 2 title 'reference',\
     3.5e-6*exp(0.183*x) lw 2
    
