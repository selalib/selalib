# PARAMETERS ARE SET FOR THE FOLLOWING GNUPLOT DIAGNOSTIC
# script to put in sim2d_gc_curvilinear.gnu
set term aqua 1 fsize 20
set logscale y
set format y "%g"
set key bottom
p  'thdiag.dat' u 1:(abs($7-1.)) w l lw 2,\
  '../post_proc/gcsim2d_cartesian_ref.dat' u (abs($1)):($8-1.) \
  w l lw 2 title 'reference',\
  2.e-6*exp(0.2605*x) lw 2
