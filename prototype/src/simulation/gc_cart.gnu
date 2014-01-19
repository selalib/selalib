set term aqua 1 fsize 20

set logscale y
set format y "%g"

set key bottom

set title "guiding center 2d cartesian"

p  'thdiag.dat' u (abs($1)):(abs($8-1.)) w l lw 3,\
  '../post_proc/gcsim2d_cartesian_ref.dat' u (abs($1)):($8-1.) \
  w l lw 2 title 'reference',\
  2.e-6*exp(0.2605*x) lw 2

#p 'thdiag.dat' u 1:5 w l lw 3 title 'mass'
