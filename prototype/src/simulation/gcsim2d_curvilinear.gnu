set term aqua 1 fsize 20

set logscale y
set format y "%g"

set xrange [1:]
set key bottom

set title "guiding center 2d curvilinear"

p 'thdiag.dat' u 1:($7-1) w l lw 3, \
  '../post_proc/gcsim2d_curvilinear_ref.dat' u (abs($1)):($8-1.) \
  w l lw 2 title 'reference',\
  2.e-2*exp(0.2605*x) lw 2
