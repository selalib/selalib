set term aqua 1 fsize 20
set logscale y
set key bottom right
p 'thdiag.dat' u ($1):(sqrt($2)) w l lw 3,\
  '../post_proc/dksim4d_polar_ref.dat' u ($1):(sqrt($2)) \
  w l lw 2 title 'reference',\
  4.e-7*exp(0.00354*x) lw 2