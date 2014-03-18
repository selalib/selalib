set term x11 
#set term aqua 1 fsize 20
set logscale y
set key bottom right
p 'thdiag_pade_0p1.dat' u ($1):(sqrt($2)) w l,\
  9.e-7*exp(0.00113*x) lw 2
#  '../post_proc/dksim4d_polar_ref.dat' u ($1):(sqrt($2)) \
#  w l lw 2 title 'reference',\
#  4.e-7*exp(0.00354*x) lw 2