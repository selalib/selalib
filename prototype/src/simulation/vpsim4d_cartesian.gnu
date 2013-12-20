set term aqua 1 fsize 20

set title "Vlasov Poisson 2Dx2D Linear Landau Damping"

set logscale y

p 'thdiag.dat' u ($1):(sqrt($2)) w lp lw 2, \
  '../post_proc/vpsim4d_cartesian_ref.dat' \
  u ($1):(sqrt($2)) w l lw 2 title 'reference', \
  7.e-3*exp(-0.394*x) lw 2