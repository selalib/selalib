set term aqua 1 fsize 20

set title "Vlasov Poisson 1Dx1D Linear Landau Damping"

set logscale y

p 'thdiag.dat' u 1:(sqrt($7)) w l lw 3, \
  '../post_proc/vpsim2d_cartesian_ref.dat'  \
  u 1:(sqrt($7)) w l lw 2 title 'reference', \
  2.6e-3*exp(-0.1533*x) lw 2