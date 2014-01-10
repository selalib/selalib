set term aqua 1 fsize 20

set title "Vlasov Poisson 1Dx1D Beam no splitting Kinetic energy"


set xrange [0:20]


p 'thdiag.dat' u 1:6 w l lw 3,\
  '../post_proc/vpsim2d_no_split_beam_ref.dat' u 1:6 w l lw 1 title 'reference'


  