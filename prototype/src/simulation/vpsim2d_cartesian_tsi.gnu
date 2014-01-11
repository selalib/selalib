set term aqua 1 fsize 20

set title "Vlasov Poisson 1Dx1D Two Stream Instability Int f dx T=35"

#p 'thdiag.dat' u 1:(sqrt($7)) w p

set xrange [0:4]

p 'intfdx.dat' u 1:2 w l lw 3,\
  '../post_proc/vpsim2d_tsi_ref.dat' u 1:2 w l lw 2 title 'reference'

#p 'intfdx_irma2.dat' u 1:2 w l lw 3,\
  'intfdx_irma.dat' u 1:2 w l lw 2

#  '../post_proc/vpsim2d_bot_ref.dat' u 1:(sqrt($7)) w l lw 2 title 'reference'
