set term aqua 1 fsize 20

set title "Vlasov Poisson 1Dx1D Bump on Tail"

#set logscale y

p 'thdiag.dat' u 1:(sqrt($7)) w p,\
  '../post_proc/vpsim2d_bot_ref.dat' u 1:(sqrt($7)) w l lw 2 title 'reference'
