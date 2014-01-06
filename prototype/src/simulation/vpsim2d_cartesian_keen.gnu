set term aqua 1 fsize 20

set title "Vlasov Poisson 1Dx1D KEEN"

set key top left

p   'thdiag.dat' u 1:($10) w l lw 3 title 'mode 1',\
  'thdiag.dat' u 1:($11) w l lw 3 title 'mode 2',\
  'thdiag.dat' u 1:($12) w l lw 3 title 'mode 3',\
  'thdiag.dat' u 1:($13) w l lw 3 title 'mode 4',\
  'thdiag.dat' u 1:($14) w l lw 3 title 'mode 5',\
  '../post_proc/vpsim2d_cartesian_keen_ref.dat' u 1:($10) w l lw 2 title 'mode 1 ref',\
  '../post_proc/vpsim2d_cartesian_keen_ref.dat' u 1:($11) w l lw 2 title 'mode 2 ref',\
  '../post_proc/vpsim2d_cartesian_keen_ref.dat' u 1:($12) w l lw 2 title 'mode 3 ref',\
  '../post_proc/vpsim2d_cartesian_keen_ref.dat' u 1:($13) w l lw 2 title 'mode 4 ref',\
  '../post_proc/vpsim2d_cartesian_keen_ref.dat' u 1:($14) w l lw 2 title 'mode 5 ref'
