set term aqua 1 fsize 20

set title "Vlasov Poisson 1Dx1D KEEN"

set key top left

#set logscale y

p   'thdiag.dat' u 1:($10) w l lw 3 title 'mode 1',\
  'thdiag.dat' u 1:($11) w l lw 3 title 'mode 2',\
  'thdiag.dat' u 1:($12) w l lw 3 title 'mode 3',\
  'thdiag.dat' u 1:($13) w l lw 3 title 'mode 4',\
  'thdiag.dat' u 1:($14) w l lw 3 title 'mode 5',\
  'thdiag.dat' u 1:($29) w l lw 3 title 'mode 20'

#  '../post_proc/vpsim2d_cartesian_keen_ref.dat' u 1:($10) w l lw 2 title 'mode 1 ref',\
  '../post_proc/vpsim2d_cartesian_keen_ref.dat' u 1:($11) w l lw 2 title 'mode 2 ref',\
  '../post_proc/vpsim2d_cartesian_keen_ref.dat' u 1:($12) w l lw 2 title 'mode 3 ref',\
  '../post_proc/vpsim2d_cartesian_keen_ref.dat' u 1:($13) w l lw 2 title 'mode 4 ref',\
  '../post_proc/vpsim2d_cartesian_keen_ref.dat' u 1:($14) w l lw 2 title 'mode 5 ref'

#f_hat_v diag
#p 'thdiag.dat' u 1:($30) w l lw 3 title 'mode 0',\

#p 'thdiag.dat' u 1:($31) w l lw 3 title 'mode 1',\
'thdiag.dat' u 1:($32) w l lw 3 title 'mode 2',\
'thdiag.dat' u 1:($33) w l lw 3 title 'mode 3',\
'thdiag.dat' u 1:($34) w l lw 3 title 'mode 4',\
'thdiag.dat' u 1:($35) w l lw 3 title 'mode 5',\
'thdiag.dat' u 1:($36) w l lw 3 title 'mode 6',\
'thdiag.dat' u 1:($37) w l lw 3 title 'mode 7',\
'thdiag.dat' u 1:($38) w l lw 3 title 'mode 8',\
'thdiag.dat' u 1:($39) w l lw 3 title 'mode 9',\
'thdiag.dat' u 1:($40) w l lw 3 title 'mode 10',\
'thdiag.dat' u 1:($41) w l lw 3 title 'mode 11',\
'thdiag.dat' u 1:($42) w l lw 3 title 'mode 12',\
'thdiag.dat' u 1:($50) w l lw 3 title 'mode 20'

#p   'thdiag.dat' u 1:($2) w l lw 3 title 'mass'
#p   'thdiag.dat' u 1:($3) w l lw 3 title 'L1norm'
#p   'thdiag.dat' u 1:($4) w l lw 3 title 'momentum'
#p   'thdiag.dat' u 1:($5) w l lw 3 title 'L2norm'
#p   'thdiag.dat' u 1:($6) w l lw 3 title 'kinetic_energy'
#p     'thdiag.dat' u 1:($7) w l lw 3 title 'potential_energy'
#p     'thdiag.dat' u 1:($8) w l lw 3 title 'total_energy'
#p   'thdiag.dat' u 1:($6) w l lw 3 title 'kinetic_energy',\
     'thdiag.dat' u 1:($7) w l lw 3 title 'potential_energy',\
     'thdiag.dat' u 1:($8) w l lw 3 title 'total_energy'
#p 'thdiag.dat' u 1:($9) w l lw 3 title 'mode 0'
