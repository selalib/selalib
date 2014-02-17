set term aqua 1 fsize 20

set title "Vlasov Poisson 1Dx1D KEEN"

set key top left

#p   'thdiag.dat' u 1:($10) w l lw 3 title 'mode 1',\
  'thdiag.dat' u 1:($11) w l lw 3 title 'mode 2',\
  'thdiag.dat' u 1:($12) w l lw 3 title 'mode 3',\
  'thdiag.dat' u 1:($13) w l lw 3 title 'mode 4',\
  'thdiag.dat' u 1:($14) w l lw 3 title 'mode 5'
p \
  '~/recherche/code_html/to5.dat' u 1:($3) w l lw 3 title 'mode 1 ref',\
  '~/recherche/code_html/to5.dat' u 1:($4) w l lw 3 title 'mode 2 ref',\
  '~/recherche/code_html/to5.dat' u 1:($5) w l lw 3 title 'mode 3 ref',\
  '~/recherche/code_html/to5.dat' u 1:($6) w l lw 3 title 'mode 4 ref',\
  '~/recherche/code_html/to5.dat' u 1:($7) w l lw 3 title 'mode 5 ref',\
  '~/recherche/code_html/to4.dat' u 1:($3) w l lw 2 title 'mode 1',\
  '~/recherche/code_html/to4.dat' u 1:($4) w l lw 2 title 'mode 2',\
  '~/recherche/code_html/to4.dat' u 1:($5) w l lw 2 title 'mode 3',\
  '~/recherche/code_html/to4.dat' u 1:($6) w l lw 2 title 'mode 4',\
  '~/recherche/code_html/to4.dat' u 1:($7) w l lw 2 title 'mode 5'
