#set term aqua 1 fsize 20

#set logscale y
set format y "%g"

set xrange [-pi:pi]

#set key bottom

set title "swirling deformation 2d cartesian"


p 'intfdx.dat' u 1:2 w l lw 3,\
#  '../afsim2d_cartesian_ref.dat' w l lw 2 title 'reference'

