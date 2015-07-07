 set term x11 1
set title 'alpha= 0.100    , k= 0.500    '
set xrange [0:15.000    ]
 set logscale y
 unset logscale x
 l_alpha=  0.10000000000000001     
set xlabel 'Time t in 150 steps, stepwidth= 0.100    '
 set ylabel 'Energy'
 E(x)=0.1*exp(-0.1533*2*x)
 plot 'pic1dresult.dat' using 1:3 with lines,\
 E(x) with lines linestyle 2;
 set autoscale x; unset logscale;unset xlabel;unset ylabel;unset xrange;unset yrange;
 set term x11 2
 set multiplot layout 2,1 rowsfirst title 'Distribution Kernel Density Estimates'
 set autoscale x; set autoscale y
 vel_pdf(x)=1/sqrt(2*pi)*exp(-0.5*(x**2))
 set title 'Velocity Distribution'
 plot './initial_phasespace.dat' using 2:3 smooth kdensity, vel_pdf(x)
 set title 'Spacial Distribution'
 plot './initial_phasespace.dat' using 1:3 smooth kdensity
 unset multiplot
 set term x11 3
 set autoscale x; set autoscale y
 set title 'Thermal velocity estimate'
 set xlabel 'time'
 set ylabel 'v_th'
 plot 'pic1dresult.dat' using 1:5 with lines
 set term x11 4
 set autoscale x; set autoscale y; unset logscale y
 set title 'Absolute Energies'
 set xlabel 'time'
 set ylabel 'energy'
 plot 'pic1dresult.dat' using 1:($2+$3) with lines\
 title 'total (kinetic+electrostatic)',\
 'pic1dresult.dat' using 1:2 with lines \
 title 'kinetic'
 set term x11 5
 set autoscale x; set logscale y
 set title 'Relative Energy Error'
 set xlabel 'time'
 set ylabel 'rel. error'
 plot 'pic1dresult-errors.dat' using 1:2 with lines
 set term x11 6
 set autoscale x; set logscale y
 set title 'Total Impulse Error'
 set xlabel 'time'
 set ylabel 'total. error'
 plot 'pic1dresult-errors.dat' using 1:3 with lines
 set term x11 8
 set multiplot layout 2,1 rowsfirst title 'Time Development of local Pusher Error in x'
 set autoscale x; set autoscale y
 set title 'Mean'
 plot 'pic1dresult.dat' using 1:8 with lines
 set title 'Variance'
 plot 'pic1dresult.dat' using 1:9 with lines
 unset multiplot
