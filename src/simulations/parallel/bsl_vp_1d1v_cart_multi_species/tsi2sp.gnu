#set term aqua 1 fsize 20

set title "Vlasov Poisson 1Dx1D TSI 2species"
#Digits:=20:kval:=2*Pi/21:A:=0.01:evalf(A**2*Pi/(2*kval**3));
set key top left

#set logscale y

p 'thdiag.dat' u 1:(sqrt($12)) w l lw 3 title 'electric energy', 'thdiag.ref' u 1:(sqrt($12)) w l lw 3 title 'electric energy',\
