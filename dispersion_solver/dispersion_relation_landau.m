clear all

k = 0.5;
vTx = 1.0;

lambda=@(w) w/(2*k*sqrt(2)*vTx);

D1=@(w) 1+ 1/(k*vTx)^2*(1.0+ lambda(w) *  zetaf( lambda(w) ) );

ztry=2.5+1i*0.08;
options = optimoptions('fsolve','Display','none');
w1=fsolve(D1,ztry,options)
T1=2*pi/real(w1)
%D1(w1)
