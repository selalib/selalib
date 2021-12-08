R0=239.8081535;
r0 = 0.1;
r1 = 14.5;
rbar = (r0+r1)/2.0;
N   = 0.992378037040461;
kr = 0.055;
Te  = 1.0;
Ti  = 1.0;
B0  = 1.0;
kti = 0.27586;
m = 15;
n = 1;
kpar = n/R0;
mi  = 1.0;
qi  = 1.0;

vTi = sqrt(Ti/mi);

lambda=@(w) w/(kpar*sqrt(2)*vTi);

D=@(w) 1 + Te/Ti *(1+lambda(w)*zetaf(lambda(w)) )  +m*lambda(w)/(rbar*B0*w) * ...
    (zetaf(lambda(w))*(kr-0.5*kti) + kti*lambda(w)*(1+lambda(w)*zetaf(lambda(w))) ) ;

ztry = 0.001 + 0.03*1i;

options = optimoptions('fsolve','Display','none');

w=fsolve(D,ztry,options)
T=2*pi/real(w);


