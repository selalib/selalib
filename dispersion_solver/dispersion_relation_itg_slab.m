N   = 1.0;
Te  = 5.0;%4.0;%1.0;%
Ti  = 1.0;
B0  = 1.0;
dTi = 0.0;%-0.06;
ky  = 0.3;%0.2;%0.8;%
kz  = 2*pi/1000;%0.002;%1/120;%
mi  = 1.0;
qi  = 1.0;
Nmax = 10;

vTi = sqrt(Ti/mi);

D1=@(w) (w^2- Te*kz^2/mi)*w + Te*ky*kz^2*dTi/(B0*mi);

lambda=@(w) w/(kz*sqrt(2)*vTi);

D2=@(w) N/Te + N/Ti *(1+lambda(w)*zetaf(lambda(w)) +ky*dTi/(2.0*B0*w) * ...
    lambda(w)*(zetaf(lambda(w)) + lambda(w)*( -2.0*(1+zetaf(lambda(w))*lambda(w)) )  ));

xi  = (ky*vTi/B0)^2;
bess0 = besseli(0,xi,1);
bess1 = besseli(1,xi,1);
dbess0 = bess1 - bess0; 

D3=@(w) N/Te + N/Ti *(1+lambda(w)*zetaf(lambda(w))*bess0 +ky*dTi/(2.0*B0*w) * ...
    (bess0*lambda(w)*(zetaf(lambda(w)) - lambda(w)*( 2.0*(1+zetaf(lambda(w))*lambda(w)) ) ) -...
    2.0*xi*(bess1-bess0)*lambda(w)*zetaf(lambda(w))  ));


D4=@(omega) dispersion_bessel( omega, N, Te, Ti, B0, dTi, mi, qi, ky, kz, Nmax ) ;


options = optimoptions('fsolve','Display','none');
%ztry = -0.001 + 0.003*1i;

ztry = -0.001 - 0.003*1i;

w1=fsolve(D1,ztry,options)
T1=2*pi/real(w1);

w2=fsolve(D2,ztry,options)
T2=2*pi/real(w2);

w3=fsolve(D3,ztry,options)
T3=2*pi/real(w3);

w4=fsolve(D4,ztry,options)
T4=2*pi/real(w4);

N1 = 10;
N2 = 10;
x0 = 0.001;
dx0r = 0.001;
dx0i = -0.001*1i;
for j=1:N1
    x0 = real(x0)+dx0r;
    for l=1:N2
        x0 = x0 + dx0i;
        [x((j-1)*N2+l), fx((j-1)*N2+l)] = fsolve(D4,x0, options);
    end
end

ind = find(abs(fx)<1e-6);
[ximag,indm]=max(imag(x(ind)))
x(ind(indm))
