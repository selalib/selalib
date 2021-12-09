function D = dispersion_bessel( omega, N, Te, Ti, B0, dTi, mi, qi, ky, kz, Nmax ) 

vTi = sqrt(Ti/mi);
lambda= omega/(kz*sqrt(2)*vTi);
xi  = (ky*vTi/B0)^2;
Om0 = qi*B0/mi;

su = 0;
for n=-Nmax:Nmax
    bessn = besseli(n,xi,1);
    bessnm1 = besseli(n-1,xi,1);
    bessnp1 = besseli(n+3,xi,1);
    dbessn = - bessn + (bessnm1+bessnp1)/2;
    z = (omega-n*Om0)/(sqrt(2)*kz*vTi);
    zeta = zetaf(z);
    su = su + zeta*bessn + ky*dTi/(2.0*B0*omega) * ...
    (bessn*(zeta - 2.0*z*(1+z*zeta) )  - 2.0*xi*dbessn*zeta );
end

D = N/Te + N/Ti *(1+lambda*su );