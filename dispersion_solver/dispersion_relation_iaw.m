options = optimoptions('fsolve','Display','none');
% %Ion acoustic wave from Chen, Chacon, Barnes 2011 JCP 
% me = 1;
% mi = 200;
% Te = 1;
% Ti = 1/10000;
% 
% k = 2*pi/(10*sqrt(Te));
% 
% 
% D1=@(w)1+ 1/(Te*k^2) * (1+w/(k*sqrt(2*Te/me)) * zetaf(w/(k*sqrt(2*Te/me))) )+...
%       1/(Ti*k^2) *(1+w/(k*sqrt(2*Ti/mi)) * zetaf(w/(k*sqrt(2*Ti/mi))) ); 
% ztry= sqrt(me/(mi*(1+1/(Te*k^2))))

% w1=fsolve(D1,ztry,options)
% T1=2*pi/real(w1);
% 
% 
% Ti = 1;
% Te = 10000;
% mi = 1;
% me = 1/200;
% 
% k = 2*pi/(10*sqrt(Te));
% 
% 
% %multispecies
% D2=@(w)1+ 1/(Te*k^2) * (1+w/(k*sqrt(2*Te/me)) * zetaf(w/(k*sqrt(2*Te/me))) )+...
%       1/(Ti*k^2) *(1+w/(k*sqrt(2*Ti/mi)) * zetaf(w/(k*sqrt(2*Ti/mi))) ); 
% ztry= sqrt(me/(mi*(1+1/(Te*k^2))));
% options = optimoptions('fsolve','Display','none');
% w2 = fsolve(D2,ztry,options)
% 
% w1-w2*sqrt(me)
% 
% %adiabatic eletrons
% rt = Te/Ti;
% rm = mi/me;
% 
% D3=@(w) 1+  rt *(1+w/(k*sqrt(2*Ti/mi)) * zetaf(w/(k*sqrt(2*Ti/mi))) ); 
% w3=fsolve(D3,ztry,options)
% T3=2*pi/real(w3);
% %k^2/abs(w2)^2



%Ion acoustic wave from Sturdevant, Parker, Chen, Hause 2016 JCP 
Ti = 1;
mi = 1;
Te = 5;

B0 = 1;
kperp = 0.3;
kpar =  2*pi/1000;
vTi = sqrt(Ti/mi);
xi  = (kperp*vTi/B0)^2;

lambda=@(w) w/(kpar*sqrt(2)*vTi);

ztry = 1.1859 + 0.0*1i;
% Stephan Brunner
D3=@(w) 1.0 + Te/Ti *(1+lambda(w)*zetaf(lambda(w))*besseli(0,xi,1) );
w3=fsolve(D3,ztry,options)
D3(w3)
T3=2*pi/real(w3)

% %Sturdevant
% D4=@(w) 1.0 + Te/Ti *( 1+lambda(w) *zetaf(lambda(w)) ) * besseli(0,xi,1) ;
% w4=fsolve(D4,ztry,options)
% D4(w4)
% T4=2*pi/real(w4)

%me = 200;
%vTe = sqrt(Te/me);
%xe  = (kperp*vTe/B0)^2;
%lambdae=@(w) w/(kpar*sqrt(2)*vTe);
% 
% %Multispecies dispersion
% 
% D4=@(w)  kperp^2 + kpar^2 + 1/Te*(1+lambdae(w)*zetaf(lambdae(w))*besseli(0,xe,1))+...
%    1/Ti*(1+lambda(w)*zetaf(lambda(w))*besseli(0,xi,1));
% w4=fsolve(D4,ztry,options)
% %D3(w3)
% T4=2*pi/real(w4)

