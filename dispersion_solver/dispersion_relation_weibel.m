clear all
%Gempic paper Kraus, Kormann, Morrison, Sonnendrücker 2015 J. Plasma Phys.
k = 1.25;
sigma1 = 0.02/sqrt(2);
sigma2 = sqrt(12)*sigma1;

lambda=@(w) w/(2*k*sqrt(2)*sigma1);

%D1=@(w) w^2 - k^2 -1+ (sigma2/sigma1)^2*(1+ w/(k*sqrt(2)*sigma1) *  plasma_dispersion_function(w/(k*sqrt(2)*sigma1)) ); 
D1=@(w) w^2 - k^2  -1.0 + (sigma2/sigma1)^2*(1.0+ lambda(w) *  zetaf( lambda(w) ) );



D3=@(w) w^2 - 2*k^2  + 1.0*(sigma2/sigma1)^2*(1.0+  lambda(w)  *  zetaf( lambda(w) ) )+ ...
 lambda(w)  *  zetaf( lambda(w) ) -1.0 *(1.0 +lambda(w)  *  zetaf( lambda(w) ) ) ;
 
ztry=0.0-1i*0.03;
options = optimoptions('fsolve','Display','none');
w1=fsolve(D1,ztry,options)
w3=fsolve(D3,ztry,options)
%D1(w1)
%T1=2*pi/real(w1)

% Chen & Chacon 2015 Comput. Phys. Commun. , 2016 JCP
%  me = 1;
%  mi = 25;%100;%1836;
%  Tex = 0.1^2;
%  Tey = 9*Tex; %4*Tex;
%  Tix = Tex;
%  Tiy = Tix; %Tey
%  k = 2; %2*pi/20;
%  
%  %D2 =@(w) 1 - k^2/(w^2) - 1/(me*w^2)*(1+Tey/(2*Tex) * plasma_dispersion_function_derivative(w/(k*sqrt(2*Tex/me))) )-...
%  %    1/(mi*w^2)*(1+Tiy/(2*Tix) * plasma_dispersion_function_derivative(w/(k*sqrt(2*Tix/mi))) ); 
%  
%   
%  D2 =@(w) 1 - k^2/(w^2) + 1/(me*w^2)*(-1+Tey/(Tex) *(1+w/(k*sqrt(2*Tex/me))*zetaf(w/(k*sqrt(2*Tex/me))) ) )+...
%      1/(mi*w^2)*(-1+Tiy/(2*Tix) * (1+ w/(k*sqrt(2*Tix/mi))*zetaf(w/(k*sqrt(2*Tix/mi))) ) ); 
%  
%  ztry=0.0+1i*0.1;
% w2=fsolve(D2,ztry,options);
% D2(w2);
% %T2=2*pi/real(w2)

