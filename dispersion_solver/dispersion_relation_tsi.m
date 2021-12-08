%Two stream instability
Te1=1;
Te2=1;
r=Te1/Te2;
k=0.2;
delta=0.5;
v0 = [-2.4, 2.4];

m=200;
T=1000;

% for i=1:2
%     zeta(i) = (omega/k - v(i))/(sqrt(2)*vth(i));
% end
% D = 1 + 1/(k^2)*(delta(1)*(1+zeta(1)*Z(1))/(vth(1)^2) +
% delta(2)*(1+zeta(2)*Z(2))/(vth(2)^2));

%single species electrons
me = 1;
D1=@(w) 1 + 1/k^2*(delta/Te1*(1+(w/k-v0(1))/sqrt(2*Te1/me).*zetaf((w/k-v0(1))/sqrt(2*Te1/me)) ) +...
    (1-delta)/Te2*(1+(w/k-v0(2))/sqrt(2*Te2/me).*zetaf((w/k-v0(2))/sqrt(2*Te2/me)) ) );

ztry=0.1;
%options = optimoptions('fsolve','Display','none','PlotFcn',@optimplotfirstorderopt);
options = optimoptions('fsolve','Display','none');

w1=fsolve(D1,ztry,options)

%multispecies for electrons as reference species
me = 1;
mi = m;
Ti = 1/T;
Te = 1;
D2=@(w)1 + 1/k^2*(delta/Te1*(1+(w/k-v0(1))/sqrt(2*Te1/me).*zetaf((w/k-v0(1))/sqrt(2*Te1/me)) ) +...
    (1-delta)/Te2*(1+(w/k-v0(2))/sqrt(2*Te2/me).*zetaf((w/k-v0(2))/sqrt(2*Te2/me)) ) ) +...
    1/(k^2*Ti)*(1+w/(k*sqrt(2*Ti/mi))*zetaf(w/(k*sqrt(2*Ti/mi))) ); 
ztry = 0.;
w2=fsolve(D2,ztry,options)


%multispecies for ions as reference species
mi = 1;
me = 1/m;
Ti = 1;
Te = T;
k = 0.2;
D3=@(w)1+ 1/(Te*k^2) * (1+w/(k*sqrt(2*Te/me)) * zetaf(w/(k*sqrt(2*Te/me))) )+...
    1/(Ti*k^2)*(delta*(1+(w/k-v0(1))/sqrt(2*Ti/mi).*zetaf((w/k-v0(1))/sqrt(2*Ti/mi)) ) +...
    (1-delta)*(1+(w/k-v0(2))/sqrt(2*Ti/mi).*zetaf((w/k-v0(2))/sqrt(2*Ti/mi)) ) );
w3=fsolve(D3,ztry,options)

%single species ions with adiabatic electrons
D4=@(w) 1 + Te/Ti*(delta*(1+(w/k-v0(1))/sqrt(2*Ti/mi).*zetaf((w/k-v0(1))/sqrt(2*Ti/mi)) ) +...
    (1-delta)*(1+(w/k-v0(2))/sqrt(2*Ti/mi).*zetaf((w/k-v0(2))/sqrt(2*Ti/mi)) ) );
w4=fsolve(D4,ztry,options)
