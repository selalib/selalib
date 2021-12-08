%Jeans instability
sigma_1=1;
k=0.2;
D1=@(x)1-1/k^2*(1+x/(sqrt(2)*sigma_1*k).*zetaf(x/(sqrt(2)*k*sigma_1)));
%D1=@(x)-k^2/k_J^2+(1+x/(sqrt(2)*sigma_1*k).*plasma_dispersion_bot(x/(sqrt(2)*k*sigma_1)));

delta=0.25;
sigma_2=1/sqrt(10);
r=sigma_1/sigma_2;
k_J=1/sigma_1;
k=0.8*k_J;
D2=@(x)-k^2/k_J^2+delta*(1+x/(sqrt(2)*k*sigma_1).*zetaf(x/(sqrt(2)*k*sigma_1)))...
    +r^2*(1-delta)*(1+x/(sqrt(2)*k*sigma_2).*zetaf(x/(sqrt(2)*k*sigma_2)));

ztry=1.0-1i*.000;
%options = optimoptions('fsolve','Display','none','PlotFcn',@optimplotfirstorderopt);
options = optimoptions('fsolve','Display','none');

w1=fsolve(D1,ztry,options)
w2=fsolve(D2,ztry,options)