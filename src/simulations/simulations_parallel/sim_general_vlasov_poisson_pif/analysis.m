 %Octave/Matlab script to analyse results of PIF
 
cd '/tmp/landau_sum'
clear all;

DATA = csvread('pif_resultB1sum.csv',2);
DATA = csvread('pif_result.csv',2);
%DATA = csvread('pif_resultB0.csv',2);
time=DATA(:,1);
l2potential=DATA(:,6);


gamma=-0.4531;
gamma=-0.4531;

semilogy(time, l2potential); hold on;
semilogy(time, l2potential(1)*exp(0.5*gamma.*time))

% close all;
% figure;

% xlabel('time')
% ylabel('electrostatic energy')
% hold on;
% %Perform analysis for Landau Damping
% %Determine the peaks
% [pks,locs]=findpeaks(log(result.electrostaticenergy),...
%      'MINPEAKDISTANCE', floor(length(result.electrostaticenergy)/10))%, 'MINPEAKHEIGHT',-0.3e-3)
% %Delete the first peaka
% pks=pks(1:end);
% locs=locs(1:end);
% plot(result.time(locs),result.electrostaticenergy(locs),'k^','markerfacecolor',[1 0 0])
% 
% 
% landau_envelope=@(t)exp(-0.1533*2*t)*(4*0.01*0.3677*sqrt(2*pi))^2 
% 
% %plot(result.time,landau_envelope(result.time))
% 
% %exp(-0.1533*2*time)*(4*0.01*0.3677*sqrt(2*pi))^2 
% 
% lin_landau_fit = fit( result.time(locs),...
%     log(result.electrostaticenergy(locs)), ...
%     'gamma*2*x +offset',...
%     'StartPoint',[-1 0])
% plot(  result.time,   exp(lin_landau_fit(result.time)),'r')
% title(sprintf('%f, %f',lin_landau_fit.gamma, exp(lin_landau_fit.gamma)))
% 
% %semilogy(result.time, result.kineticenergy)
% 
% %semilogy(result.time, result.electrostaticenergy+ result.kineticenergy)
% 
% %exp(-0.1533*2*time)*(4*0.01*0.3677*sqrt(2*pi))^2 
% result.electrostaticenergy(1)/result.kineticenergy(1)
% 
% 
% %-----------------------------------------
% figure;
% plot(result.electrostaticenergy)
% figure;
% plot(result.kineticenergy, 'r')
% close all;
% plot(result.kineticenergy+ result.electrostaticenergy*16)
% 
% 
% figure;
% plot(sqrt(result.electrostaticenergy))
% figure;
% plot(result.kineticenergy)
% close all;
% plot(result.kineticenergy.^2+result.electrostaticenergy)
% total=result.kineticenergy/2+result.electrostaticenergy
% semilogy((total-total(1)))
% 
% L=4*pi
% k=0.5
% alpha=0.1
% x=linspace(0,L,100)
% 
% f=@(x)(1+alpha*cos(k*x))/L
% 
% uprime=@(x)(alpha*sin(k*x)/k)/L
% u=@(x)(-alpha*cos(k*x)/(k^2))/L
% 
% plot(x,uprime(x)*L/2)
% plot(x,u(x))
% 
% 
% (min(uprime(x))+1.59154963650877000E-003)/min(uprime(x))
% max(uprime(x))
% 
% 
%  xmax=1.59154942883948892E-003
% 
% max(0.01*sin(1.0*x))
% 
% 
% 0.02/1.34003521107987434E-003
% 
% -1.31848257489253890E-003
% 
% 0.02
% 
% 
% Ef=
% Ekin=
% 
