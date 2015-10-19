 %Octave/Matlab script to analyse results of PIF
 
cd '/tmp/landau_prod'
cd '/tmp/landau_sum'

clear all;

cd '~/DATA/Studium/Promotion/selalib_build/'
prefix='piftest';
prefix='landau_B';

prefix='landau_sum_B10'
prefix='landau_B1'
prefix='landau_B10'

DATA = csvread([prefix,'_result.csv'],2);

% dlmread([prefix,'_result.csv'],',');


time=DATA(:,1);
l2potential=DATA(:,6);
energy_error=DATA(:,4);

gamma=-0.4531; %B=1
gamma = -0.3972; %B=10


close all;
figure;
semilogy(time, l2potential); hold on;
semilogy(time, l2potential(1)*exp(gamma.*time))

figure;
semilogy(time, energy_error)

rho=[];phi=[];
for tstep=1:length(time)
h5filename=sprintf('%s_data_%04d.h5',prefix,tstep);
rho(:,tstep)=h5read(h5filename,'/rho');
phi(:,tstep)=h5read(h5filename,'/phi');
end


for idx=1:size(phi,1)
semilogy(time,phi(idx,:));
idx
pause(2);
end

phi(:,1)

size(time)
size(MODES)
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

digits(32)
vpa(sym('1/(2-2^(1/3))'))
vpa(sym('1-2/(2-2^(1/3))'))
vpa(sym('(2^(1.0/3.0) +2^(-1/3)-1)/6'))
