 %Octave/Matlab script to analyse results of PIF


clear all;

cd '/tmp/landau_sumB'

prefix='piftest';
prefix='landau_B';

prefix='landau_sum_B10'
prefix='landau_B1'
prefix='landau_B10'

gamma=-0.4257; %B=0;
prefix='landau_sumB-0.1.nml';gamma=-0.6665; %B=0.1
prefix='landau_sumB-1.nml';gamma=-0.6740; %B=1
prefix='landau_sumB-10.nml';gamma=-0.4277; %B=10
prefix='landau_sumB-100.nml';gamma=-0.4257; %B=100


prefix='landau_prodB-5';

cd '/home/mrbobble/DATA/Studium/Promotion/helios/B1/results'
prefix='landau_diagB-1';modeidx=4;


cd '/home/mrbobble/DATA/Studium/Promotion/helios/B1Analysis'
prefix='B100Analysis';


DATA = csvread([prefix,'_result.csv'],2);

% dlmread([prefix,'_result.csv'],',');


time=DATA(:,1);
l2potential=DATA(:,6);
energy_error=DATA(:,4);

%modeidx=13;



rho=[];phi=[];
for tstep=1:length(time)
h5filename=sprintf('%s_data_%04d.h5',prefix,tstep);
rho(:,tstep)=h5read(h5filename,'/rho');
phi(:,tstep)=h5read(h5filename,'/phi');
end
unitmode=importdata([prefix,'_phi_unitmodes.dat']);

close all;
figure;
semilogy(time, l2potential); hold on;
semilogy(time, l2potential(1)*exp(gamma.*time))


semilogy(time,phi(15,:));

close all;
%------------- Fourier transform in time
phim=[fliplr(phi), phi   ];
phim=abs(fft(phim,[],2));
dt=time(2)-time(1);
 phim=phim(:,2:end/2);
omega=[0, 2*pi./( (size(phim,2):-1:2)*dt)];

T= (size(phim,2):-1:1)*dt

plot(T)

modeidx=15
plot(omega)
for idx=1:size(phim,1)
semilogy(time,phi(idx,:));
%    semilogy(phim(idx,:).')
legend(sprintf('mode= %d %d %d', unitmode(idx,:)))
pause(2);
end    
    phim=phim(modeidx,:)

semilogy(phim(15,:))

semilogy(omega,phim.')
xlabel('omega');


semilogy(phim.')

semilogy(phi(15,:).')


% 
% fphim=fft(log(phim));
% fphim(40:end)=0;
% phim=ifft(fphim,'symmetric');
% 
% semilogy(exp(phim))
% 
% semilogy(fphim(1:100));hold on;
% 
% 
% fun=@(x)exp(-1*x).*cos(1*x).^2;
% 
% 
% fphim=abs(fft((fun(time))));
% semilogy(fphim(2:40));hold on;
%---------------------------------------

semilogy(time,fun(time));

plot(abs(fft(fun(time))))



figure;
semilogy(time, energy_error)

close all;
figure;
phim=phi(modeidx,:);
semilogy(time,phim); hold on;

semilogy(time,phi);
%Suppose dt=const. we can do an fft over time
% phim_f=abs(fft(log(phim)));
% plot(phim_f(2:200));


%Determine linear growth or decay
[pks,locs]=findpeaks(log(phim),...
         'MINPEAKDISTANCE', ...
         floor(length(time)/15));
%time window
locs=locs(time(locs)>=fit_range.tmin & time(locs)<=fit_range.tmax );
pks=pks(time(locs)>=fit_range.tmin & time(locs)<=fit_range.tmax );
%plot(time(locs),phim(locs),'k^','markerfacecolor',[1 0 0])

% %Growth rate
%mean(diff(pks).'./diff(time(locs)))
omegafit = fit(time(locs),log(phim(locs)).','A + omega*x',...
    'StartPoint',[log(phim(1)), -0.5]);


omegafit.omega

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
clear all; close all;

cd '/home/mrbobble/DATA/Studium/Promotion/helios/results/landau_sumB'
prefix='landau_sumB';modeidx=4;

clear all; close all;

cd '/home/mrbobble/DATA/Studium/Promotion/helios/results/landau_prodB'
prefix='landau_prodB';modeidx=15;

cd '/home/mrbobble/DATA/Studium/Promotion/helios/results/landau2_prodB'
prefix='landau2_prodB';modeidx=15;

cd '/home/mrbobble/DATA/Studium/Promotion/helios/results/landau_diagB/'
prefix='landau_diagB';modeidx=15;

%Fit
fit_range.tmin=1;
fit_range.tmax=8;

mkdir('./plots')

%----------------------
flist=dir([prefix,'*_result.csv']);

fnum=length(flist);
B = zeros(fnum,1);
omega = zeros(fnum,1);
omega_err = zeros(fnum,1);
omega_freq = zeros(fnum,1);
for idx=1:fnum
prefix_loc=strrep(flist(idx).name,'_result.csv','');
    
B(idx)=str2double(strrep(strrep(prefix_loc,[prefix,'-'],''),'.nml',''));

DATA = csvread([prefix_loc,'_result.csv'],2);
time=DATA(:,1);
l2potential=DATA(:,6);
energy_error=DATA(:,4);
clear DATA;

figure;
semilogy(time, l2potential); hold on;


rho=[];phi=[];
for tstep=1:length(time)
h5filename=sprintf('%s_data_%04d.h5',prefix_loc,tstep);
rho(:,tstep)=h5read(h5filename,'/rho');
phi(:,tstep)=h5read(h5filename,'/phi');
end


fitidx=(time>=fit_range.tmin & time<=fit_range.tmax );
phim=phi(modeidx,:);


%Determine linear growth or decay
[pks,locs]=findpeaks(log(phim),...
         'MINPEAKDISTANCE', ...
         floor(length(time)/15));
%time window
fitidxloc=(time(locs)>=fit_range.tmin & time(locs)<=fit_range.tmax );
locs=locs(fitidxloc);
pks=pks(fitidxloc);
omega_freq(idx)=mean(diff(time(locs)));

% %Growth rate
%mean(diff(pks).'./diff(time(locs)))
% omegafit = fit(time(fitidx),log(phim(fitidx)).','A + omega*x',...
%     'StartPoint',[log(phim(1)), -0.5]);

omegafit = fit(time(locs),log(phim(locs)).','A + omega*x',...
    'StartPoint',[log(phim(1)), -0.5])



semilogy(time,phim); hold on;
plot(time, exp(omegafit.A + omegafit.omega*time) );
plot(time(locs),phim(locs),'k^','markerfacecolor',[1 0 0])
ylabel('mode');
title(sprintf('Testcase %s, |B|=%g',strrep(prefix,'_','-'),B(idx)) )
try
    unitmode=importdata([prefix_loc,'_phi_unitmodes.dat']);
    legend('Potential \Phi', ...
        ['mode ', sprintf('%d ', unitmode(modeidx,:))], ...
        'damping fit')
catch
    legend('Potential \Phi', 'selected mode', 'damping fit')
end


omega(idx)=omegafit.omega;
try
omegaconfident=confint(omegafit,0.99);
omega_err(idx)=abs(diff(omegaconfident(:,2)))/2;
catch
omega_err(idx)=0;
end

xlabel(sprintf('time,  \\omega=%g',omega(idx)));hold off;
print('-dpng', '-r300', ['plots/',prefix_loc,'.png']);

end

%Sort result
[B,I]=sort(B);
omega=omega(I);
omega_err=omega_err(I);
omega_freq=omega_freq(I);

figure;
errorbar(B(B>0),omega(B>0),omega_err(B>0),'x-')
xlabel('Strength of Magnetic Field B');
title(['Damping Rates for ',strrep(prefix,'_','-')] )
set(gca,'xscale','log');
print('-dpng', '-r300', ['plots/',prefix,'dispersion.png']);

close all;
figure;
semilogx(B(2:end), -omega(2:end),'x-')


close all;
figure;
semilogx(B(2:end), 2*pi./(omega_freq(2:end)*2),'x-')

2*pi./(omega_freq*2)
B(1)
omega(1)





%       log(phim(fitidx)).', ...
%         'A + log(cos(x*omega0  + omega1)^2) -omega3*x',F)



% F = fitoptions('method', 'NonlinearLeastSquares', ...
%       'Algorithm','Trust-Region','StartPoint',[phim(1), 1, 0, 0.5],...
%       'Robust', 'LAR','Normalize','on')
%    
%  
%  lin_landau_fit = fit( time(fitidx),...
%       phim(fitidx).', ...
%          'A*cos(x*omega0  + omega1)*exp(-omega3*x)',F)
% 
%  plot(lin_landau_fit);hold on;
%      plot(time,phim)

% fitidx=time>4 & time<15;
% 
% 
%  F = fitoptions('method', 'NonlinearLeastSquares', ...
%       'Algorithm','Trust-Region','StartPoint',[log(phim(1)), 1.2, 0.1, 0.5],...
%       'Normalize','on','MaxFunEvals',1e4)
%      
%  
%  lin_landau_fit = fit( time(fitidx),...
%       log(phim(fitidx)).', ...
%          'A + log(cos(x*omega0  + omega1)^2) -omega3*x',F)
% 
%  plot(lin_landau_fit);hold on;
%      plot(time,log(phim))

%Perform analysis for Landau Damping
%Determine the peaks
% %Delete the first peaka
% pks=pks(1:end);
% locs=locs(1:end);
% plot(result.time(locs),result.electrostaticenergy(locs),'k^','markerfacecolor',[1 0 0])
% 
% 
%landau_envelope=@(t)exp(-0.1533*2*t)*(4*0.01*0.3677*sqrt(2*pi))^2 
% 
% %plot(result.time,landau_envelope(result.time))
% 
% %exp(-0.1533*2*time)*(4*0.01*0.3677*sqrt(2*pi))^2 
% 
% lin_landau_fit = fit( result.time(locs),...
%     log(result.electrostaticenergy(locs)), ...
%     'gamma*2*x +offset',...
%     'StartPoint',[-1 0])


%unitmode(modeidx,:)
% 
% 
 for idx=1:size(phi,1)
 semilogy(time,phi(idx,:));
 title(unitmode(idx,:));
 idx
 
 pause(2);
 end
% 
% phi(:,1)
% 
% size(time)
% size(MODES)
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
% % 
% 
% digits(32)
% vpa(sym('1/(2-2^(1/3))'))
% vpa(sym('1-2/(2-2^(1/3))'))
% vpa(sym('(2^(1.0/3.0) +2^(-1/3)-1)/6'))
