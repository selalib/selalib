% Load the particles for all processes (here only two)
%load  particles_start_1_0000.dat
%a1=particles_start_1_0000;
%load  particles_start_1_0001.dat
%a2=particles_start_1_0001;

%% Concatenate the particles from all processes
%a = [a1;a2];
% Load particles form all processes np
np = 2;
a = [];
for j=0:np-1
	name = ['particles_start_1_',num2str(j,'%04d'),'.dat'];
	b=load(name);
	a =[a;b];
end

n = length(a);

nbins = 10; % resolution of histogram 
ind = 2; % index for variable: 1 for x, 2 for v_1, 3 for v_2 

% Create the weighted histogram
[d,c]= histwc(a(:,ind),a(:,4),nbins);
d = d/n;
% Plot the histogram
plot(c,d)
