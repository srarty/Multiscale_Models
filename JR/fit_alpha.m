%% FIT_ALPHA fits an alpha function to the EPSP and IPSP of the Brunel 
% populations
%
% Artemio - February 2022

% load('C:\Users\artemios\Documents\Multiscale_Models_Data\pyramidal_ipsp.mat'); psp = ipsp; % Onhibitory (GABA) on pyramidal
load('C:\Users\artemios\Documents\Multiscale_Models_Data\inhibitory_epsp.mat'); psp = epsp; % Excitatory (AMPA) on inhibitory interneurons
T = linspace(0,0.25,length(psp));


ft = fittype( 'b*t*exp(-t/a)', 'independent', 't', 'dependent', 'y' ); % Alpha function | a = tau_mn

opts = fitoptions(ft);
opts.StartPoint = [0.02 1];
opts.Lower = [0.001 -1000];     
opts.Upper = [0.1 1000];    
% fitresult = fit(T', psp', ft, opts) % With options
fitresult = fit(T', psp', ft) % No options

figure
plot(fitresult)
hold
plot(T,psp)
legend({'fit' 'brunel'})

%%
alpha = fitresult.b;
tau =  fitresult.a;

fun = @(t) alpha * t * exp(-t/tau);
output = zeros(size(T));
for i = 1:length(T)
    output(i) = fun(T(i));
end

figure
plot(T, psp);
hold
plot(T, output);
legend({'brunel', 'nmm'});
