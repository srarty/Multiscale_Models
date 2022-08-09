%% FIT_ALPHA fits an alpha function to the EPSP and IPSP of the Brunel 
% populations
%
% Artemio - February 2022

% load('C:\Users\artemios\Documents\Multiscale_Models_Data\pyramidal_ipsp.mat'); psp = -ipsp; % Inhibitory (GABA) on pyramidal
% load('C:\Users\artemios\Documents\Multiscale_Models_Data\inhibitory_epsp.mat'); psp = epsp; % Excitatory (AMPA) on inhibitory interneurons
% load('C:\Users\artemios\Documents\Multiscale_Models_Data\pyramidal_epsp.mat'); psp = epsp; % Excitatory (AMPA) on pyramidal (recursive)
% load('C:\Users\artemios\Documents\Multiscale_Models_Data\inhibitory_ipsp.mat'); psp = -ipsp; % Inhibitory (GABA) on inhibitory interneurons (recursive)

load('C:\Users\artemios\Documents\Multiscale_Models_Data\taus_fit\pyramidal_ipsp_-140.0.mat'); psp = ipsp;

% Adjuisting psp to remove the first sample, which is zero in Brian2
psp(1) = [];

T = linspace(0,0.3,length(psp));


% ft = fittype( 'b*t*exp(-t/a)', 'independent', 't', 'dependent', 'y' ); % Alpha function | a = tau_mn
ft = fittype( 'c*(exp(-t/a)-exp(-t/b))', 'independent', 't', 'dependent', 'y');

opts = fitoptions(ft);
opts.StartPoint = [0.008339 0.01001 1];
opts.Lower = [0.0001 0.0001 0];     % opts.Lower = [0.0001 0.006031 0];
opts.Upper = [0.5 1 10];    % opts.Upper = [0.1 0.006031 10];    
opts.Robust = 'Off';
fitresult = fit(T', psp', ft, opts) % With options
% fitresult = fit(T', psp', ft) % No options

figure
plot(fitresult)
hold
plot(T,psp,'b--')
legend({'fit' 'brunel'})

% %%
% alpha = fitresult.b; % 277; %
% tau = fitresult.a; % 0.01452; %
% 
% fun = @(t) alpha * t * exp(-t/tau);
% output = zeros(size(T));
% for i = 1:length(T)
%     output(i) = fun(T(i));
% end
% 
% figure
% plot(T, psp);
% hold
% plot(T, output);
% legend({'brunel', 'nmm'});
