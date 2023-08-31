%% FIT_ALPHA fits an alpha function to the EPSP and IPSP of the Brunel 
% populations
%
% Artemio - February 2022

function [fitresult, T, psp] = fit_alpha(varargin)

if nargin > 0
    synapse = varargin{1};
    PLOT = false;
else
    synapse = 'pi'; % synapses 'ab', where a=post-synaptic and b=pre-synaptic
    PLOT = true;
end

params = set_parameters('gabab');

% load('C:\Users\artemios\Documents\Multiscale_Models_Data\taus_fit\inhibitory_ipsp_-8.2.mat'); psp = ipsp;

% Loading the appropriate file and time constants
switch synapse
    case 'pi'
        load('C:\Users\artemios\Documents\Multiscale_Models_Data\pyramidal_ipsp.mat'); psp = ipsp; % Inhibitory (GABA) on pyramidal
        tau_s = params.tau_sp;
        tau_m = params.tau_mp;
    case 'ip'
        load('C:\Users\artemios\Documents\Multiscale_Models_Data\inhibitory_epsp.mat'); psp = epsp; % Excitatory (AMPA) on inhibitory interneurons
        tau_s = params.tau_si;
        tau_m = params.tau_mi;
    case 'pp'
        load('C:\Users\artemios\Documents\Multiscale_Models_Data\pyramidal_epsp.mat'); psp = epsp; % Excitatory (AMPA) on pyramidal (recursive)
        tau_s = params.tau_srp;
        tau_m = params.tau_mrp;
    case 'ii'
        load('C:\Users\artemios\Documents\Multiscale_Models_Data\inhibitory_ipsp.mat'); psp = ipsp; % Inhibitory (GABA) on inhibitory interneurons (recursive)
        tau_s = params.tau_sri;
        tau_m = params.tau_mri;
    case 'pb'
        load('C:\Users\artemios\Documents\Multiscale_Models_Data\pyramidalGabab_ipsp.mat'); psp = ipsp; % GABAb on pyramidal interneurons
        tau_s = params.tau_sb;
        tau_m = params.tau_mp;
    case 'iu'
        load('C:\Users\artemios\Documents\Multiscale_Models_Data\inhibitory_externalEPSP.mat'); psp = epsp; % Excitatory (AMPA_tha) on inhibitory (external input)
        tau_s = params.tau_si;
        tau_m = params.tau_mi;
    case 'pu'
        load('C:\Users\artemios\Documents\Multiscale_Models_Data\pyramidal_externalEPSP.mat'); psp = epsp; % Excitatory (AMPA_ext) on pyramidal (external input)
        tau_s = params.tau_srp;
        tau_m = params.tau_mp;
end
        
        

% Adjuisting psp to remove the first sample, which is zero in Brian2
psp(1) = [];
T = linspace(0,0.5,length(psp));

% ft = fittype( 'b*t*exp(-t/a)', 'independent', 't', 'dependent', 'y' ); % Alpha function | a = tau_mn
ft = fittype( 'c*(exp(-t/a)-exp(-t/b))', 'independent', 't', 'dependent', 'y');

opts = fitoptions(ft);
opts.StartPoint = [tau_m tau_s 1];
opts.Lower = [tau_m tau_s -20]; %[0.01 0.00525 0]; %0.0024 0];     % opts.Lower = [0.0001 0.006031 0]; 
opts.Upper = [tau_m tau_s 20]; %0.0024 10];    % opts.Upper = [0.1 0.006031 10];    
opts.Robust = 'Off';
fitresult = fit(T', psp', ft, opts) % With options
% fitresult = fit(T', psp', ft) % No options

if PLOT
    figure
    plot(fitresult)
    hold
    plot(T,psp,'b--')
    legend({'fit' 'brunel'}) 
end

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
