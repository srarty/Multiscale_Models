%% CALCULATE_PSP_AGONISTS does the same as fit_alpha, but altering the IPSC
% and time constants to account for GABAergic agonists Diazepam and
% Muscimol according to [1].
% 
% [1] L.Wang, M.Kloc, E.Maher, A.Erisir and A.Maffei (2019) Presynaptic
%     GABAA receptors modulate thalamocortical inputs in layer 4 of rat V1.
%     Cerebral Cortex.
%
% Artemio - April 2022


% OPTIONS
synapse = 'iu'; % synapses 'ab', where a=post-synaptic and b=pre-synaptic
agonist = 'muscimol_external'; % 'diazepam' or '', 'muscimol_external' or 'baclofen'


params = set_parameters('gabab');

% Loading the appropriate file and time constants
switch synapse
    case 'pi'
        load('C:\Users\artemios\Documents\Multiscale_Models_Data\pyramidal_ipsp.mat');
        if strcmp('muscimol',agonist)
            psp = 0.19 * ipsp; % Inhibitory (GABA) on pyramidal, the factor was found form the LIF model by multiplying the IPSC by 0.21 and modifying the time constant
            tau_s = 0.87 * params.tau_sp;
        else
            psp = 1.4 * ipsp; % Inhibitory (GABA) on pyramidal, the factor was found form the LIF model by multiplying the IPSC by 1.2 and modifying the time constant
            tau_s = 1.2 * params.tau_sp;
        end
        tau_m = params.tau_mp;
    case 'ip'
        load('C:\Users\artemios\Documents\Multiscale_Models_Data\inhibitory_epsp.mat'); 
        psp = epsp; % Excitatory (AMPA) on inhibitory interneurons
        tau_s = params.tau_si;
        tau_m = params.tau_mi;
    case 'pp'
        load('C:\Users\artemios\Documents\Multiscale_Models_Data\pyramidal_epsp.mat'); 
        psp = epsp; % Excitatory (AMPA) on pyramidal (recursive)
        tau_s = params.tau_srp;
        tau_m = params.tau_mrp;
    case 'ii'
        load('C:\Users\artemios\Documents\Multiscale_Models_Data\inhibitory_ipsp.mat'); 
        if strcmp('muscimol',agonist)
            psp = 0.18 * ipsp; % Inhibitory (GABA) on inhibitory interneurons (recursive)
            tau_s = 0.87 * params.tau_sri;
        else
            psp = 1.34 * ipsp; % Inhibitory (GABA) on inhibitory interneurons (recursive)
            tau_s = 1.2 * params.tau_sri;
        end
        tau_m = params.tau_mri;
    case 'pb'
        load('C:\Users\artemios\Documents\Multiscale_Models_Data\pyramidalGabab_ipsp.mat'); 
        if strcmp('baclofen', agonist)
            psp = 0.34 * ipsp; % GABAb on pyramidal interneurons
        else
            psp = ipsp; % GABAb on pyramidal interneurons
        end
        tau_s = params.tau_sb;
        tau_m = params.tau_mp;
    case 'iu'
        load('C:\Users\artemios\Documents\Multiscale_Models_Data\inhibitory_externalEPSP.mat'); 
        if strcmp('muscimol_external',agonist)
            psp = 0.50 * epsp; % Excitatory (AMPA_tha) on inhibitory (external input)
        else
            psp = epsp; % Excitatory (AMPA_tha) on inhibitory (external input)
        end
        tau_s = params.tau_si;
        tau_m = params.tau_mi;
    case 'pu'
        load('C:\Users\artemios\Documents\Multiscale_Models_Data\pyramidal_externalEPSP.mat'); 
        if strcmp('muscimol_external',agonist)
            psp = 0.45 * epsp; % Excitatory (AMPA_ext) on pyramidal (external input)
        else
            psp = epsp; % Excitatory (AMPA_ext) on pyramidal (external input)
        end
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