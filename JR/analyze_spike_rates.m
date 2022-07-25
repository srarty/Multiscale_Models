%% Comparison of the NMM and LIF based on spike rates.
%
% Needs to be run after NMM_diff_equations_DblExp_recursive.m (or any
% other model).
%
% Artemio - July 2022

close all;

data_file = 'C:/Users/artemios/Documents/Multiscale_Models_Data/spartan/lfp_279.mat'; % u=[0, 0.2, 0, 0.4]
% data_file = 'C:/Users/artemios/Documents/Multiscale_Models_Data/spartan/lfp_280.mat'; % u=0.1
% data_file = 'C:/Users/artemios/Documents/Multiscale_Models_Data/spartan/lfp_281.mat'; % u=0.1 , non-smooth firing rate
% data_file = 'C:/Users/artemios/Documents/Multiscale_Models_Data/spartan/lfp_282.mat'; % u=[0, 1, 0, 10]
% data_file = 'C:/Users/artemios/Documents/Multiscale_Models_Data/spartan/lfp_283.mat'; % u=0
% data_file = 'C:/Users/artemios/Documents/Multiscale_Models_Data/spartan/lfp_285.mat'; % u=[0, 1, 0, 10], non-smooth
% data_file = 'C:/Users/artemios/Documents/Multiscale_Models_Data/spartan/lfp_286.mat'; % u=10, non-smooth
% data_file = 'C:/Users/artemios/Documents/Multiscale_Models_Data/spartan/lfp_287.mat'; % u=15, non-smooth
% data_file = 'C:/Users/artemios/Documents/Multiscale_Models_Data/spartan/lfp_288.mat'; % u=[0, 1, 0, 15], non-smooth

load(data_file);
tt = lfp_dt:lfp_dt:lfp_dt*length(v_ip);

figure
plot(t,x(:,[1 3])); hold on;
plot(tt,v_pi*0.5e3,'b--');
plot(tt,v_ip*0.5e3,'r--');
legend({'x1' 'x3'});
ylabel('Membrane Potential (mV)');
xlabel('Time (s)');
legend({'NMM_{x1}' 'NMM_{x3}' 'LIF_{x1}' 'LIF_{x3}'});

figure
plot(t, f_e); hold on;
plot(t, f_i);
plot(tt,R_py,'--b');
plot(tt,R_in,'--r');
ylabel('Spike rate (Hz)');
xlabel('Time (s)');
legend({'NMM_{Py}' 'NMM_{In}' 'LIF_{Py}' 'LIF_{In}'});