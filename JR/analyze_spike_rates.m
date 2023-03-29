%% Comparison of the NMM and LIF based on spike rates.
%
% Needs to be run after NMM_diff_equations_DblExp_recursive.m (or any
% other model).
%
% Artemio - July 2022

close all;

% data_file = 'C:/Users/artemios/Documents/Multiscale_Models_Data/spartan/lfp_312.mat';
% data_file = 'C:/Users/artemios/Documents/Multiscale_Models_Data/lfp_17.mat'; % u=[0 1 3 5]
% data_file = 'C:/Users/artemios/Documents/Multiscale_Models_Data/lfp_18.mat'; % u=[10]
% data_file = 'C:/Users/artemios/Documents/Multiscale_Models_Data/lfp_20.mat'; % u=[0 0.5 1 1.5]
% data_file = 'C:/Users/artemios/Documents/Multiscale_Models_Data/lfp_39.mat'; % u=[0, 0.25, 0.75, 0.5]
% data_file = 'C:/Users/artemios/Documents/Multiscale_Models_Data/lfp_40.mat';
% data_file = 'C:/Users/artemios/Documents/Multiscale_Models_Data/lfp_41.mat'; % u=[20]
% data_file = 'C:/Users/artemios/Documents/Multiscale_Models_Data/lfp_44.mat'; % injected_py=[500]
% data_file = 'C:/Users/artemios/Documents/Multiscale_Models_Data/lfp_45.mat'; % injected_py=[300]
% data_file = 'C:/Users/artemios/Documents/Multiscale_Models_Data/lfp_46.mat'; % injected_in=[300]
% data_file = 'C:/Users/artemios/Documents/Multiscale_Models_Data/lfp_0.mat'; % u=[0, 0.25, 0.75, 0.5]

% data_file = 'C:/Users/artemios/Documents/Multiscale_Models_Data/lfp_61.mat'; % u=[0, 0.25, 0.5, 1], j_ii = 45.24, alpha_ii = -2.5 (NOTE!! I believe this recording in the LIF was computed with FIXED th and refractory period)
% data_file = 'C:/Users/artemios/Documents/Multiscale_Models_Data/lfp_63.mat'; % impulse response, j_ii = 45.24, alpha_ii = -2.5
% data_file = 'C:/Users/artemios/Documents/Multiscale_Models_Data/lfp_64.mat'; % impulse response, j_pi = 18.2577, alpha_i = -0.26
% data_file = 'C:/Users/artemios/Documents/Multiscale_Models_Data/lfp_75.mat'; % impulse response (50 pA), j_pi = 37, alpha_i = -0.5xxx
% data_file = 'C:/Users/artemios/Documents/Multiscale_Models_Data/lfp_80.mat'; % impulse response (100 pA), j_pi = 37, alpha_i = -0.5xxx
% data_file = 'C:/Users/artemios/Documents/Multiscale_Models_Data/lfp_79.mat'; % impulse response (500 pA), j_pi = 37, alpha_i = -0.5xxx
% data_file = 'C:/Users/artemios/Documents/Multiscale_Models_Data/lfp_65.mat'; % impulse response (500 pA), j_pi = 21.0666, alpha_i = -0.3
% data_file = 'C:/Users/artemios/Documents/Multiscale_Models_Data/lfp_66.mat'; % Seizure: j_pi = 21.0666, alpha_i = -0.3 (random LIF th and t_ref)

% data_file = 'C:/Users/artemios/Documents/Multiscale_Models_Data/2023/lfp_py__7_u0.mat';
data_file = 'C:/Users/artemios/Documents/Multiscale_Models_Data/2023/lfp_last.mat';

load(data_file);
tt = lfp_dt:lfp_dt:lfp_dt*length(v_ip);

figure
plot(t,x(:,[1 3])); hold on;
plot(tt,v_pi*1e3,'b--');
plot(tt,v_ip*1e3,'r--');
legend({'x1' 'x3'});
ylabel('Post-synaptic Potential (mV)');
xlabel('Time (s)');
legend({'NMM_{x1}' 'NMM_{x3}' 'LIF_{x1}' 'LIF_{x3}'});
ylim([-70 70])

figure
subplot(211)
plot(t, f_e,'Color', [0 0.45 0.74]); hold on;
subplot(212)
plot(t, f_i,'Color', [0.85 0.33 0.1]);hold on
try
    subplot(211)
    plot(tt,R_py2,'--b', 'LineWidth',2);
    subplot(212)
    plot(tt,R_in2,'--r', 'LineWidth',2);
catch
    subplot(211)
    plot(tt,R_py,'--b', 'LineWidth',2);
    subplot(212)
    plot(tt,R_in,'--r', 'LineWidth',2);
end
% r = movavg(f_i, 'triangular', floor(length(f_i)/10));
% plot(t,r)
% r = movavg(f_e, 'triangular', floor(length(f_e)/10));
% plot(t,r)
ylabel('Spike rate (Hz)');
xlabel('Time (s)');
% legend({'NMM_{Py}' 'NMM_{In}' 'LIF_{Py}' 'LIF_{In}' 'NMM_{Py}_{M.A.}' 'NNM_{In}_{M.A.}'});
subplot(211)
legend({'NMM_{Py}' 'LIF_{Py}'}); % 'NMM_{Py}_{M.A.}' 'NNM_{In}_{M.A.}'});
subplot(212)
legend({'NMM_{In}' 'LIF_{In}'}); % 'NMM_{Py}_{M.A.}' 'NNM_{In}_{M.A.}'});
