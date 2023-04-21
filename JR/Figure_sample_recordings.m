% Runs NMM and loads LIF results with Normal, Oscillations and Saturation type of LFP. It plots a histogram of the LFPs
%
% Key values: Saturation (alpha_i = 1, alpha_e = 0.5), Normal_high 
%       excitability (alpha_i = 1, alpha_e = 1), Oscillation (alpha_i = 1, 
%       alpha_e = 1.8), High Oscillation (alpha_i = 1.5, alpha_e = 0.5), 
%       with HIGH_EXC as true. i.e. alpha_e = 0.5, then it is: 
%       1.3 * 1.8 * default_parameter.
%
% LIF folder: 'C://Users/artemios/Documents/Multiscale_Models_Data/'
% LIF Files: Normal (high excitabiltiy): lfp_py__0_u0.mat
%            Oscillation: lfp_py__1_u0.mat
%            Low oscillations: lfp_py__4_u0.mat
%            Saturation: lfp_py__2_u0.mat
%            Default params: lfp_py__5_u0.mat
%
%
% Artemio - April 2023

%% Default
lif = load(['C://Users/artemios/Documents/Multiscale_Models_Data/' 'lfp_py__5_u0.mat']);
[~,~,t,f_e,f_i, params, yy]=NMM_GABA('u',0);
plot_lif_and_nmm(lif, t, f_e, f_i, yy, params, 'Normal');

%% Normal High Excitability
lif = load(['C://Users/artemios/Documents/Multiscale_Models_Data/' 'lfp_py__0_u0.mat']);
[~,~,t,f_e,f_i, params, yy]=NMM_GABA('u',0,'alpha_e',1.3, 'alpha_ri', 1.3);
plot_lif_and_nmm(lif, t, f_e, f_i, yy, params, 'Normal (high exc)');

%% Oscillations
lif = load(['C://Users/artemios/Documents/Multiscale_Models_Data/' 'lfp_py__1_u0.mat']);
[~,~,t,f_e,f_i, params, yy]=NMM_GABA('u',0,'alpha_e',(1.3 * 1.8), 'alpha_ri', 1.3);
plot_lif_and_nmm(lif, t, f_e, f_i, yy, params, 'Oscillations');

%% Low freq oscillations
lif = load(['C://Users/artemios/Documents/Multiscale_Models_Data/' 'lfp_py__4_u0.mat']);
[~,~,t,f_e,f_i, params, yy]=NMM_GABA('u',0,'alpha_e',(1.3 * 0.5), 'alpha_ri', 1.3, 'alpha_i', 1.5);
plot_lif_and_nmm(lif, t, f_e, f_i, yy, params, 'Low freq oscillations');

%% Saturation
lif = load(['C://Users/artemios/Documents/Multiscale_Models_Data/' 'lfp_py__2_u0.mat']);
[~,~,t,f_e,f_i, params, yy]=NMM_GABA('u',0,'alpha_e',(1.3 * 0.5), 'alpha_ri', 1.3);
plot_lif_and_nmm(lif, t, f_e, f_i, yy, params, 'Saturation');


