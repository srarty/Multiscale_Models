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

%% Default (no-input, i.e. Spontaneous Activity)
lif = load(['C://Users/artemios/Documents/Multiscale_Models_Data/' 'lfp_py__5_u0.mat']);
[x,~,t,f_e,f_i, params, yy]=NMM_GABA('u',0);
plot_lif_and_nmm(lif, t, f_e, f_i, yy, x, params, 'Normal');

%% Default (u=1)
lif = load(['C://Users/artemios/Documents/Multiscale_Models_Data/' 'lfp_py__4_u1.mat']);
[x,~,t,f_e,f_i, params, yy]=NMM_GABA('u',1);
plot_lif_and_nmm(lif, t, f_e, f_i, yy, x, params, 'Normal u=1');

%% Normal High Excitability
% lif = load(['C://Users/artemios/Documents/Multiscale_Models_Data/' 'lfp_py__0_u0.mat']);
% [~,~,t,f_e,f_i, params, yy]=NMM_GABA('u',0,'alpha_e',1.3, 'alpha_ri', 1.3);
% plot_lif_and_nmm(lif, t, f_e, f_i, yy, x, params, 'Normal (high exc)');

%% Oscillations
lif = load(['C://Users/artemios/Documents/Multiscale_Models_Data/' 'lfp_py__7_u0.mat']);
[x,~,t,f_e,f_i, params, yy]=NMM_GABA('u',0,'alpha_e',(1.3 * 1.8), 'alpha_ri', 1.3);
plot_lif_and_nmm(lif, t, f_e, f_i, yy, x, params, 'Oscillations');

%% Low freq oscillations
lif = load(['C://Users/artemios/Documents/Multiscale_Models_Data/' 'lfp_py__4_u0.mat']);
[x,~,t,f_e,f_i, params, yy]=NMM_GABA('u',0,'alpha_e',(1.3 * 0.5), 'alpha_ri', 1.3, 'alpha_i', 1.5);
plot_lif_and_nmm(lif, t, f_e, f_i, yy, x, params, 'Low freq oscillations');

%% Saturation
lif = load(['C://Users/artemios/Documents/Multiscale_Models_Data/' 'lfp_py__2_u0.mat']);
[x,~,t,f_e,f_i, params, yy]=NMM_GABA('u',0,'alpha_e',(1.3 * 0.5), 'alpha_ri', 1.3);
plot_lif_and_nmm(lif, t, f_e, f_i, yy, x, params, 'Saturation');

%% Changing u
lif = load(['C://Users/artemios/Documents/Multiscale_Models_Data/2023/' 'lfp_e1.00_i1.00.mat']);
[x,~,t,f_e,f_i, params, yy]=NMM_GABA();
plot_lif_and_nmm(lif, t, f_e, f_i, yy, x, params, 'Changing u');

%% Normal excitability, 0.5 concentration baclofen (high freq oscillations)
lif = load(['C://Users/artemios/Documents/Multiscale_Models_Data/2023/' 'lfp_e1.80_i1.00.mat']);
[x,~,t,f_e,f_i, params, yy]=NMM_GABA('u',0,'alpha_e', 1.8, 'alpha_b', 0.67);
plot_lif_and_nmm(lif, t, f_e, f_i, yy, x, params, 'Baclofen oscillations');

%% High excitability, no-drugs, high freq oscillations
lif = load(['C://Users/artemios/Documents/Multiscale_Models_Data/2023/' 'fast_oscillation_0.mat']);
[x,~,t,f_e,f_i, params, yy]=NMM_GABA('u',0,'alpha_e', 4.55 , 'alpha_i', 1.1, 'alpha_ri', 1.3);
plot_lif_and_nmm(lif, t, f_e, f_i, yy, x, params, 'Baclofen oscillations');



%% Loop through a folder
% folder = 'C:\Users\artemios\Documents\Multiscale_Models_Data\2023\e_vs_i\';
folder = 'C:\Users\artemios\Documents\Multiscale_Models_Data\2023\e_vs_i_fano\';
d = dir([folder '*_e*.mat']); % Load all files with _e in the name
range = round(0:0.1:2, 2, 'significant');


for i = 1:length(d)
    % Idx
    e_mult = round( str2num( d(i).name(strfind(d(i).name, '_e')+2 :  min(strfind(d(i).name, '_e')+4 , strfind(d(i).name, '_i')-1)) ) , 3, 'significant');
    i_mult = round( str2num( d(i).name(strfind(d(i).name, '_i')+2 :  min(strfind(d(i).name, '_i')+4 , strfind(d(i).name, '.mat')-1)) ) , 3, 'significant');
    
    if e_mult < 1e-5, e_mult = 0; end
    if i_mult < 1e-5, i_mult = 0; end
    
    % Load file
    lif = load([folder d(i).name]);
    
    [x,~,t,f_e,f_i, params, yy]=NMM_GABA('u',0,'alpha_e', e_mult, 'alpha_i', i_mult);
    plot_lif_and_nmm(lif, t, f_e, f_i, yy, x, params, ['Excitatory = ' num2str(e_mult) ' | Inhibitory = ' num2str(i_mult) ]);
    
end