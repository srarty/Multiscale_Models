% TO DO
% 
% Loads LIF values of balance, synchronization index, coefficient of variation and excitability (optional), 
% plots them in colormaps to compare with NMM plotted in Figure_seizure.m
%
% Artemio - May 2023

% close all
folder = 'C:\Users\artemios\Documents\Multiscale_Models_Data\Spartan\e_vs_i_excitability\';

d = dir([folder '*_e*.mat']); % Load all files with _e in the name

range = round(0:0.1:2, 2, 'significant');

params = set_parameters('gabab');

PLOT_LFP = false;
PLOT_FFT = false;

cmap =[0 0 1; 1 1 0; 1 0 0];
recovery = nan(length(range), length(range));
state = nan(length(range), length(range));
for i = 1:length(d)
    % Load file
    lif = load([folder d(i).name]);
    
    % Idx
    e_mult = round( str2num( d(i).name(strfind(d(i).name, '_e')+2 :  min(strfind(d(i).name, '_e')+4 , strfind(d(i).name, '_i')-1)) ) , 3, 'significant');
    i_mult = round( str2num( d(i).name(strfind(d(i).name, '_i')+2 :  min(strfind(d(i).name, '_i')+4 , strfind(d(i).name, '.mat')-1)) ) , 3, 'significant');
    idx_e = find(range == e_mult);
    idx_i = find(range == i_mult);
    
    lif.i_pe = lif.i_pe;
    lif.i_pi = lif.i_pi;
    lif.i_ie = lif.i_ie;
    lif.i_ii = lif.i_ii;
    
    y = -(lif.i_pe - lif.i_pi)/params.g_m_P; % LFP
    t = 1*lif.lfp_dt : lif.lfp_dt : numel(y) * lif.lfp_dt;
    
    % Calculate fft to estimate oscillatory activity
    [~, X_] = fft_plot( y(6000:10000)-mean(y(6000:10000)), t(6000:10000), [], PLOT_FFT);
    if max(X_) > 0.025
        state(idx_i, idx_e) = 1; % Oscillation
    else
        state(idx_i, idx_e) = 0; % Normal
    end
    
    % Only measure recovery time if there wasn't saturation nor
    % oscillations:
    if ~state(idx_i, idx_e)
        recovery(idx_i, idx_e) = analyze_excitability(y', t', 4899, -1.1, 2500, false, false, 0.3, 150);
    else
        recovery(idx_i, idx_e) = 0.6;
    end
    
end

%%
figure('Position', [730 400 490 300]);
colormap(jet(512));
imagesc(range, range, recovery);
caxis([0 0.6]);
xlabel('Excitatory gain');
ylabel('Inhibitory gain');
title('LIF');
ax = gca;
ax.FontSize = 12;
ax.View = ([0 -90]);
c = colorbar;
c.Label.String = 'Input current balance (nA)';
c.Label.FontSize = 12;

drawnow