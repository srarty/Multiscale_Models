% Loads LIF values of balance, synchronization index and coefficient of variation, 
% plots them in colormaps to compare with NMM plotted in Figure_seizure.m
%
% Artemio - April 2023

% close all

% folder = 'C:\Users\artemios\Documents\Multiscale_Models_Data\Spartan\e_vs_ii\';
% folder = 'C:\Users\artemios\Documents\Multiscale_Models_Data\Spartan\e_vs_i\';
% folder = 'C:\Users\artemios\Documents\Multiscale_Models_Data\2023\e_vs_i\';
folder = 'C:\Users\artemios\Documents\Multiscale_Models_Data\2023\e_vs_i_fano\';

d = dir([folder '*_e*.mat']); % Load all files with _e in the name

% range = round(0:0.05:2, 2, 'significant');
range = round(0:0.1:2, 2, 'significant');
% range = round(0:0.5:2, 2, 'significant');

params = set_parameters('gabab');

PLOT_LFP = false;
PLOT_FFT = false;

cmap =[0 0 1; 1 1 0; 1 0 0];
state = nan(length(range), length(range));
balance = nan(length(range), length(range));
CV = nan(length(range), length(range));
SI = nan(length(range), length(range));
FF = nan(length(range), length(range));
for i = 1:length(d)
    % Load file
    lif = load([folder d(i).name]);
    
    % Idx
    e_mult = round( str2num( d(i).name(strfind(d(i).name, '_e')+2 :  min(strfind(d(i).name, '_e')+4 , strfind(d(i).name, '_i')-1)) ) , 3, 'significant');
    i_mult = round( str2num( d(i).name(strfind(d(i).name, '_i')+2 :  min(strfind(d(i).name, '_i')+4 , strfind(d(i).name, '.mat')-1)) ) , 3, 'significant');
    
%     if e_mult == 0.3 && i_mult == 0.8
%         disp('dodgy debugging')
%     end
    
    idx_e = find(range == e_mult);
    idx_i = find(range == i_mult);
    
    lif.i_pe = lif.i_pe(5001:10000);
    lif.i_pi = lif.i_pi(5001:10000);
    lif.i_ie = lif.i_ie(5001:10000);
    lif.i_ii = lif.i_ii(5001:10000);
    
    y = -(lif.i_pe - lif.i_pi)/params.g_m_P; % LFP
    t = 1*lif.lfp_dt : lif.lfp_dt : numel(y) * lif.lfp_dt;
    
    % Calculate fft to estimate oscillatory activity
    [~, X_] = fft_plot( y-mean(y), t, [], PLOT_FFT);
    if max(X_) > 0.05
        state(idx_i, idx_e) = 1; % Oscillation
    elseif (mean(lif.R_py(round(end/3):round(2*end/3))) > params.naka.M) || (mean(lif.R_in(round(end/3):round(2*end/3))) > params.nakai.M)
        % If either population saturates:
        state(idx_i, idx_e) = 2; %Saturation
    else
        state(idx_i, idx_e) = 0; % Normal
    end

    if PLOT_LFP
        figure(100);
        cla;
        plot(t, y, 'Color', cmap(state(idx_i, idx_e)+1 ,:));
        ylabel('LFP (nA)');
    end

    %Balance
    nanoamps_scale = 1e9;    
    balance(idx_i, idx_e) = (mean(lif.i_pi) + mean(lif.i_pe)) * nanoamps_scale;
    
    %
    CV(idx_i, idx_e) = lif.cv_py;
    SI(idx_i, idx_e) = lif.si_py;
    FF(idx_i, idx_e) = lif.fano;
end

%%
figure('Position', [300 400 400 300]);
colormap(cmap);    
imagesc(range, range, state);
caxis([0 2]);
xlabel('Excitatory gain');
ylabel('Inhibitory gain');
title('LIF');
ax = gca;
ax.View = ([0 -90]);
ax.FontSize = 12;

figure('Position', [730 400 490 300]);
colormap(flipud(jet));    
imagesc(range, range, balance);
caxis([-0.2 0.2]);
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
%%
figure('Position',[730 400 490 300])
% colormap(flipud(jet));    
imagesc(range, range, CV);
% caxis([0 0.5]);
ax = gca;
ax.View = ([0 -90]);
ax.FontSize = 12;
xlabel('Excitatory gain');
ylabel('Inhibitory gain');
c = colorbar;
c.Label.String = 'Coefficient of variation';
c.Label.FontSize = 12;

figure('Position',[730 400 490 300])
% colormap(flipud(jet));    
% normalize fanofactor
FF_ = FF./max(FF,[],'all');
imagesc(range, range, FF_);
% caxis([0 0.5]);
ax = gca;
ax.View = ([0 -90]);
ax.FontSize = 12;
xlabel('Excitatory gain');
ylabel('Inhibitory gain');
c = colorbar;
c.Label.String = 'Fanofactor';
c.Label.FontSize = 12;


figure('Position',[730 400 490 300])
% colormap(flipud(jet));    
imagesc(range, range, SI);
% caxis([0 0.5]);
ax = gca;
ax.View = ([0 -90]);
ax.FontSize = 12;
xlabel('Excitatory gain');
ylabel('Inhibitory gain');
c = colorbar;
c.Label.String = 'Syncronization index';
c.Label.FontSize = 12;

drawnow

return
%% Code to plot any colormap of the dummy variable z_
% z_ = []; % Dummy
% z_ = nmm.state.no_drug - state; %[]; % Dummy
% z_ = nmm.balance.no_drug - balance; %[]; % Dummy
z_ = CV;


% normalize = @(x, minx, maxx) (x - minx)/(maxx - minx);
% for i = 1:numel(z_)
%     z_(i) = normalize(z_(i), min(z_,[],'all'), max(z_,[],'all'));
% end

figure('Position',[730 400 490 300])
% colormap(flipud(jet));    
imagesc(range, range, z_);
% caxis([0 0.5]);
ax = gca;
ax.View = ([0 -90]);
ax.FontSize = 12;
xlabel('Excitatory gain');
ylabel('Inhibitory gain');
c = colorbar;
c.Label.String = 'Syncronization index';
c.Label.FontSize = 12;

hold
% colormap(cmap);    
contour(range, range, state, [1,1], 'Color', 'k', 'LineWidth', 3);
contour(range, range, nmm.state.no_drug, [1,1], 'Color', 'white', 'LineWidth', 3);

