% Loads LIF values of balance, synchronization index and coefficient of variation,
% plots them in colormaps to compare with NMM plotted in Figure_seizure.m
%
% Artemio - April 2023

% close all

% folder = 'C:\Users\artemios\Documents\Multiscale_Models_Data\Spartan\e_vs_ii\';
% folder = 'C:\Users\artemios\Documents\Multiscale_Models_Data\Spartan\e_vs_i\';
% folder = 'C:\Users\artemios\Documents\Multiscale_Models_Data\2023\e_vs_i\';
% folder = 'C:\Users\artemios\Documents\Multiscale_Models_Data\2023\e_vs_i_fano\'; type_of_LIF = 'CUBN';
folder = 'C:\Users\artemios\Documents\Multiscale_Models_Data\2023\e_vs_i_fano_cobn\'; type_of_LIF = 'COBN';
% folder = 'C:\Users\artemios\Documents\Multiscale_Models_Data\Spartan\e_vs_i_highexc\';
% folder = 'C:\Users\artemios\Documents\Multiscale_Models_Data\2023\e_vs_i_highexc_2\';
% folder = 'C:\Users\artemios\Documents\Multiscale_Models_Data\2023\e_vs_i_highexc_cobn\';

d = dir([folder '*_e*.mat']); % Load all files with _e in the name

% range = round(0:0.05:2, 2, 'significant');
% range = round(0.4:0.1:4, 2, 'significant');
range = round(0.5:0.1:2, 2, 'significant');
% range = round(0:0.5:2, 2, 'significant');

params = set_parameters('gabab');

PLOT_LFP = false;
PLOT_FFT = false;

% cmap =[0 0 1; 1 1 0; 1 0 0];
cmap =[0 0 0; 0 0 1; 1 0 1; 1 0.8 0; 1 0 0];
state = nan(length(range), length(range));
balance = nan(length(range), length(range));
CV = nan(length(range), length(range));
SI = nan(length(range), length(range));
FF = nan(length(range), length(range));
PY_FR = nan(length(range), length(range));
IN_FR = nan(length(range), length(range));
for i = 1:length(d)
% for i = 1629
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
    if isempty(idx_e) || isempty(idx_i)
        continue
    end
    
    lif.i_pe = lif.i_pe(5001:10000);
    lif.i_pi = lif.i_pi(5001:10000);
    lif.i_ie = lif.i_ie(5001:10000);
    lif.i_ii = lif.i_ii(5001:10000);
    
    y = -(lif.i_pe - lif.i_pi)/params.g_m_P; % LFP
    t = 1*lif.lfp_dt : lif.lfp_dt : numel(y) * lif.lfp_dt;
    
    % Calculate fft to estimate oscillatory activity
    if PLOT_FFT, fig_101 = figure(101); cla; else, fig_101 = []; end
    [~, X_, F_] = fft_plot( y-mean(y), t, fig_101, PLOT_FFT);
    if (mean(lif.R_in(500:end)) <= params.nakai.M * 5e-4) || (mean(lif.R_py(500:end)) <= params.naka.M * 5e-4) % 0.05% of the maximum firing rate
        state(idx_i, idx_e)  = -1; % Low state
    elseif max(X_) > 5e-3%0.05
        % state(idx_i, idx_e) = 1; % Oscillation
        [~,indice] = max(X_);
        if F_(indice) < 25
            % Alpha-ish, sleep, low freq
            state(idx_i, idx_e) = 1; % Oscillation
        else
%             PLOT_LFP = true;
            state(idx_i, idx_e) = 2; % Oscillation
            % Gamma-ish, fast oscillations ~60 Hz
        end
    elseif (mean(lif.R_py(round(end/3):round(2*end/3))) > 35) && (mean(lif.R_in(round(end/3):round(2*end/3))) > 60)
        % If either population saturates:
        state(idx_i, idx_e) = 3; %Saturation
    else
        state(idx_i, idx_e) = 0; % Normal
    end

    if PLOT_LFP
        figure(100);
        cla;
        plot(t, y*1e3, 'Color', cmap(state(idx_i, idx_e)+2 ,:));
        ylabel('LFP (mV)');
        F_(indice)
%         fff=figure; plot(F_,X_); xlim([0 100]);
%         PLOT_LFP = false; close(fff);
    end

    %Balance
    nanoamps_scale = 1e9;    
    balance(idx_i, idx_e) = (mean(lif.i_pi) + mean(lif.i_pe)) * nanoamps_scale;
    
    %
    CV(idx_i, idx_e) = lif.cv_py;
    SI(idx_i, idx_e) = lif.si_py;
    IN_FR(idx_i, idx_e) = mean(lif.R_in(500:end));
    PY_FR(idx_i, idx_e) = mean(lif.R_py(500:end));
%     FF(idx_i, idx_e) = lif.fano;
end

%%
figure('Position', [300 400 400 300]);
colormap(cmap);    
imagesc(range, range, state);
caxis([-1 3]);
xlabel('Excitatory gain');
ylabel('Inhibitory gain');
title(type_of_LIF);
ax = gca;
ax.View = ([0 -90]);
ax.FontSize = 12;
% ax.XTick = [0.4:1.2:4];
% ax.YTick = [0.4:1.2:4];

figure('Position', [730 400 490 300]);
load('custom_colormap_parula.mat')
colormap(parula_custom);
imagesc(range, range, balance);
% caxis([-0.25 0.25]);
caxis([-0.5 0.5]);
xlabel('Excitatory gain');
ylabel('Inhibitory gain');
title(type_of_LIF);
ax = gca;
ax.FontSize = 12;
ax.View = ([0 -90]);
% ax.XTick = [0.4:1.2:4];
% ax.YTick = [0.4:1.2:4];
c = colorbar;
c.Label.String = 'Input current balance (nA)';
c.Label.FontSize = 12;

drawnow
%%
figure('Position',[730 400 490 300])
% colormap(flipud(jet));    
imagesc(range, range, CV);
caxis([0 1]);
ax = gca;
ax.View = ([0 -90]);
ax.FontSize = 12;
% ax.XTick = [0.4:1.2:4];
% ax.YTick = [0.4:1.2:4];
xlabel('Excitatory gain');
ylabel('Inhibitory gain');
c = colorbar;
c.Label.String = 'Coefficient of variation';
c.Label.FontSize = 12;

% figure('Position',[730 400 490 300])
% % colormap(flipud(jet));    
% % normalize fanofactor
% FF_ = FF./max(FF,[],'all');
% imagesc(range, range, FF_);
% % caxis([0 0.5]);
% ax = gca;
% ax.View = ([0 -90]);
% ax.FontSize = 12;
% xlabel('Excitatory gain');
% ylabel('Inhibitory gain');
% c = colorbar;
% c.Label.String = 'Fanofactor';
% c.Label.FontSize = 12;


figure('Position',[730 400 490 300])
% colormap(flipud(jet));    
imagesc(range, range, SI);
caxis([0 1]);
ax = gca;
ax.View = ([0 -90]);
ax.FontSize = 12;
% ax.XTick = [0.4:1.2:4];
% ax.YTick = [0.4:1.2:4];
xlabel('Excitatory gain');
ylabel('Inhibitory gain');
c = colorbar;
c.Label.String = 'Syncronization index';
c.Label.FontSize = 12;

drawnow


figure('Position', [175 324 490 612])
% Inhibitory firing rate
ax = subplot(211);   
colormap(jet)
imagesc(range, range, IN_FR);xlabel('Excitatory gain');ylabel('Inhibitory gain');
max_fun = @(a,b) max(2*median(a,'all'), b); % maximum between 3 times the median of the LIF and the first nmm's lower than the maximum firing rate    
maximum_value = max_fun(IN_FR, 0);
caxis([0 4]);
title('Interneurons');
ax.FontSize = 12;
ax.View = ([0 -90]);
% ax.XTick = [0.4:1.2:4];
% ax.YTick = [0.4:1.2:4];
c = colorbar;
c.Label.String = 'Interneurons mean firing rate';
c.Label.FontSize = 12;
drawnow
% Excitatory firing rate
ax = subplot(212);    
colormap(jet)
imagesc(range, range, PY_FR);xlabel('Excitatory gain');ylabel('Inhibitory gain');
maximum_value = max_fun(PY_FR, 0);
caxis([0 0.4]);
title('Pyramidal');
ax.FontSize = 12;
ax.View = ([0 -90]);
% ax.XTick = [0.4:1.2:4];
% ax.YTick = [0.4:1.2:4];
c = colorbar;
c.Label.String = 'Pyramidal mean firing rate';
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

