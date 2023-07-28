%% Goes through results from LIF model and plots colormaps of firing rates,
% run after colormaps_nmm to also plot colormaps from NMM,
% Set RUN_NMM = false; to skip NMM colormaps.
%
% Artemio updated 05/04/2023

% clear
close all

RUN_NMM = true;
SUBTRACTION = true;

% function to calculate maximum value of the color code
% max_fun = @(a) max(a(a < max(a)));
% max_fun = @(a) max(median(a));
max_fun = @(a,b) max(2*median(a,'all'), b); % maximum between 3 times the median of the LIF and the first nmm's lower than the maximum firing rate
% max_fun = @(a) max(max(a));

% Synapses:
% var_vec = {'py_j_AMPA'};label_vec_lif={'j^{AMPA}_p'}; label_vec_nmm = {'\alpha_{ip}');
var_vec = {'in_j_AMPA' 'py_j_GABA_' 'py_j_AMPA' 'in_j_GABA_' 'py_j_GABAb_'};
label_vec_lif = {'j^{AMPA}_i' 'j^{GABAA}_p' 'j^{AMPA}_p' 'j^{GABA}_i' 'j^{GABAB}_p'};
label_vec_nmm = {'\alpha_{ip}' '\alpha_{pi}' '\alpha_{pp}' '\alpha_{ii}' '\alpha_{pb}'};


% Define colormap
% load custom_colormap
cmap_lif = colormap(jet(512)); % use flipud(jet) to invert the colors in the colormap, jet can be any colormap
cmap_lif(1:50,:) = [];
cmap_nmm = cmap_lif;
cmap_nmm(end,:) = [0 0 0]; % set values greater than the maximum to black
close % colormap function creates a new figure if none exists. This is to close it.

% LIF
% folder = 'C:\Users\artemios\Documents\Multiscale_Models_Data\2023\excitability_spartan\no_current_pulse\'; 
% folder = 'C:\Users\artemios\Documents\Multiscale_Models_Data\2023\excitability_spartan\'; 
% folder = 'C:\Users\artemios\Documents\Multiscale_Models_Data\spartan\firing_rates\jpb\';
folder = 'C:\Users\artemios\Documents\Multiscale_Models_Data\spartan\firing_rates\';

% NMM
% Load values
% folder = 'C:\Users\artemios\Documents\Multiscale_Models_Data\Spartan\firing_rates_nmm\';
folder_nmm = 'C:\Users\artemios\Documents\Multiscale_Models_Data\Spartan\firing_rates_nmm2\';

res{1} = load([folder_nmm 'alpha_e vs u']);
res{2} = load([folder_nmm 'alpha_i vs u']);
res{3} = load([folder_nmm 'alpha_re vs u']);
res{4} = load([folder_nmm 'alpha_ri vs u']);
res{5} = load([folder_nmm 'alpha_b vs u']);


%% Run
max_values = zeros(2,5);
for ii = 1:numel(var_vec)
    %% LIF
    u_vector = 0:0.1:2.0;
    var_ = var_vec{ii};
    for jj = 1:numel(u_vector)
        j = u_vector(jj);
        d = dir([folder '*' var_ '*_u' sprintf('%0.1f',j) '*']);

        if isempty(d),continue; end

        % Sort files
        u_value = [];
        for i = 1:length(d), u_value(i) = str2double(d(i).name( strfind(d(i).name,'_u')+2 : strfind(d(i).name,'_u') + 4) ); end
        [~, idx] = sort(u_value);
        d = d(idx);

        for i = 1:length(d)
            % Load file
            a = load([folder filesep d(i).name]);
            y_lif = a.LFP_V;
            t_lif = a.lfp_dt:a.lfp_dt:length(a.LFP_V)*a.lfp_dt;
            % Analyze excitabilty
%             tr(i,jj) = analyze_excitability(y_lif', t_lif', 4899, -1.1, 2500, false);
            L = numel(a.R_in);
            fr_in(i,jj) = mean(a.R_in(1500:L-500));
            fr_py(i,jj) = mean(a.R_py(1500:L-500));
            values_(i,jj) = a.value_;
        end
        
        [~, idx] = sort(values_(:,jj));
%         tr(:,jj) = tr(idx,jj);
        fr_in(:,jj) = fr_in(idx,jj);
        fr_py(:,jj) = fr_py(idx,jj);
        values_(:,jj) = values_(idx,jj);
        
    end
    
    %% NMM
    if ~RUN_NMM, continue; end % Skip this section of the code (i.e. do not plot NMM results)
    
    r = res{ii}.results;
    

    %% Plot spike rates mesh (LIF)
    % Get maximum value:
    max_nmm_py = max( r.freqs_py(r.freqs_py < (floor(max(r.freqs_py,[],'all')))/2) ); % Highest NMM's firing rate, lower than the max firing rate of Pyramidals
%     max_nmm_py = max( r.freqs_py(r.freqs_py < 30) ); % Highest NMM's firing rate, lower than the max firing rate of Pyramidals
    max_values(1,ii) = max_fun(fr_py, max_nmm_py);
    max_nmm_in = max( r.freqs_in(r.freqs_in < floor(max(r.freqs_in,[],'all'))) ); % Highest NMM's firing rate, lower than the max firing rate of Interneurons
    max_values(2,ii) = max_fun(fr_in, 0);
    
    % if results.range2(end)<0, angle=[0 -90]; else, angle=[0 90]; end
    angle = [0 -90];
    
    f = figure('Position', [420 180 1380 590]);
    
    ax = subplot(2,1+RUN_NMM+SUBTRACTION,1);
    imagesc(values_(:,1), u_vector, fr_py');% , 'FaceColor', 'flat', 'EdgeColor', 'none')
    ylabel('u');
    xlabel(label_vec_lif{ii});
    ax.Colormap = cmap_lif;
    c = colorbar;
    c.Label.String = 'Mean firing rate (Hz)';
    caxis([0 1.2*max_values(1,ii)]);
    c.Limits = [0 max_values(1,ii)];
    title('LIF')
    ax.View = (angle);
    ax.FontSize = 12;
    c.FontSize = 12;
%     set(gca, 'ColorScale', 'log')

    ax = subplot(2,1+RUN_NMM+SUBTRACTION,2+RUN_NMM+SUBTRACTION);
    imagesc(values_(:,1), u_vector, fr_in');% , 'FaceColor', 'flat', 'EdgeColor', 'none')
    ylabel('u');
    xlabel(label_vec_lif{ii});
%     title('Inhibitory');
    ax.Colormap = cmap_lif;
    c = colorbar;
    c.Label.String = 'Mean firing rate (Hz)';
    caxis([0 1.2*max_values(2,ii)]);
    c.Limits = [0 max_values(2,ii)];
    ax.View = (angle);
    ax.FontSize = 12;
    c.FontSize = 12;
%     set(gca, 'ColorScale', 'log')

    % Plotting NMM
    angle = [0 -90];

%     f = figure();%'Position', [325 404 1112 364]);
%     f.Position = [1223 487 362 512];
    
    ax = subplot(2, 2+SUBTRACTION, 2);
    imagesc(r.range, r.range2, r.freqs_py);% , 'FaceColor', 'flat', 'EdgeColor', 'none')
    ylabel('u');
    xlabel(label_vec_nmm{ii});
    ax.Colormap = cmap_nmm;
    c = colorbar;
    c.Label.String = 'Mean firing rate (Hz)';
    caxis([0 1.2*max_values(1,ii)]);
    c.Limits = [0 max_values(1,ii)];
    title('NMM')
    ax.View = (angle);
    ax.FontSize = 12;
    c.FontSize = 12;

    ax = subplot(2, 2+SUBTRACTION, 4+SUBTRACTION);
    imagesc(r.range, r.range2, r.freqs_in);
    ylabel('u');
    xlabel(label_vec_nmm{ii});
%     title('Inhibitory');
    ax.Colormap = cmap_nmm;
    c = colorbar;
    c.Label.String = 'Mean firing rate (Hz)';
    caxis([0 1.2*max_values(2,ii)]);
    c.Limits = [0 max_values(2,ii)];
    ax.View = (angle);
    ax.FontSize = 12;
    c.FontSize = 12;
        
    %% Subtraction
    if ~SUBTRACTION, continue; end
    
    ax = subplot(2,3,3);
    pyramidal = fr_py' - r.freqs_py;
    imagesc(r.range,r.range2,pyramidal)
    angle = [0 -90];
    ax.View = (angle);
    title('Difference');
    ylabel('u');
    xlabel('Multiplier');
    c = colorbar;
    colorbar_limits = [-max(abs(caxis)) max(abs(caxis))];
    caxis(colorbar_limits);
    c.Limits = colorbar_limits;
    ax.FontSize = 12;
    c.FontSize = 12;
    
    ax = subplot(2,3,6);
    interneurons = fr_in' - r.freqs_in;
    imagesc(r.range,r.range2,interneurons);
    angle = [0 -90];
    ax.View = (angle);
    ax.FontSize = 12;
    ylabel('u');
    xlabel('Multiplier');
    c = colorbar;
    colorbar_limits = [-max(abs(caxis)) max(abs(caxis))];
    caxis(colorbar_limits);
    c.Limits = colorbar_limits;
    ax.FontSize = 12;
    c.FontSize = 12;
    
    
end