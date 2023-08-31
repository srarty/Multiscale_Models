%% Goes through results from LIF model and plots colormaps of firing rates,
% run after colormaps_nmm to also plot colormaps from NMM,
% Set RUN_NMM = false; to skip NMM colormaps.
%
% Artemio updated 05/04/2023

% clear
close all

RUN_NMM = true;
SUBTRACTION = true;
RUN_SCATTER = true;

% function to calculate maximum value of the color code
% max_fun = @(a) max(a(a < max(a)));
% max_fun = @(a) max(median(a));
max_fun = @(a,b) max(2*median(a,'all'), b); % maximum between 3 times the median of the LIF and the first nmm's lower than the maximum firing rate
% max_fun = @(a) max(max(a));

% Synapses:
% var_vec = {'py_j_AMPA'};label_vec_lif={'j^{AMPA}_p'}; label_vec_nmm = {'\alpha_{ip}');
% var_vec_old = {'in_j_AMPA' 'py_j_GABA_' 'py_j_AMPA' 'in_j_GABA_' 'py_j_GABAb_'}; % Old saving style
var_vec = {'_e' '_i' '_re' '_ri' '_b'}; % New saving style (after 15/08/2023)
% var_vec = {'_ri'};
label_vec_lif = {'j^{AMPA}_i' 'j^{GABAA}_p' 'j^{AMPA}_p' 'j^{GABA}_i' 'j^{GABAB}_p'};
label_vec_lif_co = {'g^{AMPA}_i' 'g^{GABAA}_p' 'g^{AMPA}_p' 'g^{GABA}_i' 'g^{GABAB}_p'};
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
% folder = 'C:\Users\artemios\Documents\Multiscale_Models_Data\spartan\firing_rates\'; % Old notation
folder = 'C:\Users\artemios\Documents\Multiscale_Models_Data\spartan\firing_rates_2\'; % New notation (CUBN)
folder_co = 'C:\Users\artemios\Documents\Multiscale_Models_Data\spartan\firing_rates_cobn\'; % New notation (COBN)

% NMM
% Load values
% folder = 'C:\Users\artemios\Documents\Multiscale_Models_Data\Spartan\firing_rates_nmm\';
% folder_nmm = 'C:\Users\artemios\Documents\Multiscale_Models_Data\Spartan\firing_rates_nmm2\';
folder_nmm = 'C:\Users\artemios\Documents\Multiscale_Models_Data\Spartan\firing_rates_nmm3\';

res{1} = load([folder_nmm 'alpha_e vs u']);
res{2} = load([folder_nmm 'alpha_i vs u']);
res{3} = load([folder_nmm 'alpha_re vs u']);
res{4} = load([folder_nmm 'alpha_ri vs u']);
res{5} = load([folder_nmm 'alpha_b vs u']);


%% Run
max_values = zeros(2,5);
for ii = 1:numel(var_vec)
    %% LIF CUBN
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
%             values_(i,jj) = a.value_; % Old notation
            values_(i,jj) = str2double(d(i).name(strfind(d(i).name,var_)+length(var_) : strfind(d(i).name,'_u')-1 )); % New notation (from 08/August/2023)
        end
        
        [~, idx] = sort(values_(:,jj));
%         tr(:,jj) = tr(idx,jj);
        fr_in(:,jj) = fr_in(idx,jj);
        fr_py(:,jj) = fr_py(idx,jj);
        values_(:,jj) = values_(idx,jj); 
        
    end
    
    %% LIF COBN
    for jj = 1:numel(u_vector)
        j = u_vector(jj);
        d = dir([folder_co '*' var_ '*_u' sprintf('%0.1f',j) '*']);

        if isempty(d),continue; end

        % Sort files
        u_value = [];
        for i = 1:length(d), u_value(i) = str2double(d(i).name( strfind(d(i).name,'_u')+2 : strfind(d(i).name,'_u') + 4) ); end
        [~, idx] = sort(u_value);
        d = d(idx);

        for i = 1:length(d)
            % Load file
            a = load([folder_co filesep d(i).name]);
            y_lif_co = a.LFP_V;
            t_lif_co = a.lfp_dt:a.lfp_dt:length(a.LFP_V)*a.lfp_dt;
            % Analyze excitabilty
%             tr(i,jj) = analyze_excitability(y_lif', t_lif', 4899, -1.1, 2500, false);
            L = numel(a.R_in);
            fr_in_co(i,jj) = mean(a.R_in(1500:L-500));
            fr_py_co(i,jj) = mean(a.R_py(1500:L-500));
%             values_(i,jj) = a.value_; % Old notation
            values_co_(i,jj) = str2double(d(i).name(strfind(d(i).name,var_)+length(var_) : strfind(d(i).name,'_u')-1 )); % New notation (from 08/August/2023)
        end
        
        [~, idx] = sort(values_co_(:,jj));
%         tr(:,jj) = tr(idx,jj);
        fr_in_co(:,jj) = fr_in_co(idx,jj);
        fr_py_co(:,jj) = fr_py_co(idx,jj);
        values_co_(:,jj) = values_co_(idx,jj); 
        
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
    
    ax = subplot(2,1+RUN_NMM+RUN_SCATTER,2);
    imagesc(values_(:,1), u_vector, fr_py');% , 'FaceColor', 'flat', 'EdgeColor', 'none')
    ylabel('u');
    xlabel(label_vec_lif{ii});
    ax.Colormap = cmap_lif;
    c = colorbar;
    c.Label.String = 'Mean firing rate (Hz)';
    caxis([0 1.2*max_values(1,ii)]);
    c.Limits = [0 max_values(1,ii)];
    title('CUBN')
    ax.View = (angle);
    ax.FontSize = 12;
    c.FontSize = 12;
%     set(gca, 'ColorScale', 'log')

    ax = subplot(2,1+RUN_NMM+RUN_SCATTER,3+RUN_NMM+RUN_SCATTER);
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


    % Plotting COBN
    ax = subplot(2,1+RUN_NMM+RUN_SCATTER,1);
    imagesc(values_co_(:,1), u_vector, fr_py_co');% , 'FaceColor', 'flat', 'EdgeColor', 'none')
    ylabel('u');
    xlabel(label_vec_lif_co{ii});
    ax.Colormap = cmap_lif;
    c = colorbar;
    c.Label.String = 'Mean firing rate (Hz)';
    caxis([0 1.2*max_values(1,ii)]);
    c.Limits = [0 max_values(1,ii)];
    title('COBN')
    ax.View = (angle);
    ax.FontSize = 12;
    c.FontSize = 12;
%     set(gca, 'ColorScale', 'log')

    ax = subplot(2,1+RUN_NMM+RUN_SCATTER,2+RUN_NMM+RUN_SCATTER);
    imagesc(values_co_(:,1), u_vector, fr_in_co');% , 'FaceColor', 'flat', 'EdgeColor', 'none')
    ylabel('u');
    xlabel(label_vec_lif_co{ii});
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
    
    ax = subplot(2, 1+RUN_NMM+RUN_SCATTER, 2+RUN_SCATTER);
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

    ax = subplot(2, 1+RUN_NMM+RUN_SCATTER, 4+RUN_NMM+RUN_SCATTER);
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
    %{
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
    colorbar_limits =  [-1.2*max_values(1,ii) 1.2*max_values(1,ii)]; %[-max(abs(caxis)) max(abs(caxis))];
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
    colorbar_limits = [-1.2*max_values(2,ii) 1.2*max_values(2,ii)];% [-max(abs(caxis)) max(abs(caxis))];
    caxis(colorbar_limits);
    c.Limits = colorbar_limits;
    ax.FontSize = 12;
    c.FontSize = 12;
    %}
    
    
    
    if RUN_SCATTER
        %%
        figure
        plot([0 4], [0 0], '--', 'Color', [0.5 0.5 0.5]);
        hold
        l1 = plot(r.range, mean((fr_py' - r.freqs_py)./(r.freqs_py + fr_py')),'xb', 'MarkerSize', 8, 'LineWidth', 2);
        l2 = plot(r.range, mean((fr_in' - r.freqs_in)./(r.freqs_in + fr_in')),'xr', 'MarkerSize', 8, 'LineWidth', 2);
        l3 = plot(r.range, mean((fr_py_co' - r.freqs_py)./(r.freqs_py + fr_py_co')),'+b', 'MarkerSize', 8, 'LineWidth', 2);
        l4 = plot(r.range, mean((fr_in_co' - r.freqs_in)./(r.freqs_in + fr_in_co')),'+r', 'MarkerSize', 8, 'LineWidth', 2);
        l = legend([l1(1), l2(1), l3(1), l4(1)], {'CUBN_{Py} - NMM_{Py}' 'CUBN_{In} - NMM_{In}' 'COBN_{Py} - NMM_{Py}' 'COBN_{In} - NMM_{In}'});
        l.Location = 'best';
        
        xlabel(['X ' label_vec_nmm{ii}]);
        ylabel('Normalized firing rate (a.u.)');
        ax = gca;
        ax.FontSize = 12;
        
        %%
        figure
        plot([0 1], [0 1], '--', 'Color', [0.5 0.5 0.5]);
        hold
        % interneurons
        gain_idx = [7 17 37];
        u_idx = 1:length(u_vector); %[1 6 11 16 21];%
        l1 = plot(fr_in(gain_idx,u_idx)'/max(r.freqs_in(u_idx,gain_idx),[],'all'), r.freqs_in(u_idx,gain_idx)/max(r.freqs_in(u_idx,gain_idx),[],'all'), 'xr');
        l2 = plot(fr_in_co(gain_idx,u_idx)'/max(r.freqs_in(u_idx,gain_idx),[],'all'), r.freqs_in(u_idx,gain_idx)/max(r.freqs_in(u_idx,gain_idx),[],'all'), 'or');
        % pyramidal
        l3 = plot(fr_py(gain_idx,u_idx)'/max(r.freqs_py(u_idx,gain_idx),[],'all'), r.freqs_py(u_idx,gain_idx)/max(r.freqs_py(u_idx,gain_idx),[],'all'), 'xb');
        l4 = plot(fr_py_co(gain_idx,u_idx)'/max(r.freqs_py(u_idx,gain_idx),[],'all'), r.freqs_py(u_idx,gain_idx)/max(r.freqs_py(u_idx,gain_idx),[],'all'), 'ob');
        
        l = legend([l1(1), l2(1), l3(1), l4(1)], {'CUBN_{In}' 'COBN_{In}' 'CUBN_{Py}' 'COBN_{Py}'});
        l.Location = 'best';
        
        xlabel('Normalized LIF''s firing rate (a.u.)');
        ylabel('Normalized NMM''s firing rate (a.u.)');
        ax = gca;
        ax.FontSize = 12;
        
        %%
%         figure
%         errorbar(r.range, mean(fr_py' - r.freqs_py), std(fr_py' - r.freqs_py), 'diamondb', 'MarkerSize', 8, 'CapSize', 0, 'LineWidth', 2);
%         hold
%         errorbar(r.range, mean(fr_py_co' - r.freqs_py), std(fr_py_co' - r.freqs_py), 'squareb', 'MarkerSize', 8, 'CapSize', 0, 'LineWidth', 2);
%         errorbar(r.range, mean(fr_in' - r.freqs_in), std(fr_in' - r.freqs_in), 'diamondr', 'MarkerSize', 8, 'CapSize', 0, 'LineWidth', 2);
%         errorbar(r.range, mean(fr_in_co' - r.freqs_in), std(fr_in_co' - r.freqs_in), 'squarer', 'MarkerSize', 8, 'CapSize', 0, 'LineWidth', 2);
%         ylim([-1 1]);
    end
    
    
end