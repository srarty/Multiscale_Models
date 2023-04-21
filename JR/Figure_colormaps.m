%% Goes through results from LIF model and plots colormaps of firing rates,
% run after colormaps_nmm to also plot colormaps from NMM, comment 'return'
% line at the end of the LIF for NMM colormaps.
%
% Artemio updated 05/04/2023

% clear
close all

% function to calculate maximum value of the color code
% max_fun = @(a) max(a(a < max(a)));
max_fun = @(a) max(median(a));
% max_fun = @(a) max(max(a));

%% LIF
% folder = 'C:\Users\artemios\Documents\Multiscale_Models_Data\2023\excitability_spartan\no_current_pulse\'; 
% folder = 'C:\Users\artemios\Documents\Multiscale_Models_Data\2023\excitability_spartan\'; 
% folder = 'C:\Users\artemios\Documents\Multiscale_Models_Data\spartan\firing_rates\jpb\';
folder = 'C:\Users\artemios\Documents\Multiscale_Models_Data\spartan\firing_rates\';


var_vec = {'in_j_AMPA' 'py_j_GABA_' 'py_j_AMPA' 'in_j_GABA_' 'py_j_GABAb_'};
label_vec = {'j^{AMPA}_i' 'j^{GABAA}_p' 'j^{AMPA}_p' 'j^{GABA}_i' 'j^{GABAB}_p'};
% var_vec = {'py_j_AMPA'};
max_values = zeros(2,5);
for ii = 1:length(var_vec)
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


    %% Plot spike rates mesh
    % if results.range2(end)<0, angle=[0 -90]; else, angle=[0 90]; end
    angle = [0 -90];
    
    f = figure('Position', [325 404 1112 364]);
    f.Position = [589 411 758 228];
    
    ax = subplot(1,2,1);
    imagesc(values_(:,1), u_vector, fr_py');% , 'FaceColor', 'flat', 'EdgeColor', 'none')
    ylabel('u');
    xlabel(label_vec{ii});
    % Define colormap
    load custom_colormap
%     cmap = custom_colormap;
    cmap = colormap(jet(512)); % use flipud(jet) to invert the colors in the colormap, jet can be any colormap
    cmap(1:50,:) = [];
%     cmap(end-5:end,:) = [];
%     cmap(end:end+size(cmap,1),:) = ones(1+size(cmap,1),1) * cmap(end,:);
    colormap(cmap); % apply custom colormap
    c = colorbar;
    c.Label.String = 'Mean firing rate (Hz)';
    max_values(1,ii) = 3*max_fun(fr_py); % TODO: 
    caxis([0 1.2*max_values(1,ii)]);
    c.Limits = [0 max_values(1,ii)];
    title('Pyramidal')
    ax.View = (angle);
    ax.FontSize = 12;
%     set(gca, 'ColorScale', 'log')

    ax = subplot(1,2,2);
    imagesc(values_(:,1), u_vector, fr_in');% , 'FaceColor', 'flat', 'EdgeColor', 'none')
    ylabel('u');
    xlabel(label_vec{ii});
    title('Inhibitory');
    c = colorbar;
    c.Label.String = 'Mean firing rate (Hz)';
    max_values(2,ii) = 2*max_fun(fr_in);
    caxis([0 1.2*max_values(2,ii)]);
    c.Limits = [0 max_values(2,ii)];
    ax.View = (angle);
    ax.FontSize = 12;
%     set(gca, 'ColorScale', 'log')

end

% return

%% NMM
% Load values
fr_py_e  = zeros(20,20);
fr_py_i  = zeros(20,20);
fr_py_re = zeros(20,20);
fr_py_ri = zeros(20,20);
fr_py_b  = zeros(20,20);

fr_in_e  = zeros(20,20);
fr_in_i  = zeros(20,20);
fr_in_re = zeros(20,20);
fr_in_ri = zeros(20,20);
fr_in_b  = zeros(20,20);

folder = 'C:\Users\artemios\Documents\Multiscale_Models_Data\Spartan\firing_rates_nmm\';


res{1} = load([folder 'alpha_e vs u']);
res{2} = load([folder 'alpha_i vs u']);
res{3} = load([folder 'alpha_re vs u']);
res{4} = load([folder 'alpha_ri vs u']);
res{5} = load([folder 'alpha_b vs u']);

%% Plot
label_vec = {'\alpha_{ip}' '\alpha_{pi}' '\alpha_{pp}' '\alpha_{ii}' '\alpha_{pb}'};
for j = 1:numel(res)
    r = res{j}.results;
    
    angle = [0 -90];

    f = figure('Position', [325 404 1112 364]);
    f.Position = [589 760 758 228];
    ax = subplot(1,2,1);

    imagesc(r.range, r.range2, r.freqs_py);% , 'FaceColor', 'flat', 'EdgeColor', 'none')
    ylabel('u');
    xlabel(label_vec{j});
    cmap(end,:) = [0 0 0]; % set values greater than the maximum to black
    colormap(cmap);% apply custom colormap
    c = colorbar;
    c.Label.String = 'Mean firing rate (Hz)';
    caxis([0 1.2*max_values(1,j)]);
    c.Limits = [0 max_values(1,j)];
    title('Pyramidal')
    ax.View = (angle);
    ax.FontSize = 12;

    ax = subplot(1,2,2);
    imagesc(r.range, r.range2, r.freqs_in);
    ylabel('u');
    xlabel(label_vec{j});
    title('Inhibitory');
    c = colorbar;
    c.Label.String = 'Mean firing rate (Hz)';
    caxis([0 1.2*max_values(2,j)]);
    c.Limits = [0 max_values(2,j)];
    ax.View = (angle);
    ax.FontSize = 12;
    
end