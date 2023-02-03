% Loads all the LIF files from a folder and gets the excitability according
% to a varying parameter
% 
% It is only for LIF, to get the NMM figure run Figure_excitabilty.m

folder = 'C:\Users\artemios\Documents\Multiscale_Models_Data\2023\excitability\'; 

%% external input
var_ = 'u';
d = dir([folder '*' var_ '*']);
tr_u_lif = zeros(1,length(d));
fr_in_u_lif = zeros(1,length(d));
fr_py_u_lif = zeros(1,length(d));
value = zeros(1,length(d));

% Sort files
dates = [];
for i = 1:length(d), dates(i) = d(i).datenum; end
[~, idx] = sort(dates);
d = d(idx);

% Calculate excitability
for i = 1:length(d)
    % Load file
    a = load([folder filesep d(i).name]);
    y_lif = a.LFP_V;
    t_lif = a.lfp_dt:a.lfp_dt:length(a.LFP_V)*a.lfp_dt;
    % Analyze excitabilty
    tr_u_lif(i) = analyze_excitability(y_lif', t_lif', 4899, -1.1, 2500, false);
    fr_in_u_lif(i) = mean(a.R_in(1000:4000));    
    fr_py_u_lif(i) = mean(a.R_py(1000:4000));    
    value(i) = a.input_spike_rate;
end

% NMM reults
load parameter_sweeps/u;
nmm_result = tr_u;
nmm_in = fr_in_u;
nmm_py = fr_py_u;
title_ = 'External input';

% Plot excitability
figure; 
plot(value, tr_u_lif, 'k--');
hold
plot(values, nmm_result, 'k');
title('Excitability vs external input');
l = legend({'LIF', 'NMM'});
l.Location = 'best';
xlabel('u');

% Plot firing rates
figure;
l1 = plot(value, fr_py_u_lif, '--'); hold
l2 = plot(value, fr_in_u_lif, '--');
plot(values, nmm_py, 'Color', l1.Color);
plot(values, nmm_in, 'Color', l2.Color);
title(['Firing rates vs ' title_]);
l = legend({'Pyramidal (LIF)', 'Inhibitory (LIF)', 'Pyramidal (NMM)', 'Inhibitory (NMM)'});
l.Location = 'best';
xlabel(title_);
ylabel('Firing rate (Hz)');

%% Gaba_a into Inhibitory interneurons
var_vec = {'py_j_AMPA', 'in_j_AMPA', 'py_j_GABA_', 'py_j_GABAb_', 'in_j_GABA_'};

for ii = 1:length(var_vec)
    var_ = var_vec{ii};
    d = dir([folder '*' var_ '*']);
    tr_in_GABA_lif = zeros(1,length(d)-2);
    fr_in_GABA_lif = zeros(1,length(d)-2);
    fr_py_GABA_lif = zeros(1,length(d)-2);
    values_ = zeros(1,length(d)-2);

    % Sort files
    dates = [];
    for i = 1:length(d), dates(i) = d(i).datenum; end
    [~, idx] = sort(dates);
    d = d(idx);

    % Calculate excitability
    temp_values = 0:0.1:2;
    for i = 1:length(d)
        % Load file
        a = load([folder filesep d(i).name]);
        y_lif = a.LFP_V;
        t_lif = a.lfp_dt:a.lfp_dt:length(a.LFP_V)*a.lfp_dt;
        % Analyze excitabilty
        tr_in_GABA_lif(i) = analyze_excitability(y_lif', t_lif', 4899, -1.1, 2500, false);    
        fr_in_GABA_lif(i) = mean(a.R_in(1000:4000));    
        fr_py_GABA_lif(i) = mean(a.R_py(1000:4000));    
        try
            values_(i) = a.value_;
        catch ME
            if strcmp(ME.identifier, 'MATLAB:nonExistentField')
                values_(i) = temp_values(i);
            else
                rethrow(ME);
            end
        end
    end

    % Plot
    switch var_
        case 'py_j_AMPA'
            load parameter_sweeps/alpha_re;
            nmm_result = tr_alpha_re;
            nmm_in = fr_in_alpha_re;
            nmm_py = fr_py_alpha_re;
            title_ = 'Py_{AMPA}';
        case 'in_j_AMPA'
            load parameter_sweeps/alpha_e;
            nmm_result = tr_alpha_e;
            nmm_in = fr_in_alpha_e;
            nmm_py = fr_py_alpha_e;
            title_ = 'In_{AMPA}';
        case 'py_j_GABA_'
            load parameter_sweeps/alpha_i;
            nmm_result = tr_alpha_i;
            nmm_in = fr_in_alpha_i;
            nmm_py = fr_py_alpha_i;
            title_ = 'Py_{GABA_{A}}';
        case 'py_j_GABAb_'
            load parameter_sweeps/alpha_b;
            nmm_result = tr_alpha_b;
            nmm_in = fr_in_alpha_b;
            nmm_py = fr_py_alpha_b;
            title_ = 'Py_{GABA_{B}}';
        case 'in_j_GABA_'
            load parameter_sweeps/alpha_ri;
            nmm_result = tr_alpha_ri;
            nmm_in = fr_in_alpha_ri;
            nmm_py = fr_py_alpha_ri;
            title_ = 'In_{GABA_{A}}';
        otherwise
            error('Incorrect option');
    end

    % Excitability
    figure;
    plot(values_, tr_in_GABA_lif, 'k--');
    hold
    plot(values, nmm_result, 'k');
    l = legend({'LIF', 'NMM'});
    l.Location = 'best';
%     try
%         title(['Excitability vs ' a.pop_ '_' a.parameter]);
%     catch ME 
%         if strcmp(ME.identifier, 'MATLAB:nonExistentField')
            title(['Excitability vs' title_]);
%         else
%             rethrow(ME);
%         end
%     end
    ylabel('t_{recovery} (ms)');
    xlabel(title_);

    % Firing rate
    figure;
    l1 = plot(values_, fr_py_GABA_lif, '--'); hold
    l2 = plot(values_, fr_in_GABA_lif, '--');
    plot(values, nmm_py, 'Color', l1.Color);
    plot(values, nmm_in, 'Color', l2.Color);
    title(['Firing rates vs ' title_]);
    l = legend({'Pyramidal (LIF)', 'Inhibitory (LIF)', 'Pyramidal (NMM)', 'Inhibitory (NMM)'});
    l.Location = 'best';
    xlabel(title_);
    ylabel('Firing rate (Hz)');
end