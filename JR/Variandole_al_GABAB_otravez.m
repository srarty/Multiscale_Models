%% Runs various GABA agonists on NMM_GABAB.m
% Crear una gráfica como la de excitation de Maths in the Brain para las
% tres drogas diferentes. Poner una raya vertical en el valor reportado en
% el paper Wang et al 2019
%
% For this script to make any sense, NMM_GABA needs to be ran with the
% injected current different to zero (params.options.CURRENT > 0)
%
% Artemio - May 2023

var_vec = {'diazepam', 'baclofen', 'muscimol'}; %, 'muscimol_e'};
% var_vec = {'muscimol', 'muscimol_e'};
% var_vec = {'diazepam'};
% var_vec = {'baclofen'};

HIGH_EXC = true;

% close all
figure
hold

recovery = struct;%size(range);
for j = 1:length(var_vec)
    
    drug = var_vec{j}; 
    range = 0:0.05:1;
%     range = 0.3;
    
    switch drug
        case 'diazepam'
            % Values from paper:
            % tau_sri = (1 + 0.2) * tau_sri;
            % tau_sp = (1 + 0.2) * tau_sp;
            % a_i = (1 + 0.4) * a_i
            % a_ri = (1 + 0.34) * a_ri
            %
            a_i = 1 + 0.4 * range;
            a_ri = 1 + 0.34 * range;
            a_b = 1 * ones(size(range));
            t_sp = 1 + 0.2 * range;
            t_sri = 1 + 0.2 * range;
            a_u = 1 * ones(size(range));
            a_ui = 1 * ones(size(range));
            
        case 'muscimol'
            % Values from paper:
            % tau_sri = 0.87 * tau_sri;
            % tau_sp = 0.87 * tau_sp;
            % a_i = 0.19 * a_i
            % a_ri = 0.18 * a_ri
            %
            a_i = 1 - 0.81 * range;
            a_ri = 1 - 0.82 * range;
            a_b = 1 * ones(size(range));
            t_sp = 1 - 0.13 * range;
            t_sri = 1 - 0.13 * range;
            a_u = 1 * ones(size(range));
            a_ui = 1 * ones(size(range));
%             a_u = 0.43 * range;
%             a_ui = 0.5 * range;
            
        case 'baclofen'
            % Values from paper:
            % a_b = 0.34 * a_b
            %
            a_i = 1 * ones(size(range));
            a_ri = 1 * ones(size(range));
            a_b = 1 - 0.66 * range;
            t_sp = 1 * ones(size(range));
            t_sri = 1 * ones(size(range));
            a_u = 1 * ones(size(range));
            a_ui = 1 * ones(size(range));
            
        case 'muscimol_e'
            a_u = 1 - 0.57  * range;
            a_ui = 1 - 0.5 * range;
            a_i = 1 * ones(size(range));
            a_ri = 1 * ones(size(range));
            a_b = 1 * ones(size(range));
            t_sp = 1 * ones(size(range));
            t_sri = 1 * ones(size(range));
            
        otherwise
            error('Invalid option');
    end
    
    if HIGH_EXC
        a_e = 1.3 * ones(size(range));
        a_ri = 1.3 * a_ri;
        
        title_str = 'High excitability';
    else
        a_e = ones(size(range));
        
        title_str = 'Normal state';
    end
    
    recovery.(drug) = [];
    state.(drug) = [];
    for i = 1:length(range)
        [~, ~, t, f_e, f_i, ~, y] = NMM_GABA('u', 0, 'alpha_e', a_e(i), 'alpha_i', a_i(i), 'alpha_ri', a_ri(i), 'alpha_b', a_b(i), 'alpha_u', a_u(i), 'alpha_uinterneuron', a_ui(i), 'tau_sp', t_sp(i), 'tau_sri', t_sri(i));
        
        % Calculate fft to estimate oscillatory activity
        [~, X_] = fft_plot(y(1000:end)-mean(y(1000:end)), t(1000:end),[],false);
        if (mean(f_i(500:end)) > 60) || (mean(f_e(500:end)) > 35)
            state.(drug)(i) = 2; % Saturation
        elseif max(X_) > 1e8
            state.(drug)(i) = 1; % Oscillation
        else
            state.(drug)(i) = 0; % Normal
        end
        
        [~, ~, t_rec, ~, ~, ~, y_rec] = NMM_GABA('CURRENT', 50e-12, 'u', 0, 'alpha_e', a_e(i), 'alpha_i', a_i(i), 'alpha_ri', a_ri(i), 'alpha_b', a_b(i), 'alpha_u', a_u(i), 'alpha_uinterneuron', a_ui(i), 'tau_sp', t_sp(i), 'tau_sri', t_sri(i));
        recovery.(drug)(i) = analyze_excitability(y_rec,t_rec, 1489, -5, 1000, false);
        
    end
    
%     figure;
    plot(range, 100 * recovery.(drug) / recovery.(drug)(1), '-o'); ylabel('t_{r} (%)');
%     plot(range(~state.(drug)), 1e3*recovery.(drug)(~state.(drug)), '-o'); ylabel('t_{r} (ms)');
    xlabel('Concentration (au)');

    %% Store values
    %{

    results = struct;
    results.value = value;
    results.range = range;
    results.value2 = value2;
    results.range2 = range2;
    results.freqs = freqs;
    results.freqs_py = freqs_py;
    results.freqs_in = freqs_in;

    folder = 'C:\Users\artemios\Documents\Multiscale_Models_Data\Spartan\firing_rates_nmm\';
    name = [value ' vs ' value2];
    if isempty(dir([folder name]))
        save([folder name '.mat'], 'results');
        disp('results saved');
    else
        error(['Results not saved, file exists | j = ' num2str(j)]);
    end
    %}
end

% plot(range([1 end]), 1e3 * recovery.(drug)(1) * ones([1 2]),'--', 'Color', [0.7 0.7 0.7]);
l = legend(var_vec);
l.Location = 'best';
title(title_str);
ylim(1e3*[0.1 5]);
