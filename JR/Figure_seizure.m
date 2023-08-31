% var_vec = {'no_drug', 'diazepam', 'baclofen', 'muscimol', 'muscimol_e'};
% var_vec = {'no_drug', 'diazepam', 'muscimol'};
% var_vec = {'muscimol', 'muscimol_e'};
% var_vec = {'muscimol'};
% var_vec = {'diazepam'};
% var_vec = {'baclofen'};
var_vec = {'no_drug'};

% Key values: Saturation (alpha_i = 1, alpha_e = 0.5), Normal (alpha_i = 1,
% alpha_e = 1), Oscillation (alpha_i = 1, alpha_e = 1.8) with HIGH_EXC as
% true. I.e., alpha_e = 0.5, then it is: 1.3 * 1.8 * default_parameter
%
% LIF folder: 'C://Users/artemios/Documents/Multiscale_Models_Data/'
% LIF Files: Oscillation: lfp_py__1_u0.mat
%            Saturation: lfp_py__2_u0.mat
%            Normal: lfp_py__0_u0.mat
%
% Run Figure_sample_recordings.m to plot an example of each response type


% NOTE: To run for different gains, use range_gains, to run for different
% input currents, use range_current (serch for text 'alt' to know where to
% change the code)

% range_gains = 0:0.05:2;
% range_gains = 0:0.1:2;
% range_gains = 0.4:0.3:4;
% range_gains = 0.4:0.1:4;
range_gains = 0.5:0.05:2;
% range_gains = 0.4:0.25:4;
% range_gains = 0:0.1:2;
% range_gains = 0:0.25:2;
% range_gains = 1;

% (alt)
% range_current = [0:50:500] * 1e-12;
% range_current = [0:5:250] * 1e-12;
% range_current = [-500:50:500] * 1e-12;
% range_current = [0:100:2000] * 1e-12;
range_current = 0;

HIGH_EXC = false;
NEW_PARAMS = false;
PLOT_LFP = false;
PLOT_FFT = false;

% close all
f = figure;
f.Position = [300 400 400*length(var_vec) 300];
f2 = figure;
f2.Position = [150 268 490*length(var_vec) 300];
f3 = figure;
f3.Position = [175 324 490*length(var_vec) 612];
cmap =[0 0 0; 0 0 1; 1 0 1; 1 0.8 0; 1 0 0];

clear state balance
for j = 1:length(var_vec)
    
    drug = var_vec{j}; 
    concentration = 0.5;
    
    switch drug
        case 'diazepam'
            % Values from paper:
            % tau_sri = (1 + 0.2) * tau_sri;
            % tau_sp = (1 + 0.2) * tau_sp;
            % a_i = (1 + 0.4) * a_i
            % a_ri = (1 + 0.34) * a_ri
            %
            a_i = 1 + 0.4 * concentration;
            a_ri = 1 + 0.34 * concentration;
            a_b = 1;
            t_sp = 1 + 0.2 * concentration;
            t_sri = 1 + 0.2 * concentration;
            a_u = 1;
            a_ui = 1;
            
        case 'muscimol'
            % Values from paper:
            % tau_sri = 0.87 * tau_sri;
            % tau_sp = 0.87 * tau_sp;
            % a_i = 0.19 * a_i
            % a_ri = 0.18 * a_ri
            %
            a_i = 1 - 0.81 * concentration;
            a_ri = 1 - 0.82 * concentration;
            a_b = 1;
            t_sp = 1 - 0.13 * concentration;
            t_sri = 1 - 0.13 * concentration;
            a_u = 1;
            a_ui = 1;
%             a_u = 0.43 * range;
%             a_ui = 0.5 * range;
            
        case 'baclofen'
            % Values from paper:
            % a_b = 0.34 * a_b
            %
            
            a_i = 1;
            a_ri = 1;
            a_b = 1 - 0.66 * concentration;
            t_sp = 1;
            t_sri = 1;
            a_u = 1;
            a_ui = 1;
            
        case 'muscimol_e'
            a_u = 1 - 0.57  * concentration;
            a_ui = 1 - 0.5 * concentration;
            a_i = 1;
            a_ri = 1;
            a_b = 1;
            t_sp = 1;
            t_sri = 1;
            
        otherwise
            a_u = 1;
            a_ui = 1;
            a_i = 1;
            a_ri = 1;
            a_b = 1;
            t_sp = 1;
            t_sri = 1;
    end
    
    if HIGH_EXC
        a_e = 1.3; %#ok<UNRCH>
        a_ri = 1.3 * a_ri;
        
        title_str = 'High excitability';
    elseif NEW_PARAMS
        a_b = 1.8;
        a_ri = 1.5;
        a_e = 1;
        title_str = 'New params';
    else
        a_e = 1; %#ok<UNRCH>
        
        title_str = 'Normal state';
    end
    
    state.(drug) = [];
    balance.(drug) = [];
    firing_rate_in.(drug) = [];
    firing_rate_py.(drug) = [];
    for i = 1:length(range_gains) % Varying inhibitory gain (alt)
        for ii = 1:length(range_gains) % Varying excitatory gain (alt)
%     for i = 1:length(range_current) % Varying inhibitory current (alt)
%         for ii = 1:length(range_current) % Varying excitatory current (alt)
            disp([drug ' | ' num2str(i) ' , ' num2str(ii)]);
            % Modify following line for gains or current variation (alt)
            [x, ~, t, f_e, f_i, params, y] = NMM_GABA('u', 0,...
                                            'alpha_e', a_e,... % x axis (alt)
                                            'alpha_i', range_gains(i) * a_i,... % y axis (alt)
                                            'alpha_re', 1,... 
                                            'alpha_ri', range_gains(ii) * a_ri,... 
                                            'alpha_b', a_b,...%range_gains(ii) * a_b,... 
                                            'alpha_u', a_u,... 
                                            'alpha_uinterneuron', a_ui,...
                                            'tau_sp', t_sp,... 
                                            'tau_sri', t_sri,...
                                            'CURRENT_E', 0,... %range_current(ii),... % (alt)
                                            'CURRENT_I', 0,... %range_current(i),... % (alt)
                                            ...%'u_bkg', 0,... % (alt)
                                            'CURRENT_TIME', 1:2000 );
            % Calculate fft to estimate oscillatory activity
            if PLOT_FFT, fig_101 = figure(101); cla; else, fig_101 = []; end
            [~, X_, F_] = fft_plot(y(500:end)-mean(y(500:end)), t(500:end),fig_101,PLOT_FFT);
            
            if (mean(f_i(500:end)) > 65) && (mean(f_e(500:end)) > 35)
                state.(drug)(i,ii) = 3; % Saturation
                
            elseif max(X_) > 1.5e-3 % 2.5e9
                % Oscillations, papers on Gamma band Epilepsy and Low freq:
                % + Medvedev et al., 2000, Kainic acid induces distinct 
                %   types of epileptiform discharge with differential 
                %   involvement of hippocampus and neocortex
                %
                % + Herrmann et al., 2005, Human EEG gamma oscillations in 
                %   neuropsychiatric disorders
                %
                %state.(drug)(i,ii) = 1; % Oscillation
                [~,indice] = max(X_);
                if F_(indice) < 13
                    % Alpha-ish, sleep, low freq
                    state.(drug)(i,ii) = 1; % Oscillation
                elseif F_(indice) < 200
                    state.(drug)(i,ii) = 2; % Oscillation
                    % Gamma-ish, fast oscillations ~60 Hz
                else
                    state.(drug)(i,ii) = 0; % Not really oscillation, look closely to the LFP it is likely varying very fast with a very low variance
                end
            elseif (mean(f_i(500:end)) == 0) || (mean(f_e(500:end)) == 0)
                state.(drug)(i,ii) = -1; % Low state
            else
                state.(drug)(i,ii) = 0; % Normal
            end
            
            if PLOT_LFP
                figure(100); %#ok<UNRCH>
                
                subplot(211)
                cla;
                plot(t, 1e9*y, 'Color', cmap(state.(drug)(i,ii)+2 ,:));
                ylabel('LFP (nA)');
                
                subplot(212)
                cla;
                plot(t, f_e, 'Color', cmap(state.(drug)(i,ii)+2 ,:)); hold on;
                plot(t, f_i, 'Color', cmap(state.(drug)(i,ii)+2 ,:)); hold off;
                xlabel('Time (s)');
                ylabel('Firing rates (Hz)');
                legend({'Pyramidal' 'Interneurons'});
                
                drawnow
            end
            
            %Balance
            milivolts_scale = 1e-3;
            nanoamps_scale = 1e9;
            nmm_i_pi = (x(:,1) + x(:,9)) * milivolts_scale * params.g_m_P;
            nmm_i_pe = (x(:,5) + x(:,7)) * milivolts_scale * params.g_m_P;
            nmm_i_ie = (x(:,3) + x(:,19)) * milivolts_scale * params.g_m_I;
            nmm_i_ii = (x(:,11)) * milivolts_scale * params.g_m_I;
            L = numel(nmm_i_pe);
            balance.(drug)(i,ii) = -(mean(nmm_i_pi(round(L/2):end)) + mean(nmm_i_pe(round(L/2):end)))*nanoamps_scale;
            
            % Firing rate
            firing_rate_in.(drug)(i,ii) = mean(f_i(500:end));
            firing_rate_py.(drug)(i,ii) = mean(f_e(500:end));
                
        end
    end
    
    %%
    figure(f)
    ax = subplot(1, length(var_vec), j);
    colormap(cmap);    
%     imagesc(range_gains, range_gains, state.(drug));xlabel('Excitatory gain');ylabel('Inhibitory gain'); % (alt)
    imagesc(range_gains, range_gains, state.(drug)); xlabel('Recurrent inhibitory gain');ylabel('Inhibitory gain'); % (alt)
%     imagesc(range_current*1e12, range_current*1e12, state.(drug));xlabel('Py input current (pA)');ylabel('In input current (pA)'); % (alt)
    if (concentration ~= 0) && ~strcmp('no_drug', drug)
        title([drug ' | concentration: ' num2str(concentration)]);
    else 
        title('NMM');
    end
    caxis([-1 3]);
    ax.FontSize = 12;
    ax.View = ([0 -90]);
%     ax.XTick = [0.4:1.2:4];
%     ax.YTick = [0.4:1.2:4];
    drawnow
    
    figure(f2)
    ax = subplot(1, length(var_vec), j);
%     colormap(flipud(jet));    
    load('custom_colormap_parula.mat')
    colormap(parula_custom);
%     imagesc(range_gains, range_gains, balance.(drug));xlabel('Excitatory gain');ylabel('Inhibitory gain'); % (alt)
    imagesc(range_gains, range_gains, balance.(drug)); xlabel('Recurrent inhibitory gain');ylabel('Inhibitory gain'); % (alt)
%     imagesc(range_current*1e12, range_current*1e12, balance.(drug));xlabel('Py input current (pA)');ylabel('In input current (pA)'); % (alt)
%     caxis([-0.25 0.25]);
    caxis([-0.5 0.5]);
    if (concentration ~= 0) && ~strcmp('no_drug', drug)
        title([drug ' | concentration: ' num2str(concentration)]);
    else 
        title('NMM');
    end
    ax.FontSize = 12;
    ax.View = ([0 -90]);
%     ax.XTick = [0.4:1.2:4];
%     ax.YTick = [0.4:1.2:4];
    c = colorbar;
    c.Label.String = 'Input current balance (nA)';
    c.Label.FontSize = 12;
    drawnow
    
    
    figure(f3)
    % Inhibitory firing rate
    ax = subplot(2, length(var_vec), j);
%     colormap(flipud(jet));    
    colormap(jet)
%     imagesc(range_gains, range_gains, firing_rate_in.(drug));xlabel('Excitatory gain');ylabel('Inhibitory gain'); % (alt)
    imagesc(range_gains, range_gains, firing_rate_in.(drug)); xlabel('Recurrent inhibitory gain');ylabel('Inhibitory gain'); % (alt)
%     imagesc(range_current*1e12, range_current*1e12, firing_rate_in.(drug));xlabel('Py input current (pA)');ylabel('In input current (pA)'); % (alt)
    max_fun = @(a,b) max(3*median(a,'all'), b); % maximum between 3 times the median of the LIF and the first nmm's lower than the maximum firing rate    
    maximum_value = max_fun(firing_rate_in.(drug), 0);
%     caxis([0 30]);
    caxis([0 maximum_value]);
    if (concentration ~= 0) && ~strcmp('no_drug', drug)
        title([drug ' | concentration: ' num2str(concentration)]);
    else 
        title('Interneurons');
    end
    ax.FontSize = 12;
    ax.View = ([0 -90]);
%     ax.XTick = [0.4:1.2:4];
%     ax.YTick = [0.4:1.2:4];
    c = colorbar;
    c.Label.String = 'Interneurons mean firing rate';
    c.Label.FontSize = 12;
    drawnow
    
    % Excitatory firing rate
    ax = subplot(2, length(var_vec), length(var_vec)+j);    
%     colormap(flipud(jet));    
    colormap(jet)
%     imagesc(range_gains, range_gains, firing_rate_py.(drug));xlabel('Excitatory gain');ylabel('Inhibitory gain'); % (alt)
    imagesc(range_gains, range_gains, firing_rate_py.(drug)); xlabel('Recurrent inhibitory gain');ylabel('Inhibitory gain'); % (alt)
%     imagesc(range_current*1e12, range_current*1e12, firing_rate_py.(drug));xlabel('Py input current (pA)');ylabel('In input current (pA)'); % (alt)
    maximum_value = max_fun(firing_rate_py.(drug), 0);
%     caxis([0 0.4]);
%     caxis([0 15]);
    caxis([0 maximum_value]);
    if (concentration ~= 0) && ~strcmp('no_drug', drug)
        title([drug ' | concentration: ' num2str(concentration)]);
    else 
        title('Pyramidal');
    end
    ax.FontSize = 12;
    ax.View = ([0 -90]);
%     ax.XTick = [0.4:1.2:4];
%     ax.YTick = [0.4:1.2:4];
    c = colorbar;
    c.Label.String = 'Pyramidal mean firing rate';
    c.Label.FontSize = 12;
    drawnow
    
end
