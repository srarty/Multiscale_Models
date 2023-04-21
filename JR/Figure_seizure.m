% var_vec = {'no_drug', 'diazepam', 'baclofen', 'muscimol', 'muscimol_e'};
var_vec = {'no_drug', 'diazepam', 'muscimol'};
% var_vec = {'muscimol', 'muscimol_e'};
% var_vec = {'muscimol'};
% var_vec = {'no_drug'};

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
range_gains = 0:0.05:2;
% range_gains = 0:0.25:2;
% range_gains = [0.5 1 1.8];
HIGH_EXC = false;
PLOT_LFP = false;
PLOT_FFT = false;

% close all
f = figure;
f.Position = [300 400 400*length(var_vec) 300];
f2 = figure;
cmap =[0 0 1; 1 1 0; 1 0 0];

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
        a_e = 1.3;
        a_ri = 1.3 * a_ri;
        
        title_str = 'High excitability';
    else
        a_e = 1; %#ok<UNRCH>
        
        title_str = 'Normal state';
    end
    
    state.(drug) = [];
    balance.(drug) = [];
    for i = 1:length(range_gains) % Varying inhibitory gain
        for ii = 1:length(range_gains) % Varying excitatory gain
            disp([drug ' | ' num2str(i) ' , ' num2str(ii)]);
            [x, ~, t, f_e, f_i, params, y] = NMM_GABA('u', 0,...
                                            'alpha_e', range_gains(ii)*a_e,... 
                                            'alpha_i', range_gains(i)*a_i,... 
                                            'alpha_re', 1,... 
                                            'alpha_ri', a_ri,... 
                                            'alpha_b', a_b,... 
                                            'alpha_u', a_u,... 
                                            'alpha_uinterneuron', a_ui,... 
                                            'tau_sp', t_sp,... 
                                            'tau_sri', t_sri);
            % Calculate fft to estimate oscillatory activity
            [~, X_] = fft_plot(y(1000:end)-mean(y(1000:end)), t(1000:end),[],PLOT_FFT);
            
            if (mean(f_i(500:end)) > 60) || (mean(f_e(500:end)) > 35)
                state.(drug)(i,ii) = 2; % Saturation
                
            elseif max(X_) > 1e-3 % 2.5e9
                state.(drug)(i,ii) = 1; % Oscillation
                
            else
                state.(drug)(i,ii) = 0; % Normal
            end
            
            if PLOT_LFP
                figure(100);
                
                subplot(211)
                cla;
                plot(t, 1e9*y, 'Color', cmap(state.(drug)(i,ii)+1 ,:));
                ylabel('LFP (nA)');
                
                subplot(212)
                cla;
                plot(t, f_e, 'Color', cmap(state.(drug)(i,ii)+1 ,:)); hold on;
                plot(t, f_i, 'Color', cmap(state.(drug)(i,ii)+1 ,:)); hold off;
                xlabel('Time (s)');
                ylabel('Firing rates (Hz)');
                legend({'Pyramidal' 'Interneurons'});
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
                
        end
    end
    
    figure(f)
    ax = subplot(1, length(var_vec), j);
    colormap(cmap);    
    imagesc(range_gains, range_gains, state.(drug));
    caxis([0 2]);
    ax.View = ([0 -90]);
    xlabel('Excitatory gain');
    ylabel('Inhibitory gain');
    if (concentration ~= 0) && ~strcmp('no_drug', drug)
        title([drug ' | concentration: ' num2str(concentration)]);
    else 
        title('No drugs');
    end
    drawnow
    
    figure(f2)
    ax = subplot(1, length(var_vec), j);
    colormap(jet);    
    imagesc(range_gains, range_gains, balance.(drug));
%     caxis([0 2]);
    ax.View = ([0 -90]);
    xlabel('Excitatory gain');
    ylabel('Inhibitory gain');
    if (concentration ~= 0) && ~strcmp('no_drug', drug)
        title([drug ' | concentration: ' num2str(concentration)]);
    else 
        title('No drugs');
    end
    drawnow
    
end