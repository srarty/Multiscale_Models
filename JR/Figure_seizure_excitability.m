% Test different drugs and plot colormaps of their excitability and oscillatory profile.
% 
% Colormaps vary alpha_i and alpha_e. To test the default parameters,
% run var_vec = 'no_drug' and HIGH_EXC = false.
%
% Artemio - April 2023

% var_vec = {'no_drug', 'diazepam', 'baclofen', 'muscimol', 'muscimol_e'};
% var_vec = {'no_drug', 'diazepam', 'muscimol'};
% var_vec = {'muscimol'};
% var_vec = {'diazepam'};
var_vec = {'no_drug'};

range_gains = 0:0.05:2;
% range_gains = 0:0.25:2;

HIGH_EXC = false;
PLOT_LFP = false;
PLOT_FFT = false;

% close all
f = figure;
f.Position = [300 400 400*length(var_vec) 300];
% Define colormap
% load custom_colormap
% cmap = custom_colormap;
cmap = colormap(jet(512)); % use flipud(jet) to invert the colors in the colormap, jet can be any colormap
cmap_tri =[0 0 1; 1 1 0; 1 0 0]; % Tri-color colourmap for oscillation check
colormap(cmap); % apply custom colormap
    
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
    else
        a_e = 1; %#ok<UNRCH>
        
        title_str = 'Normal state';
    end
    
    recovery.(drug) = [];
    state.(drug) = [];
    for i = 1:length(range_gains) % Varying inhibitory gain
        for ii = 1:length(range_gains) % Varying excitatory gain
            disp([drug ' | ' num2str(i) ' , ' num2str(ii)]);
            [x, ~, t, f_e, f_i, ~, y] = NMM_GABA('CURRENT', 50e-12, ...
                                            'u', 0,...
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
            [~, X_] = fft_plot(y(1600:end)-mean(y(1600:end)), t(1600:end),[],PLOT_FFT);
            
            if (mean(f_i(500:end)) > 60) || (mean(f_e(500:end)) > 35)
                state.(drug)(i,ii) = 2; % Saturation
                
            elseif max(X_) > 1e-3 % 2.5e9
                state.(drug)(i,ii) = 1; % Oscillation
                
            else
                state.(drug)(i,ii) = 0; % Normal
            end
            
            % Only measure recovery time if there wasn't saturation nor
            % oscillations:
            if ~state.(drug)(i,ii)
                recovery.(drug)(i,ii) = analyze_excitability(y, t, 1489, -3, 1000, false);
            else
                recovery.(drug)(i,ii) = 6;
            end

            if PLOT_LFP
                figure(100);
                
                subplot(211)
                cla;
                plot(t, 1e9*y, 'Color', cmap_tri(state.(drug)(i,ii)+1 ,:));
                ylabel('LFP (nA)');
                
                subplot(212)
                cla;
                plot(t, f_e, 'Color', cmap_tri(state.(drug)(i,ii)+1 ,:)); hold on;
                plot(t, f_i, 'Color', cmap_tri(state.(drug)(i,ii)+1 ,:)); hold off;
                xlabel('Time (s)');
                ylabel('Firing rates (Hz)');
                legend({'Pyramidal' 'Interneurons'});
            end
%             response = questdlg('What happened?', 'Manual input', 'Normal', 'Saturation', 'Oscillation', 'Normal');
%             switch response
%                 case 'Normal'
%                     state.(drug)(i,ii) = 0;
%                 case 'Saturation'
%                     state.(drug)(i,ii) = 1;
%                 case 'Oscillation'
%                     state.(drug)(i,ii) = 2;
%                 otherwise
%                     disp('Operation terminated by manual input');
%                     return;
%             end
                
        end
    end
    
    figure(f)
    ax = subplot(1, length(var_vec), j);
    colormap(cmap);    
    imagesc(range_gains, range_gains, recovery.(drug));
    caxis([0 0.5]);
    ax.View = ([0 -90]);
    xlabel('Excitatory gain');
    ylabel('Inhibitory gain');
    if (concentration ~= 0) && ~strcmp('no_drug', drug)
        title([drug ' | concentration: ' num2str(concentration)]);
    else 
        title('No drugs');
    end
    
end