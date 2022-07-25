%% ANALYZE_CORRELATION runs NMM and LIF for different u and calculates Spearman correlation.
% Compares the simulations of LIF (previously stored), and runs different
% NMM simulations matching the input spike rate. Compares the results with
% a Spearman correlation test.

function [u, R, P] = analyze_correlation()

    u = 1:20;%[0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9];%1:20;%
    R = nan(size(u));
    P = nan(size(u));
    for i = 1:length(u)
        data_file = ['C:/Users/artemios/Documents/Multiscale_Models_Data/spartan/lfp_u[' num2str(u(i)) '].mat'];
        [x, y, t, ~, ~] = NMM_diff_equations_DblExp_recursive('u', u(i));

        signal = 'vip'; % Options: 'vpi', 'vip', 'lfp'

        [x_nmm, x_lif, t_nmm, t_lif, v_pi, v_ip, u_lif, lfp_dt] = get_data(signal, x, y, t, data_file);

        [x_nmm, x_lif, t_nmm] = normalization(x_nmm, x_lif, t_nmm, t_lif);

        [R(i), P(i)] = do_spearmanCorr(x_nmm, x_lif);
    end % for loop
    
    do_plot(u,R,P);
    
end % function
%%
function [R,P] = do_spearmanCorr(x_nmm, x_lif)
    idxs = 1:min(length(x_nmm),length(x_lif)); % Adjust the index range so vectors are same size.
%     [R,P] = corrcoef(x_nmm(idxs), x_lif(idxs));
    [R,P] = corr(x_nmm(idxs)',x_lif(idxs)','Type','Spearman');
end

%%
function do_plot(u, R, P)
    if length(u) > 1
        figure
        RR = diag(R);
        mesh(u,u,RR, 'FaceColor', 'flat', 'EdgeColor', 'black')
        c = colorbar;
        colormap hot
        caxis([-0.2 0.2]);
        c.Limits = [-0.2 0.2];
    end
    
    figure
    stem(u, R);
end

%%
function [x_nmm, x_lif, t_nmm, t_lif, v_pi, v_ip, input_spike_rate, dt] = get_data(signal, x, y, t, data_file)
    %% NMM
    t_nmm = t;
    switch signal
        case 'vpi'
            x_nmm = x(:,1); % State 1 % NMM; % change dimm if forward model came from diff eqs instead of the nmm_toolbox
            x_nmm = x_nmm';
%             x_ = x(1,:);
        case 'vip'
            x_nmm = x(:,3); % State 3 % NMM
            x_nmm = x_nmm';
%             x_nmm = x(3,:);
        case 'lfp'
            x_nmm = y'; % NMM
%             x_nmm = y; % NMM
        otherwise
            error('Wrong options (signal)');
    end
    x_nmm = x_nmm(250:end);
    t_nmm = t_nmm(250:end);

    %% LIF
    load(data_file);
    
    trim = 2500; % Samples to remove from the beginning of the LFP_V vector
    LFP_V = LFP_V(trim:end);
    
    if exist('lfp_dt','var'), dt = lfp_dt; else, dt = 1e-4; end
    
    T = (trim + length(LFP_V)) * dt;
    t_lif = linspace(trim*dt,T,(T/dt)-trim);
    
    switch signal
        case 'vpi'
            if size(v_pi,1) > 1
                x_lif = mean(v_pi(:,trim:end),1); % State 3? % LIF
            else
                x_lif = v_pi(trim:end); % State 3? % LIF
            end
        case 'vip'
            if size(v_pi,1) > 1
                x_lif = mean(v_ip(:,trim:end),1); % State 1? % LIF
            else
                x_lif = v_ip(trim:end); % State 1? % LIF
            end
        case 'lfp'
            x_lif = LFP_V; % LIF
        otherwise
            error('Wrong options (signal)');
    end
end

%%
function [x_nmm, x_lif, t_nmm] = normalization(x_nmm, x_lif, t_nmm, t_lif)
    x_nmm = x_nmm - mean(x_nmm);
    x_nmm = x_nmm / rms(x_nmm);

    x_lif = x_lif - mean(x_lif);
    x_lif = x_lif / rms(x_lif);
    
    x_nmm = interp1(t_nmm, x_nmm, t_lif);
    t_nmm = t_lif(~isnan(x_nmm)); % Adjust t_nmm
    x_nmm = x_nmm(~isnan(x_nmm)); % Remove parsed NaN values
end