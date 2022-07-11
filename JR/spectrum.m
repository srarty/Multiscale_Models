%% Computes the FFT and Spectrogram of the states and output of each model
%
% Need to run the NMM before running spectrum.m. This is to generate the
% NMM data. The LIF data needs to be saved in the location specified by
% data_file. 
% 
% Run the NMM with either of these files:
%      + NMM_diff_equations.m
%      + NMM_diff_equations_Dbl_Exp.m
%      + example_nmm.m
%
% Artemio - March 2022

% NOTES:
% NMM reduced frequency reponse with increased u, LIF other way around
function varargout = spectrum(x, y, t, varargin)
    global PLOT

%     close all
%     data_file = 'C:/Users/artemios/Documents/Multiscale_Models_Data/lfp_3.mat';
%       data_file = 'C:/Users/artemios/Documents/Multiscale_Models_Data/spartan/lfp_96.mat'; % u=9, I_PP = 0.415
%     data_file = 'C:/Users/artemios/Documents/Multiscale_Models_Data/spartan/lfp_126.mat'; % u=20, I_PP=0.19
%     data_file = 'C:/Users/artemios/Documents/Multiscale_Models_Data/spartan/lfp_127.mat'; % u=20, I_PP=0.125

%     data_file = 'C:/Users/artemios/Documents/Multiscale_Models_Data/spartan/lfp_261.mat'; % u=9, I_PP=0.415
%     data_file = 'C:/Users/artemios/Documents/Multiscale_Models_Data/spartan/lfp_262.mat'; % u=14, I_PP=0.415
%     data_file = 'C:/Users/artemios/Documents/Multiscale_Models_Data/spartan/lfp_263.mat'; % u=20, I_PP=0.415

    data_file = 'C:/Users/artemios/Documents/Multiscale_Models_Data/spartan/lfp_267.mat'; % u=9, I_PP=0.415, 3 seconds
%     data_file = 'C:/Users/artemios/Documents/Multiscale_Models_Data/spartan/lfp_268.mat'; % u=14, I_PP=0.415, 3 seconds
%     data_file = 'C:/Users/artemios/Documents/Multiscale_Models_Data/spartan/lfp_269.mat'; % u=20, I_PP=0.415, 3 seconds
    

    if nargin > 3, PLOT = varargin{1}; else, PLOT = true; end
    if nargin > 4, data_file = varargin{2}; end

    signal = 'lfp'; % vpi, vip, lfp
    
    [x_nmm, x_lif, t_nmm, t_lif, v_pi, v_ip, u_lif, lfp_dt] = get_data(signal, x, y, t, data_file);
    
    [x_nmm, x_lif, t_nmm] = normalization(x_nmm, x_lif, t_nmm, t_lif);
    
    if PLOT
        figure(1);clf;
        plot(t_nmm, x_nmm, 'k');
        hold;
        plot(t_lif, x_lif, '--k');
        legend({'NMM', 'LIF'});
        xlabel('Time (s)');
        title(['Comparison of ' signal]);
    end
    
    
    harmonic_nmm = do_plot('nmm', signal, x_nmm, t_nmm);
    varargout = {harmonic_nmm};
    harmonic_lif = do_plot('lif', signal, x_lif, t_lif);
    varargout = {harmonic_lif};
%     varargout{end + 1} = u_lif;
    
    if PLOT
        figure(2); clf;
        subplot(2,1,1);
        plot(t,x(:,1)); hold on;
        plot(t,x(:,3));
        
        plot_lif_results(v_pi, v_ip, lfp_dt); % Won't run if lif results haven't been loaded
    end
     
%     do_correlation(x_nmm, x_lif, t_nmm, t_lif);
%     do_variance(x_nmm, x_lif, t_nmm, t_lif);
    
%     varargout = {x_nmm, x_lif, t_nmm, t_lif};
end

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
        
    LFP_V = LFP_V(2500:end);
    
    if exist('lfp_dt','var'), dt = lfp_dt; else, dt = 1e-4; end
    
    T = length(LFP_V) * dt;
    t_lif = linspace(0,T,T/dt);
    
    switch signal
        case 'vpi'
            x_lif = mean(v_pi,1); % State 3? % LIF
        case 'vip'
            x_lif = mean(v_ip,1); % State 1? % LIF
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

%%
function oscillation = do_plot(model, signal, x_, t)
    global PLOT

    dt = t(2) - t(1);

    switch signal
        case 'vpi'
            ystr = 'V_{pi}';
        case 'vip'
            ystr = 'V_{ip}';
        case 'lfp'
            ystr = 'V_m (Py)';
        otherwise
            error('Wrong options (signal)');
    end
        
	w = round(length(t)/50); % Window size
    so = round(w*0.9); % Samples overlap
    freqbins = 10 * w; %Evaluate the spectrum at (128/2)+1=65 frequencies and (length(x)−120)/(128−120)=235 time bins.    
    freqrange = [0 0.2]; % xlim values

%     disp(['freq bins = ', num2str( (freqbins/2)+1 )]);
%     disp(['time bins = ', num2str( (length(x_)-so)/(freqbins-so) )]);
            
    %% FFT

    Fs = 1/dt;
    L = length(t);
    n = 2^nextpow2(L);
    X = fft(x_,n);
    P2 = abs(X/L);
    P1 = P2(:,1:n/2+1);
    P1(:,2:end-1)=2*P1(:,2:end-1);
    
    if PLOT
        f = figure(100);
        f.Position([1 2 3 4]) = [10 64 1230 840];

        if strcmp('lif', model)
            subplot(2,2,2); 
            titlestr = 'LIF';
        else
            subplot(2,2,1); 
            titlestr = 'NMM';
        end

        plot(0:(Fs/n):(Fs/2-Fs/n), P1(1:n/2));
        title(titlestr);
        xlabel('Frequency (Hz)');
        ylabel(['Power [', ystr, ']']);
        xlim([0 200]);
    end
    % Calculate harmonic:
    freqs = 0:(Fs/n):(Fs/2-Fs/n);
    amps = P1(1:n/2);
    [~,idx] = max(amps);
    oscillation = freqs(idx);
        
    if PLOT
        %% Spectrogram
%         f = figure(101);clf;
%         f.Position([3 4]) = [1230 420];
        if strcmp('lif', model), subplot(2,2,4); else, subplot(2,2,3); end
        spectrogram(x_, w, so, freqbins, Fs, 'yaxis','power');
        [~,~,~,ps] = spectrogram(x_, w, so, freqbins, Fs, 'yaxis','power');

        title([titlestr, ' (', ystr, ')']);
        ax = gca;
    %     ax.YScale = 'log';
        ylim(freqrange);
    end
    
end

%%
function plot_lif_results(v_pi, v_ip, dt)
    t = 0:dt:(size(v_pi,2)-1) * dt;

    x1 = mean(v_pi,1);
    x1_ = [0 mean(diff(v_pi,[],2))];
    x3 = mean(v_ip,1);
    x3_ = [0 mean(diff(v_ip,[],2))];
    
    figure; 
    subplot(1,2,1)
    plot(x1*1e3, x3*1e3);
    xlabel('x1');
    ylabel('x3');

    subplot(1,2,2)
    plot(x1*1e3,x1_*1e6)
    xlabel('v');
    ylabel('z');
    
    figure(2);
    subplot(2,1,2)
    plot(t, x1*1e3); hold on;
    plot(t, x3*1e3);
    legend({'x1', 'x3'});
    xlabel('Time (s)');
    ylabel('V?');
    
%     figure; plot3(t * 1e3,x3*1e3,x3_*1e6, 'r')
%     xlabel('t');
%     ylabel('v');
%     zlabel('z');
end

function do_correlation(x_nmm, x_lif, t_nmm, t_lif)    
    [correlation, lags] = xcorr(x_nmm, x_lif);
    auto_nmm = autocorr(x_nmm, 'NumLags', 1000);
    auto_lif = autocorr(x_lif, 'NumLags', 1000);
        
    figure;
    [correlation, lags] = xcorr(x_nmm, x_nmm);
    plot(lags, correlation);
    hold
    [correlation, lags] = xcorr(x_lif, x_lif);
    plot(lags, correlation);
    
    figure;
    a = bar(auto_nmm, 'BarWidth', 0.85); hold
    bar(auto_lif, 'BarWidth', 0.75);
    legend({'NMM', 'LIF'});
end


function do_variance(x_nmm, x_lif, t_nmm, t_lif)    

    window_size = 0.1; % Size of window in seconds
    w_nmm = find(t_nmm <= window_size,1,'last');
    w_lif = find(t_lif <= window_size,1,'last');
    
    if isempty(w_nmm) || isempty(w_lif), error('The window size is too small (smaller than the sampling rate)');end
    
    variance_nmm = zeros(1,length(x_nmm) - w_nmm);
    variance_lif = zeros(1,length(x_lif) - w_lif);
    
    
    for i = 1:length(x_nmm) - w_nmm
        variance_nmm(i) = var(x_nmm(i : i + w_nmm));
    end
    
    for i = 1:length(x_lif) - w_lif
        variance_lif(i) = var(x_lif(i : i + w_lif));
    end
    
    figure; 
    plot(variance_nmm);
    hold on
    plot(variance_lif);    
end
