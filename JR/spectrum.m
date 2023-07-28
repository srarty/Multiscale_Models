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
function varargout = spectrum(x, yy, t, varargin)
    global PLOT

%     close all
%     data_file = 'C:/Users/artemios/Documents/Multiscale_Models_Data/lfp_3.mat';
%       data_file = 'C:/Users/artemios/Documents/Multiscale_Models_Data/spartan/lfp_96.mat'; % u=9, P_II = 0.451
%     data_file = 'C:/Users/artemios/Documents/Multiscale_Models_Data/spartan/lfp_126.mat'; % u=20, P_II=0.19
%     data_file = 'C:/Users/artemios/Documents/Multiscale_Models_Data/spartan/lfp_127.mat'; % u=20, P_PP=0.125

%     data_file = 'C:/Users/artemios/Documents/Multiscale_Models_Data/spartan/lfp_261.mat'; % u=9, P_II=0.451
%     data_file = 'C:/Users/artemios/Documents/Multiscale_Models_Data/spartan/lfp_262.mat'; % u=14, P_II=0.451
%     data_file = 'C:/Users/artemios/Documents/Multiscale_Models_Data/spartan/lfp_263.mat'; % u=20, P_II=0.451

%     data_file = 'C:/Users/artemios/Documents/Multiscale_Models_Data/spartan/lfp_267.mat'; % u=9, P_II=0.451, 3 seconds
%     data_file = 'C:/Users/artemios/Documents/Multiscale_Models_Data/spartan/lfp_268.mat'; % u=14, P_II=0.451, 3 seconds
%     data_file = 'C:/Users/artemios/Documents/Multiscale_Models_Data/spartan/lfp_269.mat'; % u=20, P_II=0.451, 3 seconds
    
%     data_file = 'C:/Users/artemios/Documents/Multiscale_Models_Data/spartan/lfp_282.mat'; % u = [0 1 0 10];
%     data_file = 'C:/Users/artemios/Documents/Multiscale_Models_Data/spartan/lfp_283.mat'; % u=0
%     data_file = 'C:/Users/artemios/Documents/Multiscale_Models_Data/spartan/lfp_279.mat'; % u=[0, 0.2, 0, 0.4]
%     data_file = 'C:/Users/artemios/Documents/Multiscale_Models_Data/spartan/lfp_286.mat'; % u=10
%     data_file = 'C:/Users/artemios/Documents/Multiscale_Models_Data/spartan/lfp_287.mat'; % u=15
%     data_file = 'C:/Users/artemios/Documents/Multiscale_Models_Data/spartan/lfp_288.mat'; % u = [0 1 0 15];

%     data_file = 'C:/Users/artemios/Documents/Multiscale_Models_Data/spartan/lfp_u[15].mat';
%     data_file = 'C:/Users/artemios/Documents/Multiscale_Models_Data/spartan/lfp_313.mat';

%     data_file = 'C:/Users/artemios/Documents/Multiscale_Models_Data/lfp_17.mat';
%     data_file = 'C:/Users/artemios/Documents/Multiscale_Models_Data/lfp_24.mat';
%     data_file = 'C:/Users/artemios/Documents/Multiscale_Models_Data/lfp_25.mat'; % P[P->I] = 0.2
%     data_file = 'C:/Users/artemios/Documents/Multiscale_Models_Data/lfp_40.mat';
%     data_file = 'C:/Users/artemios/Documents/Multiscale_Models_Data/lfp_45.mat';

%     data_file = 'C:/Users/artemios/Documents/Multiscale_Models_Data/lfp_61.mat'; % u=[0 0.25 0.5 1]
%     data_file = 'C:/Users/artemios/Documents/Multiscale_Models_Data/lfp_63.mat'; % step response (normal parametesr)
%     data_file = 'C:/Users/artemios/Documents/Multiscale_Models_Data/lfp_64.mat'; % step response (j_pi = 18)
% data_file = 'C:/Users/artemios/Documents/Multiscale_Models_Data/lfp_75.mat'; % impulse response (50 pA), j_pi = 37, alpha_i = -0.5xxx
% data_file = 'C:/Users/artemios/Documents/Multiscale_Models_Data/lfp_80.mat'; % impulse response (100 pA), j_pi = 37, alpha_i = -0.5xxx
% data_file = 'C:/Users/artemios/Documents/Multiscale_Models_Data/lfp_79.mat'; % impulse response (500 pA), j_pi = 37, alpha_i = -0.5xxx
% data_file = 'C:/Users/artemios/Documents/Multiscale_Models_Data/lfp_65.mat'; % impulse response (500 pA), j_pi = 21.0666, alpha_i = -0.3
% data_file = 'C:/Users/artemios/Documents/Multiscale_Models_Data/lfp_66.mat'; % Seizure: j_pi = 21.0666, alpha_i = -0.3 (random LIF th and t_ref)

% data_file = 'C:/Users/artemios/Documents/Multiscale_Models_Data/lfp_76.mat'; % GABAb | u = 5
% data_file = 'C:/Users/artemios/Documents/Multiscale_Models_Data/lfp_86.mat'; % GABAb impulse response | 500pA
% data_file = 'C:/Users/artemios/Documents/Multiscale_Models_Data/lfp_89.mat'; % GABAb impulse response | 50pA
% data_file = 'C:/Users/artemios/Documents/Multiscale_Models_Data/lfp_91.mat'; % GABAb injection at t = 1/2

% data_file = 'C:/Users/artemios/Documents/Multiscale_Models_Data/2023/lfp_py__4_u0.mat';
data_file = 'C:/Users/artemios/Documents/Multiscale_Models_Data/2023/lfp_e1.80_i1.00.mat';

    if nargin > 3, PLOT = varargin{1}; else, PLOT = true; end
    if nargin > 4, data_file = varargin{2}; end

    signal = 'lfp'; % Options: 'vpi', 'vip', 'lfp'
    
    [x_nmm, x_lif, t_nmm, t_lif, v_pi, v_ip, u_lif, lfp_dt] = get_data(signal, x, yy, t, data_file);
    
    [x_nmm, x_lif, t_nmm] = normalization(x_nmm, x_lif, t_nmm, t_lif);
    
    if PLOT
        figure(1);clf;
        plot(t_nmm, x_nmm, 'k');
        hold;
        plot(t_lif, x_lif, '--k');
        legend({'NMM', 'LIF'});
        xlabel('Time (s)');
        title(['Comparison of ' signal]);
        xlim([0.25 max(t_lif)]);
        ylabel('Normalized V_m (a.u.)')
    end
    
    
    harmonic_nmm = do_plot('nmm', signal, x_nmm, t_nmm);
    varargout = {harmonic_nmm};
    harmonic_lif = do_plot('lif', signal, x_lif, t_lif);
%     varargout = {harmonic_lif};
%     varargout{end + 1} = u_lif;
    
    if PLOT
        figure(2); clf;
        ax1 = subplot(2,1,1);
        plot(t,x(:,1)); hold on;
        plot(t,x(:,3));
        title('NMM');
        legend({'x1', 'x3'});
        xlabel('Time (s)');
        ylabel('V_m (mV)');
        box off
        
        plot_lif_results(v_pi, v_ip, lfp_dt); % Won't run if lif results haven't been loaded
        
        % Set the Y limits to the highest values among both subplots:
        figure(2)
        ax2 = subplot(2,1,2);
        limits = [min(ax1.YLim(1), ax2.YLim(1)) max(ax1.YLim(2), ax2.YLim(2))];
        ax1.YLim = limits;
        ax2.YLim = limits;
    end
   
    [w_nmm, w_lif] = do_correlation(x_nmm, x_lif, t_nmm, t_lif);
    do_variance(x_nmm, x_lif, t_nmm, t_lif);
% 	[R,P] = do_spearmanCorr(x_nmm, x_lif)
   
%     varargout = {x_nmm, x_lif, t_nmm, t_lif};
    varargout = {w_nmm, w_lif};
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
    x_nmm = x_nmm(500:end);
    t_nmm = t_nmm(500:end);

    %% LIF
    load(data_file);
    
    trim = 5000; % Samples to remove from the beginning of the LFP_V vector
    LFP_ = LFP_V(trim:end); % LFP(trim:end);
    
    if exist('lfp_dt','var'), dt = lfp_dt; else, dt = 1e-4; end
    
    T = (trim + length(LFP_)) * dt;
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
            x_lif = LFP_; % LIF
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
%     x1_ = [0 mean(diff(v_pi,[],2))];
    x1_ = [0 diff(v_pi,[],2)];
    x3 = mean(v_ip,1);
    x3_ = [0 mean(diff(v_ip,[],2))];
    
    figure; 
    subplot(1,2,1)
    plot(x1*1e3, x3*1e3);
    xlabel('x1');
    ylabel('x3');

    subplot(1,2,2)
    try
        plot(x1*1e3,x1_*1e6);
    catch 
        disp('Couldn''t plot x1 vs x1''');
    end
    xlabel('v');
    ylabel('z');
    
    figure(2);
    subplot(2,1,2)
    plot(t, x1*1e3); hold on;
    plot(t, x3*1e3);
    legend({'x1', 'x3'});
    xlabel('Time (s)');
    ylabel('V_m (mV)');
    title('LIF');
    box off;
    
%     figure; plot3(t * 1e3,x3*1e3,x3_*1e6, 'r')
%     xlabel('t');
%     ylabel('v');
%     zlabel('z');
end

function [w_nmm, w_lif] = do_correlation(x_nmm, x_lif, t_nmm, t_lif)    
    global PLOT
%     [correlation, lags] = xcorr(x_nmm, x_lif);
    [auto_nmm, lag_nmm] = autocorr(x_nmm, 'NumLags', 1000);
    [auto_lif, lag_lif] = autocorr(x_lif, 'NumLags', 1000);
%         
%     figure;
%     [correlation, lags] = xcorr(x_nmm, x_nmm);
%     plot(lags, correlation);
%     hold
%     [correlation, lags] = xcorr(x_lif, x_lif);
%     plot(lags, correlation);
    
    % Compute width of the half section
    w_nmm = 2 * (find(auto_nmm <= 0.5, 1) - 1);
    w_lif = 2 * (find(auto_lif <= 0.5, 1) - 1);
    
	if PLOT
        figure
        plty = [fliplr(auto_nmm) auto_nmm];
        pltx = [fliplr(-lag_nmm) lag_nmm];
        l_nmm = plot(pltx', plty');
        hold
        plty = [fliplr(auto_lif) auto_lif];
        pltx = [fliplr(-lag_lif) lag_lif];
        l_lif = plot(pltx', plty');

        plot([-w_lif/2 w_lif/2], [0.5 0.5], 'Color', l_lif.Color);
        plot([-w_nmm/2 w_nmm/2], [0.5 0.5], '--', 'Color', l_nmm.Color);

        legend({'NMM', 'LIF'});
        xlabel('Lag');
        ylabel('Autocorrelation');
        grid on
        box off
	end
end

function [R,P] = do_spearmanCorr(x_nmm, x_lif)
    idxs = 1:min(length(x_nmm),length(x_lif)); % Adjust the index range so vectors are same size.
%     [R,P] = corrcoef(x_nmm(idxs), x_lif(idxs));
    [R,P] = corr(x_nmm(idxs)',x_lif(idxs)','Type','Spearman');
end

function do_variance(x_nmm, x_lif, t_nmm, t_lif)    
    global PLOT

    window_size = 0.25; % Size of window in seconds
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
    
    if PLOT
        figure; 
        plot(variance_nmm);
        hold on
        plot(variance_lif);    
    end
end
