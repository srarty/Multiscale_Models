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

function spectrum(x, y, t)
    close all
%     data_file = 'C:/Users/artemios/Documents/Multiscale_Models_Data/lfp_35.mat';
    data_file = 'C:/Users/artemios/Documents/Multiscale_Models_Data/spartan/lfp_55.mat';
    signal = 'lfp'; % vpi, vip, lfp
    
    [x_nmm, x_lif, t_nmm, t_lif, v_pi, v_ip] = get_data(signal, x, y, t, data_file);
    
    [x_nmm, x_lif, t_nmm] = normalization(x_nmm, x_lif, t_nmm, t_lif);

    do_plot('nmm', signal, x_nmm, t_nmm); % After normalization, t_lif applies for both signals, because it includes an interpolation
    do_plot('lif', signal, x_lif, t_lif);
    plot_lif_results(v_pi, v_ip, t_lif); % Won't run if lif results haven't been loaded
    
    do_correlation(x_nmm, x_lif, t_nmm, t_lif);
end

function [x_nmm, x_lif, t_nmm, t_lif, v_pi, v_ip] = get_data(signal, x, y, t, data_file)
    %% NMM
    t_nmm = t;
    switch signal
        case 'vpi'
            x_nmm = x(:,1); % State 1 % NMM; % change dimm if model forward came from diff eqs instead of the nmm_toolbox
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

    %% LIF
    load(data_file);
    dt = 1e-4;
    T = length(LFP_V) * dt;
    t_lif = linspace(0,T,T/dt);
    LFP_V = LFP_V;
    
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
function varargout = do_plot(model, signal, x_, t)

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

    disp(['freq bins = ', num2str( (freqbins/2)+1 )]);
    disp(['time bins = ', num2str( (length(x_)-so)/(freqbins-so) )]);
            
    %% FFT

    Fs = 1/dt;
    L = length(t);
    n = 2^nextpow2(L);
    X = fft(x_,n);
    P2 = abs(X/L);
    P1 = P2(:,1:n/2+1);
    P1(:,2:end-1)=2*P1(:,2:end-1);

    f = figure(100); 
    f.Position([3 4]) = [1230 420];
    
    if strcmp('lif', model)
        subplot(1,2,2); 
        titlestr = 'LIF';
    else
        subplot(1,2,1); 
        titlestr = 'NMM';
    end
    
    plot(0:(Fs/n):(Fs/2-Fs/n), P1(1:n/2));
    title(titlestr);
    xlabel('Frequency (Hz)');
    ylabel(['Power [', ystr, ']']);

    xlim([0 200]);

    %% Spectrogram
    f = figure(101);
    f.Position([3 4]) = [1230 420];
    if strcmp('lif', model), subplot(1,2,2); else, subplot(1,2,1); end
    spectrogram(x_, w, so, freqbins, Fs, 'yaxis','power');
    [~,~,~,ps] = spectrogram(x_, w, so, freqbins, Fs, 'yaxis','power');

    title([titlestr, ' (', ystr, ')']);
    ax = gca;
%     ax.YScale = 'log';
    ylim(freqrange);
    
end

%%
function plot_lif_results(v_pi, v_ip, t)
    x1 = mean(v_pi,1);
    x1_ = [0 mean(diff(v_pi,[],2))];
    x3 = mean(v_ip,1);
    x3_ = [0 mean(diff(v_ip,[],2))];
    
    figure; plot(x1*1e3, x3*1e3);
    xlabel('x1');
    ylabel('x3');

%     figure; plot(x1*1e3,x1_*1e6)
%     xlabel('v');
%     ylabel('z');
%     
%     figure; plot3(t * 1e3,x3*1e3,x3_*1e6, 'r')
%     xlabel('t');
%     ylabel('v');
%     zlabel('z');
end

function do_correlation(x_nmm, x_lif, t_nmm, t_lif)    
    [correlation, lags] = xcorr(x_nmm, x_lif);
    auto_nmm = autocorr(x_nmm);
    auto_lif = autocorr(x_lif);
    
    figure;
    plot(lags, correlation);
    
    figure;
    bar(auto_nmm); hold
    bar(auto_lif);
end
