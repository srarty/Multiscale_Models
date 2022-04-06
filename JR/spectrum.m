% Computes the FFT and Spectrogram of the states and output of each model
%
% Artemio - March 2022

function spectrum(x, y, t)
    signal = 'lfp'; % vpi, vip, lfp
    do_plot('nmm', signal, x, y, t);
    [v_pi, v_ip, t] = do_plot('lif', signal);
    plot_lif_results(v_pi, v_ip, t); % Won't run if lif results haven't been loaded
end

function varargout = do_plot(model, signal, varargin)
    if nargin > 2
        x = varargin{1};
        y = varargin{2};
        t = varargin{3};
    end

    % signal = 'vip'; % vpi, vip, lfp
    % model = 'nmm'; % lif or nmm

    if strcmp('nmm', model)
        T = t(end);
        dt = T/length(x);
        t = linspace(0,T,length(x));

        titlestr = 'NMM';

        switch signal
            case 'vpi'
                x_ = x(:,1); % State 1 % NMM; % change dimm if model forward came from diff eqs instead of the nmm_toolbox
                x_ = x_';
%                 x_ = x(1,:);
                ystr = 'V_{pi}';
            case 'vip'
                x_ = x(:,3); % State 3 % NMM
                x_ = x_';
%                 x_ = x(3,:);
                ystr = 'V_{ip}';
            case 'lfp'
                x_ = y'; % NMM
%                 x_ = y; % NMM
                ystr = 'V_m (Py)';
            otherwise
                error('Wrong options (signal)');
        end

        w = 128; % Window size
        so = 120; % Samples overlap
        freqbins = 10 * 128; %Evaluate the spectrum at (128/2)+1=65 frequencies and (length(x)−120)/(128−120)=235 time bins.    
        freqrange = [0 0.1]; % xlim values

    elseif strcmp('lif', model)
        data_file = 'C:/Users/artemios/Documents/Multiscale_Models_Data/lfp_22.mat';
%         data_file = 'C:/Users/artemios/Documents/Multiscale_Models_Data/spartan/lfp_43.mat';
        load(data_file);
        dt = 1e-4;
        T = length(LFP_V) * dt;
        t = linspace(0,T,T/dt);
        LFP_V = LFP_V;
        titlestr = 'LIF';

        switch signal
            case 'vpi'
                x_ = mean(v_pi,1); % State 3? % LIF
                ystr = 'V_{pi}';
            case 'vip'
                x_ = mean(v_ip,1); % State 1? % LIF
                ystr = 'V_{ip}';
            case 'lfp'
                x_ = LFP_V; % LIF
                ystr = 'V_m (Py)';
            otherwise
                error('Wrong options (signal)');
        end


        w = 1280; % Window size
        so = 1200; % Samples overlap
        freqbins = 10 * 1280;
        freqrange = [0 0.1]; % xlim values
        
        varargout = {v_pi, v_ip, t};
    else
        error('Wrong options (model)');
    end


    disp(['freq bins = ', num2str( (freqbins/2)+1 )]);
    disp(['time bins = ', num2str( (length(x_)-so)/(freqbins-so) )]);

    x_ = x_ - mean(x_);
    x_ = x_ / rms(x_); % Normalize by something else, RMS or smth

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
    if strcmp('lif', model), subplot(1,2,2); else, subplot(1,2,1); end
    plot(0:(Fs/n):(Fs/2-Fs/n), P1(1:n/2));
    title(titlestr);
    xlabel('Frequency (Hz)');
    ylabel(['Power [', ystr, ']']);

    xlim([0 100]);

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


function plot_lif_results(v_pi, v_ip, t)
    x1 = mean(v_pi,1);
    x1_ = [0 mean(diff(v_pi,[],2))];
    x3 = mean(v_ip,1);
    x3_ = [0 mean(diff(v_ip,[],2))];
    
    figure; plot(x1*1e3, x3*1e3);
    xlabel('x1');
    ylabel('x2');

%     figure; plot(x1*1e3,x1_*1e6)
%     xlabel('v');
%     ylabel('z');
%     
%     figure; plot3(t * 1e3,x3*1e3,x3_*1e6, 'r')
%     xlabel('t');
%     ylabel('v');
%     zlabel('z');
end
