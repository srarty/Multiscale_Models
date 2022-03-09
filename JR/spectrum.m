
signal = 'lfp'; % vpi, vip, lfp
model = 'lif'; % lif or nmm

T = 1;

if strcmp('nmm', model)
    dt = 0.001;
    t = linspace(0,T,1/dt);
    
    switch signal
        case 'vpi'
            x_ = x(1,:); % State 1 % NMM
        case 'vip'
            x_ = x(3,:); % State 3 % NMM
        case 'lfp'
            x_ = y; % NMM
        otherwise
            error('Wrong options (signal)');
    end
    
    w = 128; % Window size
    so = 120; % Samples overlap
    freqbins = 128; %Evaluate the spectrum at (128/2)+1=65 frequencies and (length(x)−120)/(128−120)=235 time bins.    
    xrange = [0 500]; % xlim values
    
    disp(['freq bins', num2str( (freqbins/2)+1 )]);
    disp(['time bins', num2str( (length(x_)-so)/(freqbins-so) )]);
    
elseif strcmp('lif', model)
    data_file = 'C:/Users/artemios/Documents/Multiscale_Models_Data/lfp_4.mat';
%     data_file = 'C:/Users/artemios/Documents/Multiscale_Models_Data/spartan/lfp_22.mat';
    load(data_file);
    dt = 0.0001;
    t = linspace(0,T,1/dt);
        
    switch signal
        case 'vpi'
            x_ = mean(v_pi,1); % State 3? % LIF
        case 'vip'
            x_ = mean(v_ip,1); % State 1? % LIF
        case 'lfp'
            x_ = LFP_V; % LIF
        otherwise
            error('Wrong options (signal)');
    end
    
    
    w = 1280; % Window size
    so = 1200; % Samples overlap
    freqbins = 1280;
    xrange = [0 0.5]; % xlim values
    
    disp(['freq bins', num2str( (freqbins/2)+1 )]);
    disp(['time bins', num2str( (length(x_)-so)/(freqbins-so) )]);
    
else
    error('Wrong options (model)');
end

x_ = x_ - mean(x_);

%% FFT

Fs = 1/dt;
L = length(t);
n = 2^nextpow2(L);
X = fft(x_,n);
P2 = abs(X/L);
P1 = P2(:,1:n/2+1);
P1(:,2:end-1)=2*P1(:,2:end-1);

figure; 
plot(0:(Fs/n):(Fs/2-Fs/n), P1(1:n/2));
xlabel('Frequency (Hz)');
ylabel('Power');

xlim([0 500]);

%% Spectrogram
figure;
spectrogram(x_, w, so, freqbins, T/dt);
xlim(xrange);
