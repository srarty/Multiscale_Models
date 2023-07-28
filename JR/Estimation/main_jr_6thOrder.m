% this code repricates the structure of the Voss / Schiff FN simulation /
% estimation

clear
% close all
%clc
% dbstop if error

tic

global dT dq dt nn dx log_sig analytic_KF% sampling time step as a global variable

NStates = 6;                               % in the non-augmented JR model
dq = 5; dx = dq + NStates; dy = 1; % dimension of parameter, augmented state and observation
time_vary = 1;
log_sig = 2; % 0:erf | 1:logarithmic_sigmoid | 2:naka_rushton | 3:gompertz
estimate_real_data = 1; % Fit the simulation (0) or an actual EEG (1)
analytic_KF = 0; % If 0 -> Unscented

fct = 'JRfct_6thOrder';                              % Voss Fitzhugh Nagumo function for filtering (propagating sigma points)
obsfct = 'JRobsfct_6thOrder';                    % Voss Fitzhugh Nagumo observation function

N = 100000;                                % number of samples
dT = 0.001;                               % sampling time step (global)
dt = 1*dT;                                 % integration time step
nn = fix(dT/dt);                         % (used in for loop for forward modelling) the integration time step can be small that the sampling (fix round towards zero)

t = 0:dt:(N-1)*dt;

nonlinearity = {'Error function', 'Log sigmoid', 'Naka-Rusthon', 'Gompertz'}; % Labels, not important for the KF or NMM
if estimate_real_data, str = 'Fitting iEEG'; else, str = 'Estimating fwd model'; end

% preallocate for speed
x0 = zeros(NStates,N);              % true state trajectory
%x0b = zeros(NStates,N);              % true state trajectory
xhat = zeros(dx,N);                   % estimate states (and parameters)
xhat2 = zeros(dx,N);                   % estimate states (and parameters)

Pxx = zeros(dx,dx,N);               % estimation covariance
errors = zeros(dx,N);
Ks = zeros(dx,dy,N);                 % Kalman gains

% Initial conditions
x0(:,1) = zeros(NStates,1);                   % initial state

% intial true parameters value
SetParameters_6thOrder;

parameters = {dt,...
    mu,...
    e_0,...
    v_0,...
    r,...
    A,...        % excitatory gain
    a,...        % excitatory time constant
    B,...        % inhibitory gain
    b,...        % inhibitory time constant
    C1,...       % connectivity parameter - excitatory feedback PSP
    C2,...       % connectivity parameter - excitatory feedback firing
    C3,...       % connectivity parameter - inhibitory feedback PSP
    C4,...       % connectivity parameter - excitatory feedback firing
    nrM,...
    nra,...
    nrb,...
    nrs,...
    ga,...
    gb,...
    gc,...
    gd};       

% define input
rng(1)
e = sqrt(dt)*A*a*sigma*randn(N,1);

if time_vary == 1
    % define time varying B and b
    indN = 1:N;
    slopeparam = 10*1/N;
    
    B_series = B+11 ./ (1 + exp(slopeparam*(N/2 - indN)));
    A_series = A*ones(size(B_series)); % A+3 ./ (1 + exp(slopeparam*(N/2 - indN)));
    b_series = b-(b/2) ./ (1 + exp(slopeparam*(N/2 - indN)));
    a_series = a*ones(size(B_series)); %a-(a/2) ./ (1 + exp(slopeparam*(N/2 - indN)));    
    mu_series = mu+110 ./ (1 + exp(0.0001*(N/2 - indN)));
end

% Euler-Maruyama integration
for n=1:N-1
    if time_vary == 1
        parameters = {dt,...
        mu_series(n),...
        e_0,...
        v_0,...
        r,...
        A_series(n),...         % excitatory gain
        a_series(n),...          % excitatory time constant
        B_series(n),...        % inhibitory gain
        b_series(n),...        % inhibitory time constant
        C1,...       % connectivity parameter - excitatory feedback PSP
        C2,...       % connectivity parameter - excitatory feedback firing
        C3,...       % connectivity parameter - inhibitory feedback PSP
        C4,...       % connectivity parameter - excitatory feedback firing
        nrM,...
        nra,...
        nrb,...
        nrs,...            
        ga,...
        gb,...
        gc,...
        gd};   
    end
    x0(:,n+1) = JRint_6thOrder(x0(:,n), parameters)  + [0; 0; 0; e(n); 0; 0];
    %x0(:,n+1) = JRint_6thOrder_Matrix(x0(:,n),parameters)   + [0; 0; 0; e(n); 0; 0];
end

% figure
% for i = 1:6
% subplot(6,3,3*(i-1)+1)
% plot(x0(i,:))
% subplot(6,3,3*(i-1)+2)
% plot(x0b(i,:))
% subplot(6,3,3*(i-1)+3)
% plot(x0(i,:)-x0b(i,:))
% end
% diffx0m = mean(abs(x0(:)-x0b(:)))
% diffx0s = sum(abs(x0(:)-x0b(:)))
% dsfsdfsdfsdf

C = [0 0 1 0 -1 0];           % observation function
v_p = C*x0;

Fs = 1/dt;
NFFT = 2^nextpow2(numel(t)); % Next power of 2 from length of y
Y = fft(v_p-mean(v_p),NFFT)/numel(t);
f = Fs/2*linspace(0,1,NFFT/2+1);

% % Plot single-sided amplitude spectrum.
% f_i_100 = find(f < 100);
% figure
% plot(f(1:f_i_100(end)),2*abs(Y(1:f_i_100(end)))) 
% title('Single-Sided Amplitude Spectrum of y(t)')
% xlabel('Frequency (Hz)')
% ylabel('|Y(f)|')

%%

% ~~~~~~~~~~~~~~~~~~~
% this is estimation from here on
% ~~~~~~~~~~~~~~~~~~~

x = [zeros(dq,N) ; x0];           % form augmented state vector

% observation stuff
R = 0.2^2 * var(JRobsfct_6thOrder(x))*eye(dy,dy);        % observation noise covariance matrix
%R = 0.2^2 * eye(dy,dy)* 7.8242;
y = feval(obsfct,x) + sqrtm(R)*randn(dy,N);             % create noisey data
% toc;

% Simulated data:
% ~~~plot data and noisy data
FS = 16;
fig_width = 40;
fig_height = 10;
figure('units','centimeters',...
    'color','white',...
    'papersize',[fig_width fig_height],...
    'PaperPositionMode','auto',...
    'renderer','painters')

subplot(211),plot(t,y,'k'),axis tight,box off,set(gca,'fontsize',FS),ylabel('iEEG','fontsize',FS)%,ylim([-5 20])
subplot(212),plot(t,y,'r.',t,v_p,'k'),xlim([49.5,50.5]),xlabel('Time (s)','fontsize',FS),box off,set(gca,'fontsize',FS),ylabel('iEEG','fontsize',FS)%,ylim([-5 20])


%% LOAD DATA 
if estimate_real_data
    data = load('Seizure_1.mat');
    N = length(data.Seizure);
    dt = length(data.Seizure)/data.T;
    t = 0:dt:(N-1)*dt;
    y = data.Seizure(:,13)';
else
    if time_vary
        x(1,:) = B_series; 
        x(2,:) = b_series;
        x(3,:) = mu_series; %mu*factor;
        x(4,:) = A_series;
        x(5,:) = a_series;        
    else
        x(1,:) = B * ones(1,size(x,2)); 
        x(2,:) = b * ones(1,size(x,2));
        x(3,:) = mu * ones(1,size(x,2)); %mu*factor;
        x(4,:) = A * ones(1,size(x,2));
        x(5,:) = a * ones(1,size(x,2));
    end
end

% ~~~plot data and noisy data
FS = 16;
fig_estimation = figure('units','centimeters',...
    'color','white',...
    'papersize',[fig_width fig_height],...
    'PaperPositionMode','auto',...
    'renderer','painters');

subplot(211),plot(t,y,'k'),axis tight,box off,set(gca,'fontsize',FS),ylabel('iEEG','fontsize',FS)%,ylim([-5 20])
handle = title(['Forward model (' nonlinearity{log_sig + 1} ')']);
handle.FontSize = 12;
subplot(212),plot(t,y,'k'),xlim([49.5,50.5]),xlabel('Time (s)','fontsize',FS),box off,set(gca,'fontsize',FS),ylabel('iEEG','fontsize',FS)%,ylim([-5 20])

% Initialize states
xhat(dq+1:end, 1) = mean( x0(: , size(x0,2)/2:end) ,2);

% set initial conditions on parameters
if time_vary == 1
    factor = 1;%1.84; % 0.24;
else
    factor = 1.2;%1.6;
end

if dq == 5
    xhat(1,1) = B*factor; 
    xhat(2,1) = b*factor;     % true value = 3.25
    xhat(3,1) = mu*factor; %mu*factor;
    xhat(4,1) = A*factor;
    xhat(5,1) = a*factor;

    % Covariances
    scale = 1e-3;
    Q1 = (scale*B)*B;
    Q2 = (scale*b)*b;
    Q3 = (scale*mu)*mu;
    Q4 = (scale*A)*A;
    Q5 = (scale*a)*a;

    % initialise state error covariance matrix
    Pxx(:,:,1) = diag([Q1,Q2,Q3,Q4,Q5,R*ones(1,6)]);
    Sigma_scale = 1e-7;                                                         % this is our trust parameter for our model
    Sigma = diag([Q1,Q2,Q3,Q4,Q5,Sigma_scale*ones(1,6)]);
elseif dq == 3
    xhat(1,1) = B*factor; 
    xhat(2,1) = b*factor;     % true value = 3.25
    xhat(3,1) = mu*factor; %mu*factor;
    
    % Covariances
    scale = 1e-3;
    Q1 = (scale*B)*B;
    Q2 = (scale*b)*b;
    Q3 = (scale*mu)*mu;%Q3 = sigma*sigma%(scale*mu)*mu;
        
    % initialise state error covariance matrix
    Pxx(:,:,1) = diag([Q1,Q2,Q3,R*ones(1,6)]);
    Sigma_scale = 1e-7;                                                         % this is our trust parameter for our model
    Sigma = diag([Q1,Q2,Q3,Sigma_scale*ones(1,6)]);    
end

Sigma(dx-2,dx-2) = dt*(A*a*sigma)^2;

% Main loop for the recursive estimation
% tic
for k=2:N%2000%N
    
    if analytic_KF == 0 
        [xhat(:,k), Pxx(:,:,k), Ks(:,:,k), Xa] = UKFJR(xhat(:,k-1), Pxx(:,:,k-1), y(:,k),fct,obsfct,parameters,R,Sigma);
    elseif analytic_KF == 1
        [xhat(:,k), Pxx(:,:,k), Ks(:,:,k)] = sigmoid_KF( xhat(:,k-1),Pxx(:,:,k-1),y(:,k),parameters,R,Sigma );
    end

    if dq == 5
        Pxx(1,1,k) = Q1;
        Pxx(2,2,k) = Q2;
        Pxx(3,3,k) = Q3;
        Pxx(4,4,k) = Q4;
        Pxx(5,5,k) = Q5;
    elseif dq == 3
        Pxx(1,1,k) = Q1;
        Pxx(2,2,k) = Q2;
        Pxx(3,3,k) = Q3;
    elseif dq == 2
        Pxx(1,1,k) = Q1;
        Pxx(2,2,k) = Q2;
    elseif dq == 1
        Pxx(1,1,k) = Q1;    
    end

    errors(:,k) = sqrt(diag(Pxx(:,:,k)));
    
end

%%
% Plot parameter estimates
figure
labels = {'Inhibitory gain' 'Inhibitory time constant' 'External input' 'Excitatory gain' 'Excitatory time constant'};
for i = 1:dq
    subplot(dq,1,i)
    plot(t, xhat(i,:),'--r');
    if ~estimate_real_data
        hold;
        plot(t, x(i,:), '-k');
    end
    ylabel(labels{i});
end
subplot(dq,1,1)
title(['Estimation (' nonlinearity{log_sig + 1} ')']);

%% RMSE
y_ = feval(obsfct,xhat);

rmse = sqrt(mean((y(N/2:end) - y_(N/2:end)).^2,2));
disp(['RMSE ' str '(' nonlinearity{log_sig + 1} ') = ' num2str(rmse)]);
% figure
% bar(rmse_x)

%% Plot observation estimate
figure(fig_estimation)
subplot(211), hold; plot(t,y_,'--r'),axis tight,box off,set(gca,'fontsize',FS),ylabel('iEEG','fontsize',FS)%,ylim([-5 20])
handle = title([str ' (' nonlinearity{log_sig + 1} '), RMSE = ' num2str(rmse)]);
handle.FontSize = 12;
subplot(212), hold; plot(t,y_,'--r'),xlim([49.5,50.5]),xlabel('Time (s)','fontsize',FS),box off,set(gca,'fontsize',FS),ylabel('iEEG','fontsize',FS)%,ylim([-5 20])

%% Get estimated firing rates
if estimate_real_data
    [~, f_p, f_e, f_i] = feval(fct, parameters, xhat);
    
    figure;
    subplot(311)
    plot(t,f_p);
    ylabel('f_p');
    title(['Firing rates - ' nonlinearity{log_sig + 1}]);
    
    subplot(312)
    plot(t,f_e);
    ylabel('f_e');
    
    subplot(313)
    plot(t,f_i);
    ylabel('f_i');
    xlabel('Time (s)');
    
end