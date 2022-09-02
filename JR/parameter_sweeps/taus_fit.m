%% taus_fit
%% Goes through the files generated locally by the LIF file 'find_taus.py'.
% Files saved in "...\Multiscale_Models_Data\taus_fit\" are simulated on a
% single neuron model, where there is a single post-synaptic spike, by 
% varying alpha_xy (the gain, called 'j' in the LIF). _xy refers to the
% synapses, where x is the pre-synaptic population and y is the
% post-synaptic population.
%
% Artemio - August 2022

PRE = 'i'; % PRE = {'i', 'e'} % 'i'=inhibitory, 'e'=excitatory
POST = 'p'; % POST = {'i', 'p'} % 'i'=interneuron, 'p'=pyramidal 

root = 'C:\Users\artemios\Documents\Multiscale_Models_Data\taus_fit\';
d = dir([root POST '*_' PRE '*.mat']);

% fix taus
switch [PRE POST]
    case 'ii'
        lower_limits = [0.0001 0.0001 0];%
        upper_limits = [0.5 0.5 10];
    case 'ip'
        lower_limits = [0.0001 0.0001 0]; % b= 0.005608
        upper_limits = [0.5 1 10];
    case 'ep'
        lower_limits = [0.0001 0.002665 0];
        upper_limits = [0.1 0.002665 10];
    case 'ei'
        lower_limits = [0.002 0.001472 0];
        upper_limits = [0.1 0.001472 10];
end

% Declare empty array
j = nan(size(d));
tau_s = nan(size(d));
tau_m = nan(size(d));
alpha_nmm = nan(size(d));
gf = nan(size(d)); % goodness of fit
for i = 1:length(d)
    % Retrieve current name
    n = d(i).name;
    
    % Find the position alpha in the file name
    alpha_idx = findstr(n,'psp_') + 4;
    end_idx = findstr(n, '.mat') - 1;
    
    % Retrieve alpha
    j(i) = -1*str2num(n(alpha_idx : end_idx));    
    
    % Load the file
    data = load([root d(i).name]);

    
    %% Fit a double exponential function to the PSP
    % Load the ipsp or epsp into a single variable psp:
    if isfield(data,'ipsp'), psp = -data.ipsp; negative_gain = true; elseif isfield(data,'epsp'), psp = data.epsp; negative_gain = false; else, error('The file didn''t have an epsp or ipsp.'); end
    % Adjuisting psp to remove the first sample, which is zero in Brian2:
    psp(1) = [];
    T = linspace(0,0.3,length(psp));

    ft = fittype( 'c*(exp(-t/a)-exp(-t/b))', 'independent', 't', 'dependent', 'y');

    opts = fitoptions(ft);
    opts.StartPoint = [0.008339 0.01001 1];
    opts.Lower = lower_limits; %[0.0001 0.0001 0];     
    opts.Upper = upper_limits; %[0.1 0.2 10];    
    opts.Robust = 'Off';
    [fitresult_alpha, G] = fit(T', psp', ft, opts); % With options
    
    tau_s(i) = fitresult_alpha.b;
    tau_m(i) = fitresult_alpha.a;
    alpha_nmm(i) = ((-1)^negative_gain) * fitresult_alpha.c; % The -1^neg_gain makes the gain negative when ipsp and positive when epsp
    gf(i) = G.rsquare;
    
%    % Uncomment to see the fit of each curve:
%     figure(1);
%     cla;
%     plot(fitresult);
%     hold
%     plot(T,psp,'b--');
%     hold
    
end

[j,idx] = sort(j);
tau_s = tau_s(idx);
tau_m = tau_m(idx);
alpha_nmm = alpha_nmm(idx);
gf = gf(idx);

%% Fit alpha to a rect
ft = fittype( 'm*x + b', 'independent', 'x', 'dependent', 'y');
L = length(j);
if negative_gain
    fit_range = 2:L;
else
    fit_range = 1:L-1;
end
fitresult_alpha = fit(j(fit_range), alpha_nmm(fit_range), ft) % With options

%% Plot
figure;
subplot(211)
plot(j, tau_m*1e3,'o');
hold; plot(j, tau_s*1e3,'x');
legend({'\tau_{m}' '\tau_{s}'});
ylabel('Time constants (ms)');
box off

subplot(212)
plot(j, alpha_nmm,'ok');
hold
plot(fitresult_alpha);
legend('off');
% yyaxis right;
% plot(alpha, gf);
xlabel('j (conductance gain in LIF)');
ylabel('\alpha (gain in NMM)');
box off

subplot(211);
title([PRE '\rightarrow' POST]);