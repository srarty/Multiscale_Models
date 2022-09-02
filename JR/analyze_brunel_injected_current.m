% Code to analyze Brunel data: Find Brunel's nonlinearity
%
%
% Artemio - August 2021

clear
% close all


%% Options ----------------------------------------------------------------

POPULATION = 'Py'; % 'Py' or 'In'
FUNCTION = 'G'; % 'G' (Gompertz) or 'S' (Sigmoid)

% -------------------------------------------------------------------------

g_m_P = 25e-9;
g_m_I = 20e-9;

%% NMM sigmoid
params = set_parameters('recursive');       % Chose params.u from a constant value in set_params
if strcmp(POPULATION, 'Py'), max_firing_rate = params.e0; elseif strcmp(POPULATION, 'In'), max_firing_rate = params.e0i; else, error('Wrong POPULATION'); end
% if strcmp(POPULATION, 'Py'), max_firing_rate = 120; elseif strcmp(POPULATION, 'In'), max_firing_rate = 185; else, error('Wrong POPULATION'); end
    
x = -20:0.1:100;
nonlinearity = nan(size(x));
for i = 1:numel(x)
    nonlinearity(i) = max_firing_rate * non_linear_sigmoid(x(i), params.r, params.v0);
end
figure
plot(x,nonlinearity, 'LineWidth', 2);
box off
grid on
ylabel('Output');
xlabel('Input');
hold;
plot([min(x) max(x)],[max_firing_rate*0.5 max_firing_rate*0.5],'--k');
plot([params.v0 params.v0], [0 max_firing_rate*1],'--k');
xlabel('Membrane potential (mV)');
ylabel('Spike rate');

%% Load the data
% folder = 'C:\Users\artemios\Documents\GitHub2\mycroeeg\simulations\CUBN\injected_current\twice_ref\';
% folder = 'C:\Users\artemios\Documents\GitHub2\mycroeeg\simulations\CUBN\recurrent_connections\no_inhibitory_input\';
% folder = 'C:\Users\artemios\Documents\GitHub2\mycroeeg\simulations\CUBN\recurrent_connections\inhibitory_input\';
% folder = 'C:\Users\artemios\Documents\GitHub2\mycroeeg\simulations\CUBN\recurrent_inhibition\';
% folder = 'C:\Users\artemios\Documents\GitHub2\mycroeeg\simulations\CUBN\recurrent_excitation\';
% folder = 'C:\Users\artemios\Documents\Multiscale_Models_Data\nonlinearity\';
% folder = 'C:\Users\artemios\Documents\Multiscale_Models_Data\nonlinearity background activity (external input for baseline)\';
% folder = 'C:\Users\artemios\Documents\Multiscale_Models_Data\nonlinearity\tau_e_13ms\';
% folder = 'C:\Users\artemios\Documents\Multiscale_Models_Data\nonlinearity\double_exp\';
% folder = 'C:\Users\artemios\Documents\Multiscale_Models_Data\nonlinearity\double_exp_v2\';
% folder = 'C:\Users\artemios\Documents\Multiscale_Models_Data\Spartan\nonlinearity\';
folder = 'C:\Users\artemios\Documents\Multiscale_Models_Data\Spartan\nonlinearity_connected\';

d = dir([folder '*.mat']);
no_files = numel(d);

membrane_potentials = zeros(1, no_files);
firing_rates = zeros(1, no_files);
inhibitory_rates = zeros(1, no_files);
for ii = 1:no_files
    data_file = [folder d(ii).name];    
    load(data_file);
    L = length(LFP);
    
    if strcmp(POPULATION, 'Py')
        try
            membrane_potentials(ii) = 1000*Vm;
        catch
            membrane_potentials(ii) = 1000 * double(input_current)*1e-12 / g_m_P;
        end
        firing_rates(ii) = mean(R_py(0.2*L : 0.8*L));
        
    elseif strcmp(POPULATION, 'In')
        try
            membrane_potentials(ii) = 1000*Vm_interneurons;
        catch
            membrane_potentials(ii) = 1000 * double(input_current)*1e-12 / g_m_I;
        end
        firing_rates(ii) = mean(R_in(0.2*L : 0.8*L)); %/max_firing_rate;
        
    else
        error('Wrong POPULATION');
        
    end
end

% When there's no spike, R_py and R_in are NaN. Fix it:
firing_rates(isnan(firing_rates)) = 0;

% Force the sigmoid in the data (dodgy)
% if strcmp(FUNCTION, 'S')
%     firing_rates(firing_rates > max_firing_rate) = max_firing_rate;
% end

%%
fig = figure;
scatter(membrane_potentials, firing_rates, 5, 'filled');
xlabel('Membrane potential (mV)');
ylabel('Firing rate');

%% Fit
if strcmp(FUNCTION, 'S')
    % ft = fittype( '2.^-(a.^-(x-b))', 'independent', 'x', 'dependent', 'y' );                          % Double exponential sigmoid
    % ft = fittype( 'c/(1+exp(-a*(x-b)))+d', 'independent', 'x', 'dependent','y' );                     % sigmoid
    ft = fittype( 'alpha*(0.5*erf((x - v0) / (sqrt(2) * r)) + 0.5)', 'independent', 'x', 'dependent', 'y' ); % Error function | a = v0 | b = r
    if strcmp(POPULATION, 'In')
        % Error (In)
        opts = fitoptions(ft);
        opts.StartPoint =  [max_firing_rate 10 5]; % [alpha r v0]
        opts.Lower =   [max_firing_rate 0.1 0];
        opts.Upper =  [max_firing_rate 50 40];    
    else
        % Error (Py)
        opts = fitoptions(ft);
        opts.StartPoint =  [max_firing_rate 10 5];
        opts.Lower =   [1 1 0];
        opts.Upper =  [30 30 100];    
    end
        
elseif strcmp(FUNCTION, 'G')
    ft = fittype( 'a*exp(-b*exp(-d*(x-c)))', 'independent', 'x', 'dependent', 'y' );                  % Gompertz sigmoid
    
    if strcmp(POPULATION, 'In')
        % Gompertz (In)
        opts = fitoptions(ft);
        opts.StartPoint =  [1 1 1 0];
        opts.Lower =  [max_firing_rate 4.0 0.2 0.01];%[1 3.5 0.8 0.1]
        opts.Upper =  [max_firing_rate 100 100 0.15];%[1 100 100 0.1]
    else
        % Gompertz (Py)
        opts = fitoptions(ft);        
        opts.StartPoint =  [1 5 0.001 0.2];
        opts.Lower =  [max_firing_rate  10    5    0.2];%[1 3.5 0.5 0.01];
        opts.Upper =  [max_firing_rate  10    50    10];%[1 100 100 0.125];
    end    
else
    error('Wrong POPULATION');
end

fitresult_Py = fit(membrane_potentials', firing_rates', ft, opts) % With options
% fitresult_Py = fit(membrane_potentials', firing_rates', ft) % No options
    
%% Define nonlinear function according to population
if strcmp(FUNCTION, 'S')
    f_nonlinearity = @(x) fitresult_Py.alpha * non_linear_sigmoid(x, fitresult_Py.r, fitresult_Py.v0);
elseif strcmp(FUNCTION, 'G')
    params.gompertz.a = fitresult_Py.a;
    params.gompertz.b = fitresult_Py.b;
    params.gompertz.c = fitresult_Py.c;
    params.gompertz.d = fitresult_Py.d;
    
    f_nonlinearity = @(x) gompertz(x, params);
end

%%
figure(fig);
yyaxis right
hold
nonlinearity_post = nan(size(x));
for i = 1:numel(x)
    nonlinearity_post(i) = f_nonlinearity(x(i));
end
% plot(x, max_firing_rate * nonlinearity_post);xlabel('Membrane potential (mV)');
plot(x, nonlinearity_post);xlabel('Membrane potential (mV)');
ylabel('Firing rate (spikes/s)');
ylim([0 max_firing_rate]);
yyaxis left
ylim([0 max_firing_rate]);
