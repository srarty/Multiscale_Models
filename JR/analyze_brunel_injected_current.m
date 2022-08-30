% Code to analyze Brunel data: Find Brunel's nonlinearity
%
%
% Artemio - August 2021

clear
close all


%% Options ----------------------------------------------------------------

POPULATION = 'Py'; % 'Py' or 'In'

% -------------------------------------------------------------------------

%% NMM sigmoid
params = set_parameters('recursive');       % Chose params.u from a constant value in set_params
% if strcmp(POPULATION, 'Py'), max_firing_rate = params.e0; elseif strcmp(POPULATION, 'In'), max_firing_rate = params.e0i; else, error('Wrong POPULATION'); end
if strcmp(POPULATION, 'Py'), max_firing_rate = 120; elseif strcmp(POPULATION, 'In'), max_firing_rate = 250; else, error('Wrong POPULATION'); end
    
x = -20:0.1:50;
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
folder = 'C:\Users\artemios\Documents\Multiscale_Models_Data\nonlinearity\double_exp_v2\';

d = dir([folder '*.mat']);
no_files = numel(d);

membrane_potentials = zeros(1, no_files);
pyramidal_rates = zeros(1, no_files);
inhibitory_rates = zeros(1, no_files);
for ii = 1:no_files
    data_file = [folder d(ii).name];    
    load(data_file);
    
    if strcmp(POPULATION, 'Py')
        membrane_potentials(ii) = 1000*Vm;
        pyramidal_rates(ii) = mean(R_py);
        
    elseif strcmp(POPULATION, 'In')
        membrane_potentials(ii) = 1000*Vm_interneurons;
        pyramidal_rates(ii) = mean(R_in); %/max_firing_rate;
        
    else
        error('Wrong POPULATION');
        
    end
end
% % Uncomment the following two lines to fit the Inhibitory population to a
% % gompertz function instead of an sigmoid error function
POPULATION = 'Py';
warning('Changin to Py');

% When there's no spike, R_py is NaN. Fix it:
pyramidal_rates(isnan(pyramidal_rates)) = 0;
% pyramidal_rates(pyramidal_rates > max_firing_rate) = max_firing_rate;

%%
fig = figure;
% yyaxis left
scatter(membrane_potentials, pyramidal_rates, 5, 'filled');
xlabel('Membrane potential (mV)');
ylabel('Pyramidal firing rate');

%% Fit
if strcmp(POPULATION, 'In')
    % ft = fittype( '2.^-(a.^-(x-b))', 'independent', 'x', 'dependent', 'y' );                          % Double exponential sigmoid
    % ft = fittype( 'c/(1+exp(-a*(x-b)))+d', 'independent', 'x', 'dependent','y' );                     % sigmoid
    ft = fittype( '(0.5*erf((x - a) / (sqrt(2) * b)) + 0.5)', 'independent', 'x', 'dependent', 'y' ); % Error function | a = v0 | b = r
    
    % Error (In)
    opts = fitoptions(ft);
    opts.StartPoint =  [10 5];
    opts.Lower =   [16 10];
    opts.Upper =  [16 10];    
        
elseif strcmp(POPULATION, 'Py')
    ft = fittype( 'a*exp(-b*exp(-d*(x-c)))', 'independent', 'x', 'dependent', 'y' );                  % Gompertz sigmoid
    
%     % Gompertz (Py)
    opts = fitoptions(ft);
    opts.StartPoint =  [1 1 1 0];
    opts.Lower =  [1 3.5 0.5 0.01];%[1 3 1 0.15];
    opts.Upper =  [1 100 100 0.125];%[1 100 100 0.15];

    % Gompertz (In)
%     opts = fitoptions(ft);
%     opts.StartPoint =  [1 1 1 0];
%     opts.Lower =  [1 3.5 0.8 0.1];%[1 3 1 0.15];
%     opts.Upper =  [1 100 100 0.1];%[1 100 100 0.15];
    
else
    error('Wrong POPULATION');
end

fitresult_Py = fit(membrane_potentials', pyramidal_rates', ft, opts) % With options
% fitresult_Py = fit(membrane_potentials', pyramidal_rates', ft) % No options
    
% Define nonlinear function according to population
if strcmp(POPULATION, 'In')
    f_nonlinearity = @(x) non_linear_sigmoid(x, fitresult_Py.b, fitresult_Py.a);
elseif strcmp(POPULATION, 'Py')
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
plot(x, max_firing_rate * nonlinearity_post);xlabel('Membrane potential (mV)');
ylabel('Firing rate (spikes/s)');
yyaxis left
ylim([0 max_firing_rate]);
