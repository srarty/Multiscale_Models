% Code to analyze Brunel data: Find Brunel's nonlinearity
%
%
% Artemio - August 2021

clear

%% NMM sigmoid
params = set_parameters('allen');       % Chose params.u from a constant value in set_params
x = -20:0.1:50;
nonlinearity = nan(size(x));
for i = 1:numel(x)
    nonlinearity(i) = params.e0 * non_linear_sigmoid(x(i), params.r, params.v0);
end
figure
plot(x,nonlinearity, 'LineWidth', 2);
box off
grid on
ylabel('Output');
xlabel('Input');
hold;
plot([min(x) max(x)],[params.e0*0.5 params.e0*0.5],'--k');
plot([params.v0 params.v0], [0 params.e0*1],'--k');
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
folder = 'C:\Users\artemios\Documents\Multiscale_Models_Data\nonlinearity\double_exp\';

d = dir([folder '*.mat']);
no_files = numel(d);

input_rate = 0.1:0.1:5;

membrane_potentials = zeros(1, no_files);
pyramidal_rates = zeros(1, no_files);
inhibitory_rates = zeros(1, no_files);
for ii = 1:no_files
    data_file = [folder d(ii).name];    
    load(data_file);
    
%     membrane_potentials(ii) = 1000*Vm;
%     pyramidal_rates(ii) = mean(R_py);
    membrane_potentials(ii) = 1000*Vm_interneurons;
    pyramidal_rates(ii) = mean(R_in); %/params.e0;
end

% % Append zeros in the beginning for a better fit
% membrane_potentials = [0 0 membrane_potentials 40 50];
% pyramidal_rates = [0 0 pyramidal_rates 70 70];

% When there's no spike, R_py is NaN. Fix it:
pyramidal_rates(isnan(pyramidal_rates)) = 0;


% warning('remove next line, it''s dodgy');
% pyramidal_rates(1:19) = 0;
%%
fig = figure;
yyaxis left
scatter(membrane_potentials, pyramidal_rates, 5, 'filled');
xlabel('Membrane potential (mV)');
ylabel('Pyramidal firing rate');

%% Fit
ft = fittype( '(0.5*erf((x - a) / (sqrt(2) * b)) + 0.5)', 'independent', 'x', 'dependent', 'y' ); % Error function | a = v0 | b = r
% ft = fittype( '2.^-(a.^-(x-b))', 'independent', 'x', 'dependent', 'y' );                          % Double exponential sigmoid
% ft = fittype( 'a*exp(-b*exp(-d*(x-c)))', 'independent', 'x', 'dependent', 'y' );                  % Gompertz sigmoid
% ft = fittype( 'c/(1+exp(-a*(x-b)))+d', 'independent', 'x', 'dependent','y' );                     % sigmoid

opts = fitoptions(ft);
% Gompertz (Py)
% opts.StartPoint =  [1 1 1 0];
% opts.Lower =   [1 3 1 0.15];
% opts.Upper =  [1 100 100 0.15];

% Error (In)
opts.StartPoint =  [10 5];
opts.Lower =   [15 8];
opts.Upper =  [15 8];

fitresult_Py = fit(membrane_potentials', pyramidal_rates', ft, opts) % With options
% fitresult_Py = fit(membrane_potentials', pyramidal_rates', ft) % No options


figure(fig);
yyaxis right
hold
plot(fitresult_Py);
ylim([0 1])
xlabel('Membrane potential (mV)');
ylabel('Pyramidal firing rate (spikes/s)');
yyaxis left
% ylim([0 1]);
% ylim([0 250]);
% ylim([0 params.e0]);

% %% Plot firing rates of both populations
% figure
% scatter(input_rate(1:length(pyramidal_rates)-1), pyramidal_rates(2:end),5,'filled');
% hold
% scatter(input_rate(1:length(inhibitory_rates)-1), inhibitory_rates(2:end), 5, 'r','filled');
% xlabel('External input (spikes/ms)/cell')
% ylabel('Spike rate (spikes/s)')
% legend({'Pyramidal', 'Inhibitory'})
