% Code to analyze Brunel data: Find Brunel's nonlinearity
%
%
% Artemio - August 2021

% clear
% close all


%% Options ----------------------------------------------------------------

POPULATION = 'B'; % 'Py' or 'In'
FUNCTION = 'G'; % 'G' (Gompertz) or 'S' (Sigmoid) or 'Ga' (Gaussian) or 'B' (Bas-Jan Zandt 2014)

% -------------------------------------------------------------------------

% LIF parameters
g_m_P = 25e-9;
g_m_I = 20e-9;

C_P = 0.5e-9;
C_I = 0.2e-9;

%% NMM sigmoid
params = set_parameters('default');       % Chose params.u from a constant value in set_params
if strcmp(POPULATION, 'Py'), max_firing_rate = params.e0; elseif strcmp(POPULATION, 'In')||strcmp(POPULATION, 'B'), max_firing_rate = params.e0i; else, error('Wrong POPULATION'); end
% if strcmp(POPULATION, 'Py'), max_firing_rate = 10; elseif strcmp(POPULATION, 'In'), max_firing_rate = 10; else, error('Wrong POPULATION'); end
    
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
% folder = 'C:\Users\artemios\Documents\Multiscale_Models_Data\nonlinearity\';
% folder = 'C:\Users\artemios\Documents\Multiscale_Models_Data\nonlinearity background activity (external input for baseline)\';
% folder = 'C:\Users\artemios\Documents\Multiscale_Models_Data\nonlinearity\tau_e_13ms\';
% folder = 'C:\Users\artemios\Documents\Multiscale_Models_Data\nonlinearity\double_exp\';
% folder = 'C:\Users\artemios\Documents\Multiscale_Models_Data\nonlinearity\double_exp_v2\';
% folder = 'C:\Users\artemios\Documents\Multiscale_Models_Data\Spartan\nonlinearity\';
% folder = 'C:\Users\artemios\Documents\Multiscale_Models_Data\Spartan\nonlinearity_disconnected\';
% folder = 'C:\Users\artemios\Documents\Multiscale_Models_Data\Spartan\nonlinearity_disconnected_noTha\';
% folder = 'C:\Users\artemios\Documents\Multiscale_Models_Data\Spartan\nonlinearity_connected\';
% folder = 'C:\Users\artemios\Documents\Multiscale_Models_Data\Spartan\nonlinearity_connected_noTha\';
% folder = 'C:\Users\artemios\Documents\Multiscale_Models_Data\Spartan\nonlinearity_I_Tha_disconnected\';
% folder = 'C:\Users\artemios\Documents\Multiscale_Models_Data\Spartan\nonlinearity_I_Tha\';
% folder = 'C:\Users\artemios\Documents\Multiscale_Models_Data\Spartan\nonlinearity_distribution\';
folder = 'C:\Users\artemios\Documents\Multiscale_Models_Data\Spartan\nonlinearity_three_pop\';


d = dir([folder '*.mat']);
no_files = numel(d);

membrane_potentials = zeros(1, no_files);
firing_rates = zeros(1, no_files);
potential_integral = zeros(1, no_files);
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
        potential_integral(ii) = -(sum(I_py - I_py_tha)/C_P)/1000;
        firing_rates(ii) = mean(R_py(0.2*L : 0.8*L));
        
    elseif strcmp(POPULATION, 'In') || strcmp(POPULATION, 'B')
        try
            membrane_potentials(ii) = 1000*Vm_interneurons;
        catch
            membrane_potentials(ii) = 1000 * double(input_current)*1e-12 / g_m_I;
        end
        potential_integral(ii) = -(sum(I_in - I_in_tha)/C_I)/2000;
%         firing_rates(ii) = mean(R_in(0.2*L : 0.8*L));
        firing_rates(ii) = mean(R_b(0.2*L : 0.8*L));
        
    else
        error('Wrong POPULATION');
        
    end
end

% When there's no spike, R_py and R_in are NaN. Fix it:
firing_rates(isnan(firing_rates)) = 0;

% Ignore higher values to find a reasonable max_firing_rate (dodgy)
warning('Remove the following three dodgy lines')
% membrane_potentials(firing_rates > 35) = [];
% potential_integral(firing_rates > 35) = [];
% firing_rates(firing_rates > 35) = [];
membrane_potentials(firing_rates > 55) = [];
potential_integral(firing_rates > 55) = [];
firing_rates(firing_rates > 55) = [];

% Sort values
[potential_integral, idx] = sort(potential_integral);
membrane_potentials = membrane_potentials(idx);
firing_rates = firing_rates(idx);

% For Gaussian, concatenate a flipped firing_rate vector
if strcmp(FUNCTION, 'Ga')
    firing_rates = [firing_rates flip(firing_rates(1:end-1))];
    potential_integral = [potential_integral potential_integral(1:end-1)+1+(max(potential_integral)-min(potential_integral))];
    membrane_potentials = [membrane_potentials membrane_potentials(1:end-1)+1+(max(membrane_potentials)-min(membrane_potentials))];
end


%%
fig = figure;
scatter(membrane_potentials, firing_rates, 5, 'filled');hold
scatter(potential_integral, firing_rates, 5, 'x');
xlabel('Membrane potential (mV)');
ylabel('Firing rate');

%% Fit
if strcmp(FUNCTION, 'S') % Sigmoid
%     ft = fittype( '2.^-(a.^-(x-b))', 'independent', 'x', 'dependent', 'y' );                          % Double exponential sigmoid
%     ft = fittype( 'c/(1+exp(-a*(x-b)))+d', 'independent', 'x', 'dependent','y' );                     % sigmoid
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
        opts.Lower =   [25.87 1 0];
        opts.Upper =  [25.87 30 100];
    end
        
elseif strcmp(FUNCTION, 'G') % Gompertz
    ft = fittype( 'a*exp(-b*exp(-d*(x-c)))', 'independent', 'x', 'dependent', 'y' ); % Gompertz sigmoid
    
    if strcmp(POPULATION, 'In')
        % Gompertz (In)
        opts = fitoptions(ft);
        opts.StartPoint =  [10 1 1 0];
        opts.Lower =  [max_firing_rate-15 -100 -100 -10];
        opts.Upper =  [max_firing_rate+15 100 100 10];
    else
        % Gompertz (Py)
        opts = fitoptions(ft);        
        opts.StartPoint =  [10 1 1 0];
        opts.Lower =  [max_firing_rate-15 -100 -100 -10];
        opts.Upper =  [max_firing_rate+15 100 100 10];
    end    
    
elseif strcmp(FUNCTION, 'Ga')% Gaussian
    ft = fittype( 'a*exp(-((x-b)/c).^2)', 'independent', 'x', 'dependent', 'y' ); % Gompertz sigmoid
    
    if strcmp(POPULATION, 'In')
        % Gaussian (In)
        opts = fitoptions(ft);
        opts.StartPoint =  [24 10 7];
        opts.Lower =  [max_firing_rate-10 -100 -10];
        opts.Upper =  [max_firing_rate+10 100 100 10];
    else
        % Gaussian (Py)
        opts = fitoptions(ft);        
        opts.StartPoint =  [24 10 7];
        opts.Lower =  [max_firing_rate-10 10 6];
        opts.Upper =  [max_firing_rate+10 10 6];
    end    

elseif strcmp(FUNCTION, 'B')    % Bas-Jan Zandt et al 2014
    ft = fittype( '( 1/(s*sqrt(2*pi)) ) * exp( -0.5*((x-m)/s)^2 )', 'independent', 'x', 'dependent', 'y' ); % Gompertz sigmoid
    
    if strcmp(POPULATION, 'In')
        % Bas-Jan Zandt (In)
        opts = fitoptions(ft);
        opts.StartPoint =  [0 0];
        opts.Lower =  [0 0];
        opts.Upper =  [0 0];
    else
        % Bas-Jan Zandt (In)
        opts = fitoptions(ft);        
        opts.StartPoint =  [0 0];
        opts.Lower =  [0 0];
        opts.Upper =  [0 0];
    end    
else
    error('Wrong POPULATION');
end

fitresult_Py = fit(membrane_potentials', firing_rates', ft, opts) % With options
% fitresult_Py = fit(potential_integral', firing_rates', ft, opts) % With options
% fitresult_Py = fit(potential_integral', firing_rates', ft) % No options
    
%% Define nonlinear function according to population
if strcmp(FUNCTION, 'S')
    f_nonlinearity = @(x) fitresult_Py.alpha * non_linear_sigmoid(x, fitresult_Py.r, fitresult_Py.v0);
elseif strcmp(FUNCTION, 'G')
    f_nonlinearity = @(x) gompertz(x, fitresult_Py.a, fitresult_Py.b, fitresult_Py.c, fitresult_Py.d);
elseif strcmp(FUNCTION, 'Ga')
    f_nonlinearity = @(x) fitresult_Py.a*exp(-((x-fitresult_Py.b)/fitresult_Py.c).^2);
elseif strcmp(FUNCTION, 'B')
    f_nonlinearity = @(x) ( 1/(fitresult_Py.s*sqrt(2*pi)) ) * exp( -0.5*((x-fitresult_Py.m)/fitresult_Py.s)^2 );
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
ylim([0 max_firing_rate+15]);
yyaxis left
ylim([0 max_firing_rate+15]);
% ylim([0 40]);
% yyaxis left
% ylim([0 40]);
