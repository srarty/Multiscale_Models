% Set parameters for the JR model
%
function params = set_parameters(mode, varargin)

mu = 1; % Default external input value

if nargin < 1
    mode = 'default';
elseif nargin > 1
    % If the function is called with 2 input arguments, the input (u) is 
    % defined by the second argument
    mu = varargin{1};
end

%% Set default parameters:
% Maximum firing rates
params.e0 =  27.81; % Pyramidal % 35.58;% Distribution = 27.81    | non_dsitribution = 35.58
params.e0i = 60; % Inhibitory interneurons %55.53;% Distribution = 60    | non_dsitribution = 55.53

% Sigmoid error function params
% Pyramidal
params.v0_p = 4.074; % Firing Threshold (a)
params.r_p = 3.478; % Sigmoid slope (b)
% Inhibitory
params.v0 = 6.011; % Firing Threshold (a)
params.r = 3.443; % Sigmoid slope (b)

% Gaussian nonlinearity params
% Pyramidal
params.gaussian.a = 26.2;
params.gaussian.b = 10;
params.gaussian.c = 6;
params.gaussian.d = 10;
% Inhibitory
params.gaussiani.a = 50;
params.gaussiani.b = 12.37;
params.gaussiani.c = 7.329;
params.gaussiani.d = 12.4;

% Gompertz nonlinearity params:
% Pyramidal:
params.gompertz.b = 1.512; %1.445; % Distribution = 1.512    | non_dsitribution = 1.445
params.gompertz.c = 1.508; %1.476; % Distribution = 1.508    | non_dsitribution = 1.476
params.gompertz.d = 0.2805; %0.1557;% Distribution = 0.2805    | non_dsitribution = 0.1557
% Interneurons:
params.gompertzi.b = 2.124;%2.017;%% Distribution = 2.124    | non_dsitribution = 2.017
params.gompertzi.c = 1.738;%1.78;%% Distribution = 1.738    | non_dsitribution = 1.78
params.gompertzi.d = 0.1964;%0.2067;%% Distribution = 0.1964    | non_dsitribution = 0.2067

% time constants
params.tau_mi = 0.01; %0.009738; % Membrane time constant - Interneurons (decay)
params.tau_si = 0.0012; %0.001472; % Synaptic time constant - Py->In
params.tau_mp = 0.02;%0.01966;  % Decay tau (membrane time constant) - Pyramidal
params.tau_sp = 0.00525;%0.005608; % Rising tau Pyramidal - In -> Py
params.tau_mrp = 0.02; %0.01975;  % Decay tau recursive excitation
params.tau_srp = 0.0024; %0.002665; % Rising tau recursive excitation
params.tau_mri = 0.01; %0.009229; % Decay tau recursive inhibition
params.tau_sri = 0.00525; %0.006031; % Rising tau recursive inhibition

% Gains
params.alpha_i = -0.5269;       %Inhibitory gain into pyramidal (Interneuron -> Py) %-0.26 <- oscillation
params.alpha_e = 1.124;         % Excitatory gain into interneuron (Py -> Interneuron) (bifurcation: alpha_ei > 392.6)
params.alpha_re = 0.4009;       % Recursive excitatory gain (bifuration: alpha_re > 4.7)
params.alpha_ri = -2.5; %-0.9698;% Recursive inhibitory gain, -2.5 increases Py resting membrane potential
params.alpha_u = 0.06152;        % External excitatory gain into pyramidal (U -> Py)

% Connectivity parameters:
params.c_constant = 1; % Connectivity constant
params.P_pyTOin = 0.395; 	% Probability of connection between Py -> Interneuron
params.P_inTOpy = 0.411;    % Probability of connection between In -> Pyramidal
params.P_pyTOpy = 0.16;    % Probability of connection between Py -> Pyramidal
params.P_inTOin = 0.451;   % Probability of connection between In -> Interneuron

% Integrate and Fire Morphological parameters:
params.g_m_P = 25e-9; % Membrane conductance Pyramidal cells
params.g_m_I = 20e-9; % Membrane conductance Inhibitory cells

params.C_P = 0.5e-9; % Membrane conductance Pyramidal cells
params.C_I = 0.2e-9; % Membrane conductance Inhibitory cells

% External input
params.u = mu; % mean input firing rate

params.dt = 0.001;     % sampling time step         
params.scale = 1; % Scale to fix mismatch in state amplitudes. Not to be confused with the scale in analytic_kalman_filter_2
    
%% Mode specific parameters
switch mode
    case 'default'
        % do nothing
    case 'seizure'
        params.alpha_i = -0.3;
    case 'gabab'
        % Gains:
        params.alpha_i = -0.2635;           % Inhibitory -> Pyramidal
        params.alpha_e = 0.2813;            % Pyramidal -> Inhibitory
        params.alpha_ri = -1.25; %-0.2501;  % Inhibitory -> GABAb
        params.alpha_re = 0.2005;%2*0.2005; % Pyramidal -> Pyramidal
        params.alpha_b = -0.1143;           % GABAb -> Pyramidal %  Note, this is (non-intuituvely) positive during the fit, but should be negative in here.
        params.alpha_eb = 0.2813; %0.1364;  % Py -> GABAb
        params.alpha_u = 0.7075; %0.1364;   % Py -> GABAb
        
        % Probabilities
        params.c_constant = 1;
        params.P_inTOin = 0.5 * params.P_inTOin;
        params.P_inTOpy = 0.5 * params.P_inTOpy;
        
        % Time constants:
        params.tau_sb = 20 * params.tau_sp;
        
        % Nonlinearity
        % Pyramidal:
        params.e0 =  33.32; % 30.8; % Maximum firing rates
        params.gompertz.b = 0.9314;% 1.618;
        params.gompertz.c = 1.024;% 1.457;
        params.gompertz.d = 0.1726;% 0.2382;
        
        % Interneurons:
        params.e0i = 51.62;% 49.3; % Maximum firing rates
        params.gompertzi.b = 2.143;% 2.083;
        params.gompertzi.c = 1.563;% 1.623;
        params.gompertzi.d = 2.449;% 0.2538;
        
        % GABAb:
        params.e0b = 54.44;% 64.3; % Maximum firing rates
        params.gompertzb.b = 2.182;% 2.19;
        params.gompertzb.c = 1.607;% 1.719;
        params.gompertzb.d = 0.2373;% 0.3111;
        
    otherwise
        error('%s rythm not implemented, sorry!', mode);
end
