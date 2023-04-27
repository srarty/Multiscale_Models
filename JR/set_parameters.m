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
params.tau_mi = 0.01; %0.009738; % Membrane time constant - Interneurons (decay) [Galarreta 2001]
params.tau_si = 0.0012; %0.001472; % Synaptic time constant - Py->In [Galarreta 2001]
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
params.P_pyTOin = 0.43; % 0.395; 	% Probability of connection between Py -> Interneuron
params.P_inTOpy = 0.437; % 0.411;    % Probability of connection between In -> Pyramidal
params.P_pyTOpy = 0.243; % 0.16;    % Probability of connection between Py -> Pyramidal
params.P_inTOin = 0.451; % 0.451;   % Probability of connection between In -> Interneuron

% Integrate and Fire Morphological parameters:
params.g_m_P = 25e-9; % Membrane conductance Pyramidal cells
params.g_m_I = 20e-9; % Membrane conductance Inhibitory cells

params.C_P = 0.5e-9; % Membrane conductance Pyramidal cells
params.C_I = 0.2e-9; % Membrane conductance Inhibitory cells

% External input
params.u = mu; % mean input firing rate
params.u_bkg = 1; % Background input firing rate

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
        params.alpha_i = -0.705;%-1.128;%-0.235;% -0.2635;           % Inhibitory -> Pyramidal
        params.alpha_e = 2.472;% 0.2813;            % Pyramidal -> Inhibitory
        params.alpha_ri = -3.316;% -1.25; %-0.2501;  % Inhibitory -> GABAb
        params.alpha_re = 1.255; % 0.2005;%2*0.2005; % Pyramidal -> Pyramidal
        params.alpha_b = -0.32;%-0.128;%-0.5332;% -0.1143;           % GABAb -> Pyramidal %  Note, this is (non-intuituvely) positive during the fit, but should be negative in here.
        params.alpha_u = 1.637;%0.7075; %0.1364;   % GABAb
        params.alpha_uinterneuron = 2.591;
        
        % Probabilities
        params.c_constant = 1;
%         params.P_inTOin = 0.5 * params.P_inTOin;
%         params.P_inTOpy = 0.5 * params.P_inTOpy;
        
        % Time constants:
        params.tau_sb = 10 * params.tau_sp;
        
        % Nonlinearity
        % Pyramidal:
        params.e0 =         37.4;% 37.87; % 25.65; % Maximum firing rates
        params.gompertz.b = 4.817;% 0.0882; % 0.8736;
        params.gompertz.c = 4.072;% 0.3734; % 0.9686;
        params.gompertz.d = 0.102;% 0.09833; % 0.4128;
        
        params.naka.M = 39.35;  % Thalamic input: 42.52;    % No thalamic input: 41.18; <- no bias % exponent: 2 -> 42.34
        params.naka.a = 3.397;  % Thalamic input: 1.988;    % No thalamic input: 2.919; <- no bias % exponent: 2 -> 2
        params.naka.b = 0;      % Thalamic input: -22.13;   % No thalamic input: 0; <- no bias % exponent: 2 -> 7.558
        params.naka.s = 23.94;   % Thalamic input: 17.5;     % No thalamic input: 24.8; <- no bias % exponent: 2 -> 17.12
       
        % Interneurons:
        params.e0i =         63.23;% 63.15; % 65; % Maximum firing rates
        params.gompertzi.b = 4.337;% 0.1725; % 2.381;
        params.gompertzi.c = 3.964;% 0.9534; % 1.79;
        params.gompertzi.d = 0.09888;% 0.09868; % 0.2661;
        
        params.nakai.M = 66.81; % Thalamic input: 71.03;    % No thalamic input: 69.29; <- no bias % exponent: 2 -> 70.6
        params.nakai.a = 3.186; % Thalamic input: 1.903;    % No thalamic input: 2.734; <- no bias % exponent: 2 -> 2
        params.nakai.b = 0;     % Thalamic input: -16.32;   % No thalamic input: 0; <- no bias % exponent: 2 -> 6.957
        params.nakai.s = 23.2; % Thalamic input: 16.35;    % No thalamic input: 23.94; <- no bias % exponent: 2 -> 16.79
       
        
    otherwise
        error('%s rythm not implemented, sorry!', mode);
end
