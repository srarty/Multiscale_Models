% Set parameters for the JR model
%
function params = set_parameters(mode, varargin)

mu = 1; % Default external input value

if nargin < 1
    mode = 'recursive';
elseif nargin > 1
    % If the function is called with 2 input arguments, the input (u) is 
    % defined by the second argument
    mu = varargin{1};
end

switch mode
    case 'recursive'
        % Maximum firing rates
        params.e0 =  27.81; %25.87; % Pyramidal
        params.e0i = 60;%58.85; % % Inhibitory interneurons
        
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
        params.gompertz.b = 1.512; %1.555;
        params.gompertz.c = 1.508; %1.245;
        params.gompertz.d = 0.2805; %0.3065;
        % Interneurons:
        params.gompertzi.b = 2.124; % 2.208; %2.291; %
        params.gompertzi.c = 1.738; % 1.667; %1.551; %
        params.gompertzi.d = 0.1964; % 0.209; %0.2315; %
        
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
        params.alpha_i = -0.5269;       %Inhibitory gain into pyramidal (Interneuron -> Py)
        params.alpha_e = 1.124;         % Excitatory gain into interneuron (Py -> Interneuron) (bifurcation: alpha_ei > 392.6)
        params.alpha_re = 0.4009;       % Recursive excitatory gain (bifuration: alpha_re > 4.7)
        params.alpha_ri = -2.5;%-0.9698;% Recursive inhibitory gain, -2.5 increases Py resting membrane potential
        params.alpha_u = 0.0615;        % External excitatory gain into pyramidal (U -> Py)
        
        % Connectivity parameters:
        params.c_constant = 1; % Connectivity constant
        params.P_pyTOin = 0.395; 	% Probability of connection between Py -> Interneuron
        params.P_inTOpy = 0.411;    % Probability of connection between In -> Pyramidal
        params.P_pyTOpy = 0.16;    % Probability of connection between Py -> Pyramidal
        params.P_inTOin = 0.451;   % Probability of connection between In -> Interneuron
        
        % Integrate and Fire Morphological parameters:
        params.g_m_P = 25e-9; % Membrane conductance Pyramidal cells
        params.g_m_I = 20e-9; % Membrane conductance Inhibitory cells
        
        % External input
        params.u = mu; % mean input firing rate
        
        params.dt = 0.001;     % sampling time step         
        params.scale = 1; % Scale to fix mismatch in state amplitudes. Not to be confused with the scale in analytic_kalman_filter_2
    
    otherwise
        error('%s rythm not implemented, sorry!', mode);
end
