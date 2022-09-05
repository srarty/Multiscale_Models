% Set parameters for the JR model
%
function params = set_parameters(mode, varargin)

mu = 20; %11; % Default external input value

if nargin < 1
    mode = 'alpha';
elseif nargin > 1
    % If the function is called with 2 input arguments, the input (u) is 
    % defined by the second argument
    mu = varargin{1};
end

switch mode
    case 'recursive'
        params.e0 =  30; %1.5;%120;% 10; %120;%30;%10;%30; %15; % max firing rate
        params.e0i = 30; %20;%250;% 30; %250;%35;%40;%35; % max firing rate
        
        params.v0 = 19.72; % Firing Threshold (a)
        params.r = 34.34; % Sigmoid slope (b)
        
        params.v0_p = 16; %15; % Firing Threshold (a)
        params.r_p = 10; %8; % Sigmoid slope (b)
        
        params.gompertz.a = 1;
        params.gompertz.b = 3.824;
        params.gompertz.c = -9.142;
        params.gompertz.d = 0.03179;
        
        params.gompertzi.a = 1;
        params.gompertzi.b = 4;%3;
        params.gompertzi.c = 0.2;%1;
        params.gompertzi.d = 0.15;%0.15;
        
        % time constants
        params.tau_mp = 0.009738; % Membrane time constant - Pyramidal (decay)
        params.tau_sp = 0.001472; % Synaptic time constant - AMPA on Py (rising)
        
        params.tau_mi = 0.01966;  % Decay tau Interneurons
        params.tau_si = 0.005608; % Rising tau Inhibitory
        
        params.tau_mrp = 0.01975;  % Decay tau recursive excitation
        params.tau_srp = 0.002665; % Rising tau recursive excitation
        
        params.tau_mri = 0.009229; % Decay tau recursive inhibition
        params.tau_sri = 0.006031; % Rising tau recursive inhibition
        
        % External input
        params.u = mu; %11;%220;%15;%11;1/0.0        % mean input mem potential
        
        % Gains
        params.alpha_ie = -0.5533; %Inhibitory gain into pyramidal (Interneuron -> Py)
        params.alpha_ei = 1.1975; % Excitatory gain into interneuron (Py -> Interneuron) (bifurcation: alpha_ei > 392.6)
        params.alpha_re = 0.4132; % Recursive excitatory gain (bifuration: alpha_re > 4.7)
        params.alpha_ri = -1.441; % Recursive inhibitory gain
        params.alpha_u = 0.06338; % External excitatory gain into pyramidal (U -> Py)
        
        params.c_constant = 1000; % Connectivity constant
        
        params.P_pyTOin = 0.395; %0.1975;%0.395; % Probability of connection between Py -> Interneuron
        params.P_inTOpy = 0.411;%0.411; % Probability of connection between In -> Pyramidal
        params.P_pyTOpy = 0.16; %0.16;  % Probability of connection between Py -> Pyramidal
        params.P_inTOin = 0.451;%0.451;%0.451; %0.250; %0.451; % Probability of connection between In -> Interneuron
        

        params.dt = 0.001;     % sampling time step         
        params.scale = 1;%1e2;% 1e-3; % Scale to fix mismatch in state amplitudes. Not to be confused with the scael in analytic_kalman_filter_2
        
    case 'allen'
        params.e0 = 30; % max firing rate
        params.r = 8; %11;%1;   % Sigmoid slope
        params.v0 = 15; %15;%21; % Firing Threshold
        
        params.gompertz.a = 1; % 70;
        params.gompertz.b = 3; % 14;
        params.gompertz.c = 1; % 14.49;
        params.gompertz.d = 0.15; % 0.3292;
        
        % inverse time constants
        params.decay_e = 270.2703; % bifurcation(u=1) = 76.3601; %270.2703;%268.4672;   % Excitatory synapse (AMPA) into inhibitory interneurons 1/tau_e
        params.decay_i = 106.3830; %106.1121;   % Inhibitory synapse (GABA) into pyramidal cells (1/tau_i)
        
        params.u = mu; %11;%220;%15;%11;1/0.0        % mean input mem potential
        
        ratio = -3.106472;
        params.alpha_ie = -0.479; %-128.8; % Inhibitory gain into pyramidal (Interneuron -> Py)
        params.alpha_ei = ratio * params.alpha_ie;% 1.488; %1011 <- fit result; % Excitatory gain into interneuron (Py -> Interneuron)
        
        params.P_pyTOin = 0.395; % Probability of connection between Py -> Interneuron
        params.P_inTOpy = 0.411; % Probability of connection between Py -> Interneuron
        
        params.dt = 0.001;     % sampling time step         
        params.scale = 1;%1e2;% 1e-3; % Scale to fix mismatch in state amplitudes. Not to be confused with the scael in analytic_kalman_filter_2
        
    case 'brunel'
        params.e0 = 30; % max firing rate
        params.r = 1;%1;   % Sigmoid slope
        params.v0 = 21;%21; % Firing Threshold
        
        params.gompertz.a = 30;
        params.gompertz.b = 14;
        params.gompertz.c = 14.49;
        params.gompertz.d = 0.3292;
        
        % inverse time constants
        params.decay_e = 223.7074;   % Excitatory synapse (AMPA) into inhibitory interneurons 1/tau_e
        params.decay_i = 88.4151;   % Inhibitory synapse (GABA) into pyramidal cells (1/tau_i)
        params.tau_e = 1/params.decay_e;
        params.tau_i = 1/params.decay_i;
        
        params.u = mu; %11;%220;%15;%11;        % mean input mem potential
        params.alpha_ei = 171.2; % Excitatory gain into interneuron (Py -> Interneuron)
        params.alpha_ie = -238.6; % Inhibitory gain into pyramidal (Interneuron -> Py)

        params.P_pyTOin = 0.2; % Probability of connection between Py -> Interneuron
        params.P_inTOpy = 0.2; % Probability of connection between Py -> Interneuron

        
        params.dt = 0.001;     % sampling time step         
        params.scale =  1;%1e-3; % Scale to fix mismatch in state amplitudes. Not to be confused with the scael in analytic_kalman_filter_2
        
    case 'alpha' % 8 Hz - 12 Hz
        params.e0 = 2.5;  % max firing rate
        params.r = 3.0285;  % erf sigmoid slope % params.r = 0.56;  % logistic sigmoid slope
        params.v0 = 6; % Firing Threshold
        % inverse time constants
        params.decay_e = 50; % 100;% (1/ tau_e)
        params.decay_i = 100; % 50; % (1/tau_i)
        params.alpha_ei = 3.25;% 3.25;     % Gains (a_ei = excitatory), lumped parameter will look like: % alpha_i = 162500
        params.alpha_ie = 6.25;%22;%12.5;  % (a_ie = inhibitory), % alpha_e = 440000
        params.u = mu;%11;%220;%15;%11;        % mean input mem potential
        params.dt = 0.001;     % sampling time step   
        params.scale = 1; % Scale to fix mismatch in state amplitudes. Not to be confused with the scael in analytic_kalman_filter_2
        
    case 'beta' % 12 Hz - 30 Hz
        params.e0 = 2.5;  % max firing rate
        params.r = 3.0285;  % erf sigmoid slope
        params.v0 = 6; % Firing Threshold
        % inverse time constants
        params.decay_e = 150;% (1/ tau_e)
        params.decay_i = 150; % (1/tau_i)        
        params.alpha_ei = 22;% 3.25;     % Gains (a_ei = excitatory), lumped parameter will look like: % alpha_i = 162500
        params.alpha_ie = 12;%6.25;%6.25;%22;%12.5;  % (a_ie = inhibitory), % alpha_e = 440000
        params.u = mu;%11;%220;%15;%11;        % mean input mem potential
        params.dt = 0.001;     % sampling time step           
        params.scale = 1; % Scale to fix mismatch in state amplitudes. Not to be confused with the scael in analytic_kalman_filter_2
        
    case 'gamma' % > 30 Hz
        params.e0 = 2.5;  % max firing rate
        params.r = 3.0285;  % erf sigmoid slope
        params.v0 = 6; % Firing Threshold
        % inverse time constants
        params.decay_e = 150*2;% (1/ tau_e)
        params.decay_i = 150*4; % (1/tau_i)
        params.alpha_ei = 22;% 3.25;     % Gains (a_ei = excitatory), lumped parameter will look like: % alpha_i = 162500
        params.alpha_ie = 12;%6.25;%6.25;%22;%12.5;  % (a_ie = inhibitory), % alpha_e = 440000
        params.u = mu;%11;%220;%15;%11;        % mean input mem potential
        params.dt = 0.001;     % sampling time step
        params.scale = 1; % Scale to fix mismatch in state amplitudes. Not to be confused with the scael in analytic_kalman_filter_2
        
    case 'delta' % 1 Hz - 4 Hz
        params.e0 = 2.5;  % max firing rate
        params.r = 3.0285;  % erf sigmoid slope
        params.v0 = 6; % Firing Threshold
        % inverse time constants
        params.decay_e = 100/4;% (1/ tau_e)
        params.decay_i = 50/3; % (1/tau_i)
        params.alpha_ei = 3.25;% 3.25;     % Gains (a_ei = excitatory), lumped parameter will look like: % alpha_i = 162500
        params.alpha_ie = 22;%6.25;%6.25;%22;%12.5;  % (a_ie = inhibitory), % alpha_e = 440000
        params.u = mu;%11;%220;%15;%11;        % mean input mem potential
        params.dt = 0.001;     % sampling time step
        params.scale = 1; % Scale to fix mismatch in state amplitudes. Not to be confused with the scael in analytic_kalman_filter_2
    otherwise
        error('%s rythm not implemented, sorry!', mode);
end
