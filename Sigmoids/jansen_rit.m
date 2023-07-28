% Jansen and Rit Neural Mass Model

close all

% addpath ../JR/ % Load nonlinearity

% Parameters
A = 3.25;   % Excitatory time constant [ms]
B = 22;     % Inhibitory time constant [ms]
a = 100;    % ms^-1
b = 50;     % ms^-1
C = 135;    % Connectivity scaling parameter [ms^-1]

I = 3.1;    % External input [mV]

% Simulation parameters
dt = 1e-3;      % Time step [ms]
t = 0:dt:1000;   % Time span [ms]

% Generate noise signal
noise_mean = 0;         % Mean of the noise [mV]
noise_amplitude = 0.1;  % Amplitude of the noise [mV]
noise_sigma = 0.1;      % Standard deviation of the noise [mV]
noise = noise_mean + noise_amplitude * randn(size(t));

% Apply noise to input current
I = I + noise_sigma * noise;

% Connectivity constants
C1 = C;
C2 = 0.8 * C;
C3 = 0.25 * C;
C4 = 0.25 * C;

alpha_ep = A * C1 * 2 * e0 * a;
alpha_pe = A * C2 * 2 * e0 * a;
alpha_ip = A * C3 * 2 * e0 * b;
alpha_pi = B * C4 * 2 * e0 * a;

x0 =[0;
    0;
    0;
    0;
    0;
    0;
    0;
    0;
    u;
    alpha_i;
    alpha_e;
    ];

[t,x,y] = ode45(@(t,x) ode(t,x,params, dt), [min(t) max(t)], x0);
for i = 1:size(x,1)
    y(i) = x(i,1) + u;
end
do_plot(x,t,y)


function dx = ode(t,x,params,dt)
    alpha_e = params.alpha_ei;
    alpha_i = params.alpha_ie;
    tau_e = 1/params.decay_e;
    tau_i = 1/params.decay_i;
    
    b = params.gompertz.b;
    c = params.gompertz.c;
    d = params.gompertz.d;
    
    v0 = params.v0;
    r = params.r;
    e_0 = params.e0; 
    u = params.u;
    
    % Synaptic functions
    S1 = @(x3) 0.5*erf((x3 - v0) / (sqrt(2)*r)) + 0.5;               
    S2 = @(x1) exp(-b*exp(-d*(x1+u-c)));                              
    
    c_constant = 1000;
    c1 = 4 * c_constant * params.P_pyTOin; % Excitatory synapses into inhibitory population
    c2 = 1 * c_constant * params.P_inTOpy; % Inhibitory synapses into pyramidal population
    lump_i = c2 * 2 * e_0 * params.decay_i; %2.6167e6;
    lump_e = c1 * 2 * e_0 * params.decay_e; %2.5450e7;
    
    dx = zeros(7,1);
    dx(1) = x(2);
    dx(2) = - 2*x(2)/tau_i - x(1)/tau_i^2 + lump_i*alpha_i*S1(x(3));
    dx(3) = x(4);
    dx(4) = - 2*x(4)/tau_e - x(3)/tau_e^2 + lump_e*alpha_e*S2(x(1));
    dx(5) = 0;
    dx(6) = 0;
    dx(7) = 0;
end