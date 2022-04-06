%% Nullclines

%% Fetch parameters
u = 11;
params = set_parameters('allen', u);

% Population parameters
alpha_e = params.alpha_ei;
alpha_i = params.alpha_ie;
tau_e = 1/params.decay_e;
tau_i = 1/params.decay_i;

% Gompertz nonlinearity parameters
b = params.gompertz.b;
c = params.gompertz.c;
d = params.gompertz.d;

% Nonlinear error function parameters
v0 = params.v0;
r = params.r;
e_0 = params.e0; 
u = params.u;

% Synaptic functions
S1 = @(x3) 0.5*erf((x3 - v0) / (sqrt(2)*r)) + 0.5;               
S2 = @(x1) exp(-b*exp(-d*(x1+u-c)));                              

% Connectivity constants
c_constant = 1000;
c1 = 4 * c_constant * params.P_pyTOin; % Excitatory synapses into inhibitory population
c2 = 1 * c_constant * params.P_inTOpy; % Inhibitory synapses into pyramidal population
lump_i = c2 * 2 * e_0 * params.decay_i; %2.6167e6;
lump_e = c1 * 2 * e_0 * params.decay_e; %2.5450e7;                            

% Time vector
dt = params.dt;

%% Define nullclines
nullclineX1 = @(x3) tau_i^2 * (lump_i*alpha_i*S1(x3)); 
nullclineX3 = @(x1) tau_e^2 * (lump_e*alpha_e*S2(x1)); 

%% Plot
x1 = [-100:0.01:0];
x3 = [0:0.01:100];
figure
plot(x1,nullclineX3(x1));hold
plot(nullclineX1(x3),x3);



