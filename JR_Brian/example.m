close all

addpath ../JR/ % Load nonlinearity

N = 3000; % Number of samples: 1 sample = 1 milisecond
u = 6;%1.5;

params = set_parameters('allen', u);

dt = params.dt;
t = 0:dt:(N-1)*dt;

c_constant = 1000;
c1 = 4 * c_constant * params.P_pyTOin; % Excitatory synapses into inhibitory population
c2 = 1 * c_constant * params.P_inTOpy; % Inhibitory synapses into pyramidal population

alpha_i = params.alpha_ie;% * c2 * 2 * params.e0 * params.decay_i; % lumped constant (inhibitory input to pyramidal cells)
alpha_e = params.alpha_ei;% * c1 * 2 * params.e0 * params.decay_e; % lumped constant (excitatory input to inhibitory interneurons)

tau_e = 1/params.decay_e;
tau_i = 1/params.decay_i;


x0 =[0;
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
    
    u = params.u;
    
    % Synaptic functions
    S1 = @(x3) 0.5*erf((x3 - v0) / (sqrt(2)*r)) + 0.5;               
    S2 = @(x1) exp(-b*exp(-d*(x1+u-c)));                              
    
    dx = zeros(7,1);
    dx(1) = x(1) + x(2);                                                  
    dx(2) = x(2) - 2*x(2)/tau_i - x(1)/tau_i^2 + 2.6167e6*alpha_i*S1(x(3));
    dx(3) = x(3) + x(4);                                                  
    dx(4) = x(4) + - 2*x(4)/tau_e - x(3)/tau_e^2 + 2.5450e7*alpha_e*S2(x(1));
    dx(5) = x(5);
    dx(6) = x(6);
    dx(7) = x(7);
end

function do_plot(y,t, Vm)    
    figure
    plot(t,y(:,[1 3]));
    legend({'x1' 'x3'});

    figure
    plot(y(:,1), y(:,2)); hold on
    plot(y(:,3), y(:,4));
    xlabel('v');
    ylabel('z');
    legend({'Py' 'Inh'});
    
    figure
    plot(t, Vm, 'k');
end