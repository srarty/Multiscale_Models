% close all

addpath ../JR/ % Load nonlinearity

N = 3000; % Number of samples: 1 sample = 1 milisecond
u = 20; % 1.5;

params = set_parameters('allen', u);

dt = params.dt;
t = 0:dt:(N-1)*dt;

alpha_i = params.alpha_ie;% * c2 * 2 * params.e0 * params.decay_i; % lumped constant (inhibitory input to pyramidal cells)
alpha_e = params.alpha_ei;% * c1 * 2 * params.e0 * params.decay_e; % lumped constant (excitatory input to inhibitory interneurons)

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
    
    tau_r_e = 0.001001;
    tau_d_e = 0.008339;
    tau_m_e = 0.001;
    
    tau_r_i = 0.004378;
    tau_d_i = 0.01668;
    tau_m_i = 0.002;
    
    alpha_e = 2.25;
    alpha_i = -1.054;
    
    dx = zeros(7,1);
    dx(1) = x(1) + x(2) - x(1)/tau_r_i;
    dx(2) = x(2) - x(2)/tau_d_i + c2 * 2 * e_0 * alpha_i * (1/(tau_d_i + tau_r_i)) * S1(x(3));
    dx(3) = x(3) + x(4) - x(3)/tau_r_e;
    dx(4) = x(4) - x(4)/tau_d_e + c1 * 2 * e_0 * alpha_e * (1/(tau_d_e + tau_r_e)) * S2(x(1));
    dx(5) = x(5);
    dx(6) = x(6);
    dx(7) = x(7);
end

%     dx(2) = x(2) - x(2)/(tau_d_i) - x(1)/tau_r_i^2 + c2 * 2 * e_0 * alpha_i * (1/(tau_d_i)) * S1(x(3));       
%     dx(4) = x(4) - x(4)/(tau_d_e) - x(3)/tau_r_e^2 + c1 * 2 * e_0 * alpha_e * (1/(tau_d_e)) * S2(x(1));

function do_plot(y,t, Vm)    
    figure
    plot(t,y(:,[1 3]));
    legend({'x1' 'x3'});
    ylabel('mV');
    xlabel('time (s)');
    
    figure
    plot(y(:,1),y(:,3));
    xlabel('x1');
    ylabel('x3');

    figure
    plot(t, Vm, 'k');
    title('Vm_{Py}');
    ylabel('mV');
    xlabel('time (s)');
end
