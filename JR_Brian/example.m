
addpath ../JR/ % Load nonlinearity

N = 5000; % Number of samples: 1 sample = 1 milisecond
u = 5;%1.5;

params = set_parameters('allen', u);

dt = params.dt;
t = 0:dt:(N-1)*dt;

c_constant = 1000;
c1 = 4 * c_constant * params.P_pyTOin; % Excitatory synapses into inhibitory population
c2 = 1 * c_constant * params.P_inTOpy; % Inhibitory synapses into pyramidal population

alpha_i = params.alpha_ie * c2 * 2 * params.e0 * params.decay_i; % lumped constant (inhibitory input to pyramidal cells)
alpha_e = params.alpha_ei * c1 * 2 * params.e0 * params.decay_e; % lumped constant (excitatory input to inhibitory interneurons)

y0 =[0;
    0;
    0;
    0;
    u;
    alpha_i;
    alpha_e;
    ];

[t,y] = ode45(@(t,y) ode(t,y,params,alpha_i,alpha_e, dt), [min(t) max(t)], y0);
% [t,y] = ode23(@(t,y) ode(t,y,params,alpha_i,alpha_e, dt), [min(t) max(t)], y0);
% [t,y] = ode15s(@(t,y) ode(t,y,params,alpha_i,alpha_e, dt), [min(t) max(t)], y0);

figure
plot(t,y(:,[1 3]));
% ylim([-1e5 0.5e6]);
% legend({'x1' 'x2' 'x3' 'x4'});
legend({'x1' 'x3'});

function dy = ode(t,y, params, alpha_i, alpha_e, dt)
    dy = zeros(7,1);
    dy(1) = y(2);
    dy(2) = alpha_i * params.decay_i * non_linear_sigmoid(y(3),params.r, params.v0) - 2*params.decay_i*y(2) - (params.decay_i^2) * y(1);
    dy(3) = y(4);
    dy(4) = alpha_e * params.decay_e * gompertz(y(5)+ y(1), params) - 2*params.decay_e*y(4) - (params.decay_e^2) * y(3);
    dy(5) = y(5);
    dy(6) = y(6);
    dy(7) = y(7);

end
