
addpath ../JR/ % Load nonlinearity

N = 1000000; % Number of samples: 1 sample = 1 milisecond
u = 1.5;%1.5;

params = set_parameters('brunel', u);

t = 0:params.dt:(N-1)*params.dt;

C = 100*0.2;
c1 = 4*C;
c2 = 1*C;

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

% [t,y] = ode45(@(t,y) ode(t,y,params,alpha_i,alpha_e), [min(t) max(t)], y0);
% [t,y] = ode23(@(t,y) ode(t,y,params,alpha_i,alpha_e), [min(t) max(t)], y0);
[t,y] = ode15s(@(t,y) ode(t,y,params,alpha_i,alpha_e), [min(t) max(t)], y0);

figure
plot(t,y(:,[1:4]));
% ylim([-1e5 0.5e6]);
legend({'x1' 'x2' 'x3' 'x4'});

function dy = ode(t,y, params, alpha_i, alpha_e)
    dy = zeros(7,1);
    dy(1) = y(2);
    dy(2) = alpha_i*params.decay_i*phi(y(3),params) - 2*params.decay_i*y(2) - 1*params.decay_i*y(1);
    dy(3) = y(4);
    dy(4) = alpha_e*params.decay_e*phi(y(5)+ y(1), params) - 2*params.decay_e*y(4) - 1*params.decay_e*y(3);
    dy(5) = y(5);
    dy(6) = y(6);
    dy(7) = y(7);

end

function out = phi(x, params)
	out = non_linear_sigmoid(x,params.r, params.v0);
end