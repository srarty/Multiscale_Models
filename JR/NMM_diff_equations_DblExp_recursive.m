% Same as NMM_diff_equations_DblExp.m but with recursive connections in the
% Pyramidal and Ihibitory populations.
function [x, y, t] = NMM_diff_equations_DblExp_recursive(varargin)
    option = ''; value = [];
    if nargin >= 2
        option = varargin{1};
        value = varargin{2};
    end
    if nargin == 4
        option2 = varargin{3};
        value2 = varargin{4};
    end
    
% 	close all

    N = 2000; % Number of samples: 1 sample = 1 milisecond
    u = 9;%20; % 1.5;

    params = set_parameters('allen', u);
    % Parse inputs --------------------------------------------------------
    switch option
        case "AmplitudeE"
            params.AmplitudeE = value;
        case "AmplitudeI"
            params.AmplitudeI = value;
        case "tau_s_e"
            params.tau_s_e = value;
        case "tau_s_i"
            params.tau_s_i = value;
        case "tau_m_e"
            params.tau_m_e = value;
        case "tau_m_i"
            params.tau_m_i = value;
        case "u"
            params.u = value;
    end
    if exist('option2','var')
        switch option2
            case "AmplitudeE"
                params.AmplitudeE = value2;
            case "AmplitudeI"
                params.AmplitudeI = value2;
            case "tau_s_e"
                params.tau_s_e = value2;
            case "tau_s_i"
                params.tau_s_i = value2;
            case "tau_m_e"
                params.tau_m_e = value2;
            case "tau_m_i"
                params.tau_m_i = value2;
        end
    end
    % --------------------------------------------------- End input parsing
    
    dt = params.dt;
    t = 0:dt:(N-1)*dt;

    alpha_i = params.alpha_ie;% * c2 * 2 * params.e0 * params.decay_i; % lumped constant (inhibitory input to pyramidal cells)
    alpha_e = params.alpha_ei;% * c1 * 2 * params.e0 * params.decay_e; % lumped constant (excitatory input to inhibitory interneurons)

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
    % [t,x,y] = dde23(@(t,x) ode(t,x,params, dt), [min(t) max(t)], x0);
    
    for i = 1:size(x,1)
        y(i) = x(i,1) + u + x(i,5);
    end
    
    do_plot(x,t,y)
end

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
    c3 = 4 * c_constant * 0.16;
    c4 = 1 * c_constant * 0.125; % Allen(LIF) = 0.451; % Same oscillatory freq (with u=9) = 0.125
    lump_i = c2 * 2 * e_0 * params.decay_i; %2.6167e6;
    lump_e = c1 * 2 * e_0 * params.decay_e; %2.5450e7;
    
    tau_s_e = 0.001227; % 0.001001;
    tau_m_e = 0.008115; % 0.008339;
    
    tau_s_i = 0.004674; % 0.004378;
    tau_m_i = 0.01638; % 0.01668;
    
    tau_s_re = 0.002221;
    tau_m_re = 0.01646;
    
    tau_s_ri = 0.005021;
    tau_m_ri = 0.007698;
    
    alpha_e = 2.395; % 2.25;
    alpha_i = -1.107; % -1.054;
    
    alpha_re = 0.8264;
    alpha_ri = -2.869;
    
%     j_e = 14e-12;
%     j_i = -74e-12;

    % Parse the optional inputs -------------------------------------------
    if ~isfield(params,'AmplitudeI')
        AmplitudeI = c2 * 2 * e_0 * alpha_i * (1/(tau_m_i + tau_s_i));
    else
        AmplitudeI = params.AmplitudeI;
    end
    
    if ~isfield(params,'AmplitudeE')
        AmplitudeE = c1 * 2 * e_0 * alpha_e * (1/(tau_m_e + tau_s_e));
    else
        AmplitudeE = params.AmplitudeE;
    end
    
    AmplitudeRE = c3 * 2 * e_0 * alpha_re * (1/(tau_m_re + tau_s_re));
    AmplitudeRI = c4 * 2 * e_0 * alpha_ri * (1/(tau_m_ri + tau_s_ri));
    
    if isfield(params,'tau_s_e'), tau_s_e = params.tau_s_e; end
    if isfield(params,'tau_s_i'), tau_s_i = params.tau_s_i; end
    if isfield(params,'tau_m_e'), tau_m_e = params.tau_m_e; end
    if isfield(params,'tau_m_i'), tau_m_i = params.tau_m_i; end    
    % -------------------------------------------------- End parsing inputs

    
    % Diff equations ------------------------------------------------------
    dx = zeros(9,1);
    
%     dx(1) = x(1) - x(1)/tau_m_i - x(2)/tau_m_i;
%     dx(2) = x(2) - x(2)/tau_s_i + c2 * e_0 * alpha_i * (tau_m_i/(tau_m_i - tau_s_i)) * S1(x(3));
%     dx(3) = x(3) - x(3)/tau_m_e - x(4)/tau_m_e;
%     dx(4) = x(4) - x(4)/tau_s_e + c1 * e_0 * alpha_e * (tau_m_e/(tau_m_e - tau_s_e)) * S2(x(1));
      
%     dx(1) = x(1) + (-x(1) -x(2)/25e-9)/tau_m_i; % see if dt is needed
%     dx(2) = x(2) - j_i * x(2) / tau_s_i + c2 * e_0 * alpha_i * S1(x(3));
%     dx(3) = x(3) + (-x(3) -x(4)/20e-9)/tau_m_e;
%     dx(4) = x(4) - j_e * x(4) / tau_s_e + c1 * e_0 * alpha_e * S2(x(1));

    % Double exponential from Nicola-Campbell (2013):
    % TODO: Check the coefficients of the convolution, specifically 1/(tau_m+tau_s)
    % I->P
    dx(1) = x(1) + x(2) - x(1)/tau_m_i;
    dx(2) = x(2) - x(2)/tau_s_i + AmplitudeI * S1(x(3) + x(7));
    % P->I
    dx(3) = x(3) + x(4) - x(3)/tau_m_e;
    dx(4) = x(4) - x(4)/tau_s_e + AmplitudeE * S2(x(1) + x(5));    
    % Recurrent Pyramidal P->P
    dx(5) = x(5) + x(6) - x(5)/tau_m_re;                        
    dx(6) = x(6) - x(6)/tau_s_re + AmplitudeRE * S2(x(1) + x(5));    
    % Recurrent Inhibition I->I
    dx(7) = x(7) + x(8) - x(7)/tau_m_ri;                       
    dx(8) = x(8) - x(8)/tau_s_ri + AmplitudeRI * S1(x(3) + x(7));    
    % Parameters:
    dx(9) = x(9);
    dx(10) = x(10);
    dx(11) = x(11);
    
%     % Recurrent
%     % TODO
%     dx(1) = x(1) + x(2) - x(1)/tau_m_i;
%     dx(2) = x(2) - x(2)/tau_s_i + AmplitudeI * S1(x(3));
%     dx(2) = x(2) - x(2)/tau_s_i + AmplitudeI * S1(x(3));
%     dx(3) = x(3) + x(4) - x(3)/tau_m_e;
%     dx(4) = x(4) - x(4)/tau_s_e + AmplitudeE * S2(x(1));
%     dx(4) = x(4) - x(4)/tau_s_e + AmplitudeE * S2(x(1));    
%     dx(5) = x(5);
%     dx(6) = x(6);
%     dx(7) = x(7);

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
