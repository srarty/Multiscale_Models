% Same as NMM_diff_equations_DblExp.m but with recursive connections in the
% Pyramidal and Ihibitory populations.
%
% Spike train synchrony (Alex): pubmed.ncbi.nlm.nih.gov/17628690/
%
% Machine learning to find ideal parameters for NMM taht match the LIF.
% Having one NMM for each different LIF.
function [x, y, t, f_e, f_i] = NMM_diff_equations_DblExp_recursive(varargin)
    clear option option2
    if nargin >= 2
        option = varargin{1};
        value = varargin{2};
    end
    if nargin == 4
        option2 = varargin{3};
        value2 = varargin{4};
    end
    
    N = 2000; % Number of samples: 1 sample = 1 milisecond
    u = 0;

    params = set_parameters('recursive', u);
    
    % Options
    params.options.ADD_NOISE = 0; % External input noise (0 = no noise, 1 = noise)
    params.options.CHANGE_U = 0; % 0: U doesn't change during simulation. Anyother value of CHANGE_U: U changes.
    
    INPUT_CURRENT_PY = 1000 * 0e-12 / params.g_m_P;
    INPUT_CURRENT_IN = 1000 * 0e-12 / params.g_m_I;
    % Parse inputs --------------------------------------------------------
    if exist('option','var')
        try
            params.(option) = value;
        catch
            error(['Couldn''t assign value: ' num2str(value) ' to the parameter: ' option]);
        end
    end
    
    if exist('option2','var')
        try
            params.(option2) = value2;
        catch
            error(['Couldn''t assign value: ' num2str(value2) ' to the parameter: ' option2]);
        end
    end
    % --------------------------------------------------- End input parsing
    
    % Phi (nonlinearity) functions 
    b = params.gompertz.b;
    c = params.gompertz.c;
    d = params.gompertz.d;
    
    ib = params.gompertzi.b;
    ic = params.gompertzi.c;
    id = params.gompertzi.d;
    
    v0 = params.v0;
    r = params.r;
    e_0 = params.e0; 
    e_0i = params.e0i;
    
    % Synaptic functions
    S1 = @(x) gompertz_io(x + INPUT_CURRENT_IN, e_0i, ib, ic ,id);% sigmoid_io(x, e_0i, v0, r); %
    S2 = @(x) gompertz_io(x + INPUT_CURRENT_PY, e_0, b, c, d);  % sigmoid_io(x, e_0, v0, r); % 
        
    dt = params.dt;
    t = 0:dt:(N-1)*dt;

    alpha_i = params.alpha_i;
    alpha_e = params.alpha_e; 
    alpha_re = params.alpha_re; 
    alpha_ri = params.alpha_ri; 
    alpha_u = params.alpha_u; 

    % Steady state
    eq =   [-3.8824;
            -197.4781;
            0.7247;
            74.4335;
            0.2090;
            10.5810;
            -2.5466;
            -275.9404;
            0;
            0;];
        
%     eq =[0;
%          0;
%          0;
%          0;
%          0;
%          0;
%          0;
%          0;
%          0;
%          0;];
        
    x0 =[eq(1);
        eq(2);
        eq(3);
        eq(4);
        eq(5);
        eq(6);
        eq(7);
        eq(8);
        eq(9);
        eq(10);
        alpha_i;
        alpha_e;
        alpha_re;
        alpha_ri;
        alpha_u;
        ];
    
    [t,x,y] = ode45(@(t,x) ode(t,x,params,dt, S1, S2), t, x0); % use "t", instead of "[min(t) max(t)]" fix the output time vector
    %     [t,x,y] = ode23(@(t,x) ode(t,x,params,dt), [min(t) max(t)], x0);
    %     [t,x,y] = ode113(@(t,x) ode(t,x,params,dt), [min(t) max(t)], x0);
    
    for i = 1:size(x,1)
        y(i) = x(i,1) + x(i,9) + x(i,5) + INPUT_CURRENT_PY;
    end
    
    % Calculate firing rate
    f_i = S1(x(:,3) + x(:,7)); 
    f_e = S2(x(:,1) + x(:,5) + x(:,9));
    
    if nargin == 0        
%         close all
        do_plot(x,t,y, f_e, f_i);
    end
end

function dx = ode(t,x,params,dt, S1, S2)
    
    u = params.u;
    
    % Following lines are meant to change the input mid simulation, comment
    % them to run it with constant input.
    if isfield(params,'options') & isfield(params.options,'CHANGE_U') & params.options.CHANGE_U
        if t >= 1.5
            u = 0.5;
        elseif t >= 1
            u = 0.75;
        elseif t >= 0.5
            u = 0.25;
        end
    end

    Tau_coeff = @(m, s) 1/(m*s);% Nicola Campbell
    
    c_constant = params.c_constant;
    c1 = 70.5575    * c_constant * params.P_pyTOin; % Excitatory synapses into inhibitory population
    c2 = 29.7217    * c_constant * params.P_inTOpy; % Inhibitory synapses into pyramidal population
    c3 = 141.1007   * c_constant * params.P_pyTOpy; % Recursive excitation to pyramidal cells
    c4 = 9.5339     * c_constant * params.P_inTOin; % Recursive inhibition to inhibitory cells
    c5 = 17.8288    * c_constant;                   % External excitatory synapses into pyramidal population
    
    tau_sp = params.tau_sp;
    tau_mp = params.tau_mp;
    
    tau_si = params.tau_si;
    tau_mi = params.tau_mi;
    
    tau_srp = params.tau_srp;
    tau_mrp = params.tau_mrp;
    
    tau_sri = params.tau_sri;
    tau_mri = params.tau_mri;
    
    % Lumped gain parameters
    AmplitudeI  = c2 * x(11) * Tau_coeff(tau_mp,  tau_sp);
    AmplitudeE  = c1 * x(12) * Tau_coeff(tau_mi,  tau_si);
    AmplitudeRE = c3 * x(13) * Tau_coeff(tau_mrp, tau_srp);
    AmplitudeRI = c4 * x(14) * Tau_coeff(tau_mri, tau_sri);
    AmplitudeU  = c5 * x(15) * Tau_coeff(tau_mrp, tau_srp);
    
    % Diff equations ------------------------------------------------------
    dx = zeros(15,1);

    % Double exponential from Nicola-Campbell (2013):
    % I->P
    dx(1) = x(2) - x(1)/tau_mp;
    dx(2) = - x(2)/tau_sp + AmplitudeI * S1(x(3) + x(7));
    % P->I
    dx(3) = x(4) - x(3)/tau_mi;
    dx(4) = - x(4)/tau_si + AmplitudeE * S2(x(1) + x(5) + x(9));
    % Recurrent Pyramidal P->P
    dx(5) = x(6) - x(5)/tau_mrp;
    dx(6) = - x(6)/tau_srp + AmplitudeRE * S2(x(1) + x(5) + x(9));
    % Recurrent Inhibition I->I
    dx(7) = x(8) - x(7)/tau_mri;                       
    dx(8) = - x(8)/tau_sri + AmplitudeRI * S1(x(3) + x(7));    
    % External input u->P
    dx(9) = x(10) - x(9)/tau_mrp;
    dx(10) = - x(10)/tau_srp + AmplitudeU * (u + (params.options.ADD_NOISE * (sqrt(u).*randn(1,1))));
    % Parameters:
    dx(11) = 0; % alpha_i
    dx(12) = 0; % alpha_e
    dx(13) = 0; % alpha_re
    dx(14) = 0; % alpha_ri
    dx(15) = 0; % alpha_u

end

function do_plot(x,t, Vm, f_e, f_i)    
    figure
    plot(t,x(:,[1 3]));
    legend({'x1' 'x3'});
    ylabel('mV');
    xlabel('Time (s)');
    
    figure
    plot(t, x(:,1) + x(:,5) + x(:,9));
    hold
    plot(t, x(:,3) + x(:,7));
    legend({'V_{mp}' 'V_{mi}'});
    ylabel('mV');
    xlabel('Time (s)');

    figure
    plot(t, Vm, 'k');
    title('Vm_{Py}');
    ylabel('mV');
    xlabel('Time (s)');
    
    figure
    plot(t, f_e); hold on;
    plot(t, f_i);
    ylabel('Spike rate (Hz)');
    xlabel('Time (s)');
    
    figure
    plot(t, x(:,9));    
    
    figure
    subplot(1,2,1)
    plot(x(:,1), x(:,3));
    subplot(1,2,2)
    plot(x(:,1), x(:,2));
end

function out = sigmoid_io(x, e0, v0, r)
    out = e0 * (0.5 * erf((x - v0) / (sqrt(2)*r)) + 0.5);
end

function out = gompertz_io(x, a, b, c, d)
    out = a * exp(-b*exp(-d*(x-c)));
end