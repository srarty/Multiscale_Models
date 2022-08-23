% Same as NMM_diff_equations_DblExp.m but with recursive connections in the
% Pyramidal and Ihibitory populations.
%
% Alex's: pubmed.ncbi.nlm.nih.gov/17628690/
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
    
    N = 4000; % Number of samples: 1 sample = 1 milisecond
    u = 0; %20; % 1.5;

    params = set_parameters('recursive', u);
    
    % Options
    params.options.ADD_NOISE = 0; % External input noise (0 = no noise, 1 = noise)
    params.options.CHANGE_U = 1; % 0: U doesn't change during simulation. Anyother value of CHANGE_U: U changes.
    
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
    
    dt = params.dt;
    t = 0:dt:(N-1)*dt;

    alpha_i = params.alpha_ie;
    alpha_e = params.alpha_ei; 
    alpha_re = params.alpha_re; 
    alpha_ri = params.alpha_ri; 

% %     Equilibrium point (calculated with matcont, u=0)
%     eq =   [-3.23168169898182;
%             -197.294592134328;
%             0.466141540947026;
%             57.4505344246254;
%             0.119624560623965;
%             7.26758948851708;
%             -0.766353825579389;
%             -99.5524432048123;];
        
    eq =[0;
         0;
         0;
         0;
         0;
         0;
         0;
         0;];
        
    x0 =[eq(1);
        eq(2);
        eq(3);
        eq(4);
        eq(5);
        eq(6);
        eq(7);
        eq(8);
        params.u;
        alpha_i;
        alpha_e;
        alpha_re;
        alpha_ri;
        ];
    
    [t,x,y] = ode45(@(t,x) ode(t,x,params,dt), [min(t) max(t)], x0);
%     [t,x,y] = ode23(@(t,x) ode(t,x,params,dt), [min(t) max(t)], x0);
%     [t,x,y] = ode113(@(t,x) ode(t,x,params,dt), [min(t) max(t)], x0);
    % Try stochastic solver instead. euler murayama
    for i = 1:size(x,1)
        y(i) = x(i,1) + x(i,9) + x(i,5);
    end
    
    % Calculate firing rate
    f_i = params.e0i * sigmoid_io(x(:,3) + x(:,7), params.v0, params.r);
    f_e = params.e0 * gompertz_io(x(:,1) + x(:,5) + x(:,9), params.gompertz.b, params.gompertz.c, params.gompertz.d);
    
    if nargin == 0        
        close all
        do_plot(x,t,y, f_e, f_i);
    end
end

function dx = ode(t,x,params,dt)
    b = params.gompertz.b;
    c = params.gompertz.c;
    d = params.gompertz.d;
    
    v0 = params.v0;
    r = params.r;
    e_0 = params.e0; 
    e_0i = params.e0i; 
    u = params.u;
    
    % Following lines are meant to change the input mid simulation, comment
    % them to run it with constant input.
    if isfield(params,'options') & isfield(params.options,'CHANGE_U') & params.options.CHANGE_U
        if t >= 3000 * 1e-3
            u = 5;
        elseif t >= 2000 * 1e-3
            u = 3;
        elseif t >= 1000 * 1e-3
            u = 1;
        end
        x(9) = u;
    end
    
    % Synaptic functions
    S1 = @(x) sigmoid_io(x, v0, r);
    S2 = @(x) sigmoid_io(x, v0, r);%gompertz_io(x, b, c, d);   
    Tau_coeff = @(m, s) 1/(m+s); % 1/(m+s)
    
    c_constant = 2000;%params.c_constant; %1000;
    c1 = 4 * c_constant * params.P_pyTOin; % Excitatory synapses into inhibitory population
    c2 = 1 * c_constant * params.P_inTOpy; % Inhibitory synapses into pyramidal population
    c3 = 4 * c_constant * params.P_pyTOpy;
    c4 = 1 * c_constant * params.P_inTOin; % Allen(LIF) = 0.451; % Same oscillatory freq (with u=9) = 0.125; % Same osc freq (u=20)=0.08; % (u=14)=0.06;
    
    tau_se = params.tau_se;
    tau_me = params.tau_me;
    
    tau_si = params.tau_si;
    tau_mi = params.tau_mi;
    
    tau_sre = params.tau_sre;
    tau_mre = params.tau_mre;
    
    tau_sri = params.tau_sri;
    tau_mri = params.tau_mri;
    
    % Lumped parameters
    AmplitudeI  = c2 * 2 * e_0i * x(10) * Tau_coeff(tau_mi,  tau_si);
    AmplitudeE  = c1 * 2 * e_0  * x(11) * Tau_coeff(tau_me,  tau_se);
    AmplitudeRE = c3 * 2 * e_0  * x(12) * Tau_coeff(tau_mre, tau_sre);
    AmplitudeRI = c4 * 2 * e_0i * x(13) * Tau_coeff(tau_mri, tau_sri);
    
    % Diff equations ------------------------------------------------------
    dx = zeros(9,1);

    % Double exponential from Nicola-Campbell (2013):
    % TODO: Check the coefficients of the convolution, specifically 1/(tau_m+tau_s)
    % I->P
    dx(1) = x(2) - x(1)/tau_mi;
    dx(2) = - x(2)/tau_si + AmplitudeI * S1(x(3) + x(7));
    % P->I
    dx(3) = x(4) - x(3)/tau_me;
    dx(4) = - x(4)/tau_se + AmplitudeE * S2(x(1) + x(5) + x(9));
    % Recurrent Pyramidal P->P
    dx(5) = x(6) - x(5)/tau_mre;
    dx(6) = - x(6)/tau_sre + AmplitudeRE * S2(x(1) + x(5) + x(9));
    % Recurrent Inhibition I->I
    dx(7) = x(8) - x(7)/tau_mri;                       
    dx(8) = - x(8)/tau_sri + AmplitudeRI * S1(x(3) + x(7));    
    % Parameters:
    dx(9) = 0;%-x(9) + u + (params.options.ADD_NOISE * (sqrt(u).*randn(1,1))); % Random number % 1.1; % <- steady increase of 1.1 spike/ms/cell/s
    dx(10) = 0; % alpha_i
    dx(11) = 0; % alpha_e
    dx(12) = 0; % alpha_re
    dx(13) = 0; % alpha_ri

end

function do_plot(x,t, Vm, f_e, f_i)    
    figure
    plot(t,x(:,[1 3]));
    legend({'x1' 'x3'});
    ylabel('mV');
    xlabel('Time (s)');
    
    figure
    plot(x(:,1),x(:,3));
    xlabel('x1');
    ylabel('x3');

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
end

function out = sigmoid_io(x, v0, r)
    out = 0.5 * erf((x - v0) / (sqrt(2)*r)) + 0.5;
end

function out = gompertz_io(x, b, c, d)
    out = exp(-b*exp(-d*(x-c)));
end