% NMM_GABAb.m It is the NMM with 3 populations (1 pyramidal and 2
% inhibitory: fast (GABA_A) and slow (GABA_B)).
%
function [x, y, t, f_e, f_i, f_b, params, yy] = NMM_GABAb(varargin)
    clear option
    if nargin >= 2
        option = cell(1,nargin/2);
        value = cell(1,nargin/2);
        for i = 2:2:nargin
            option{round(i/2)} = varargin{i-1};
            value{round(i/2)} = varargin{i};
        end
    end
    
    N = 2000; % Number of samples: 1 sample = 1 milisecond
    u = 0;

%     params = set_parameters('seizure', u);
    params = set_parameters('gabab', u);
%     params = set_parameters('wendling', u);
    params.time = N * params.dt;
    
    % Options  ------------------------------------------------------------
    params.options.ADD_NOISE = 1; % External input noise (0 = no noise, 1 = noise)
    params.options.CHANGE_U = 1; % 0: U doesn't change during simulation. Any other value of CHANGE_U: U changes.
    params.options.CHANGE_AGONIST = 0; % 1.5 X GABA_A gain at mid-simulation
    
    CURRENT = 0e-12; %50e-12;
    params.options.CURRENT_TIME = 1490:1500;
    params.options.INPUT_CURRENT_PY = 1000 * CURRENT / params.g_m_P; % 1000 for milivolts, then xe-12 A, where x is the amplitude in pA
    params.options.INPUT_CURRENT_IN = 1000 * CURRENT / params.g_m_I;    
    params.options.INPUT_CURRENT_B = 1000 * CURRENT / params.g_m_I;    
    % --------------------------------------------------------- End Options
    
    % Parse inputs --------------------------------------------------------
    if exist('option','var')
        for i = 1:numel(option)
            try
                params.(option{i}) = params.(option{i}) * value{i};
            catch
                error(['Couldn''t assign value: ' num2str(value{i}) ' to the parameter: ' option{i}]);
            end
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
    
    bb = params.gompertzb.b;
    bc = params.gompertzb.c;
    bd = params.gompertzb.d;
    
    e_0 = params.e0; 
    e_0i = params.e0i;
    e_0b = params.e0b;
    
    % Synaptic functions
    S1 = @(x) gompertz_io(x, e_0i, ib, ic ,id);% sigmoid_io(x, e_0i, v0, r); %gaussian_io(x, params.gaussiani.a, params.gaussiani.b, params.gaussiani.c, params.gaussiani.d);% 
    S2 = @(x) gompertz_io(x, e_0, b, c, d);  % sigmoid_io(x, e_0, v0, r); % gaussian_io(x, params.gaussian.a, params.gaussian.b, params.gaussian.c, params.gaussian.d);% 
    S3 = @(x) gompertz_io(x, e_0b, bb, bc ,bd);
        
    dt = params.dt;
    t = 0:dt:(N-1)*dt;

    alpha_i = params.alpha_i;
    alpha_e = params.alpha_e; 
    alpha_re = params.alpha_re; 
    alpha_ri = params.alpha_ri; 
    alpha_u = params.alpha_u; 
    alpha_b = params.alpha_b; 
    alpha_eb = params.alpha_eb; 
        
    eq = zeros(21,1);
%     eq = [-1.7257; -86.2826; 0.0008; 0.0766; 0.0009; 0.0441; -0.7756; -77.5645; 0; 0; -5.4558; -272.7924; 0.0008; 0.0766; -0.2635; 0.2813; 0.4010; -0.2501; 0.0615; -0.1143; 0.2813];
        
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
        eq(11);
        eq(12);
        eq(13);
        eq(14);
        alpha_i;
        alpha_e;
        alpha_re;
        alpha_ri;
        alpha_u;
        alpha_b;
        alpha_eb;
        ];
     
    
    [t,x,y] = ode45(@(t,x) ode(t,x,params,dt, S1, S2, S3), t, x0); % use "t", instead of "[min(t) max(t)]" fix the output time vector
    %     [t,x,y] = ode23(@(t,x) ode(t,x,params,dt), [min(t) max(t)], x0);
    %     [t,x,y] = ode113(@(t,x) ode(t,x,params,dt), [min(t) max(t)], x0);
    
    % Create injected current vector
    I_py = zeros(size(x,1),1);
    I_in = zeros(size(x,1),1);
    I_b = zeros(size(x,1),1);
    I_py(params.options.CURRENT_TIME) = params.options.INPUT_CURRENT_PY;
    I_in(params.options.CURRENT_TIME) = params.options.INPUT_CURRENT_IN;
    I_b(params.options.CURRENT_TIME) = params.options.INPUT_CURRENT_IN;
    
    yy = zeros(size(y));
    for i = 1:size(x,1)
        y(i) = x(i,1) + x(i,5) + x(i,9) + x(i,11) + I_py(i);
        yy(i) = (x(i,6) + x(i,10) - x(i,2) - x(i,12) + I_py(i))/params.g_m_P; % Current based LFP
    end
    
    % Calculate firing rate
    f_e = S2(x(:,1) + x(:,5) + x(:,9) + x(:,11) + I_py);
    f_i = S1(x(:,3) + x(:,7) + I_in); 
    f_b = S3(x(:,13) + I_b); 
    
    if nargin == 0        
        close all
        do_plot(x,t,y, yy, f_e, f_i, f_b);
    end
end

function dx = ode(t,x,params,dt, S1, S2, S3)
    
    u = params.u;
    
    % Following lines are meant to change the input mid simulation, comment
    % them to run it with constant input.
    if isfield(params,'options') & isfield(params.options,'CHANGE_U') & params.options.CHANGE_U
        if t >= 2*params.time/3
            u = 1;
        elseif t >= 1*params.time/3
            u = 0.5;
%         elseif t >= 1*params.time/4
%             u = 0.25;
        end
    elseif isfield(params,'options') & isfield(params.options,'CHANGE_AGONIST') & params.options.CHANGE_AGONIST
        if t >= params.time/2
            x(15) = 1.5 * params.alpha_i;
        end
    end
    
    if t >=  params.options.CURRENT_TIME(end) * 1e-3
        INPUT_CURRENT_PY = 0;
        INPUT_CURRENT_IN = 0;
        INPUT_CURRENT_B = 0;
    elseif t >= params.options.CURRENT_TIME(1) * 1e-3
        INPUT_CURRENT_PY = params.options.INPUT_CURRENT_PY;
        INPUT_CURRENT_IN = params.options.INPUT_CURRENT_IN;
        INPUT_CURRENT_B = params.options.INPUT_CURRENT_B;
    else
        INPUT_CURRENT_PY = 0;
        INPUT_CURRENT_IN = 0;
        INPUT_CURRENT_B = 0;
    end
        

    Tau_coeff = @(m, s) 1/(m*s);% Nicola Campbell
    
    c_constant = params.c_constant;
    c1 = 29.7217  * c_constant * params.P_inTOpy;   % Inhibitory synapses into pyramidal population
    c2 = 70.6482 * c_constant * params.P_pyTOin ;% Excitatory synapses into inhibitory population
    c3 = 141.1007 * c_constant * params.P_pyTOpy;   % Recursive excitation to pyramidal cells
    c4 = 8.9828   * c_constant * params.P_inTOin;   % Recursive inhibition to inhibitory cells
    c5 = 17.5913;% * c_constant;%* 2.5;                                   % External excitatory synapses into pyramidal population
    c6 = 151.7094 * c_constant * params.P_inTOpy;      % Inhibitory synapses into pyramidal population % 166.4415
    c7 = 70.6482 * c_constant * params.P_pyTOin;    % Excitatory to GABAb
    
    tau_sp = params.tau_sp;
    tau_mp = params.tau_mp;
    
    tau_si = params.tau_si;
    tau_mi = params.tau_mi;
    
    tau_srp = params.tau_srp;
    tau_mrp = params.tau_mrp;
    
    tau_sri = params.tau_sri;
    tau_mri = params.tau_mri;
    
    tau_sb = params.tau_sb;
    
    % Lumped gain parameters
    AmplitudeI  = c1 * x(15) * Tau_coeff(tau_mp,  tau_sp);
    AmplitudeE  = c2 * x(16) * Tau_coeff(tau_mi,  tau_si);
    AmplitudeRE = c3 * x(17) * Tau_coeff(tau_mrp, tau_srp);
    AmplitudeRI = c4 * x(18) * Tau_coeff(tau_mri, tau_sri);
    AmplitudeU  = c5 * x(19) * Tau_coeff(tau_mrp, tau_srp);
    AmplitudeB  = c6 * x(20) * Tau_coeff(tau_mp,  tau_sb);
    AmplitudeEB = c7 * x(21) * Tau_coeff(tau_mi,  tau_si);
    
    % Diff equations ------------------------------------------------------
    dx = zeros(21,1);

    % Double exponential from Nicola-Campbell (2013):
    % GABAa -> P
    dx(1) = x(2) - x(1)/tau_mp;
    dx(2) = - x(2)/tau_sp + AmplitudeI * S1(x(3) + x(7) + INPUT_CURRENT_IN);
    % P -> GABAa
    dx(3) = x(4) - x(3)/tau_mi;
    dx(4) = - x(4)/tau_si + AmplitudeE * S2(x(1) + x(5) + x(9) + x(11) + INPUT_CURRENT_PY);
    % Recurrent Pyramidal P -> P
    dx(5) = x(6) - x(5)/tau_mrp;
    dx(6) = - x(6)/tau_srp + AmplitudeRE * S2(x(1) + x(5) + x(9)+ x(11) + INPUT_CURRENT_PY );
    % Recurrent GABAb -> GABAa
    dx(7) = x(8) - x(7)/tau_mri;
	dx(8) = - x(8)/tau_sri + AmplitudeRI * S3(x(13) + INPUT_CURRENT_IN);
    % External input u->P
    dx(9) = x(10) - x(9)/tau_mrp;
    dx(10) = - x(10)/tau_srp + AmplitudeU * (u + (params.options.ADD_NOISE *(sqrt(u).*randn(1,1))));
    % GABAb -> P
    dx(11) = x(12) - x(11)/tau_mp;
    dx(12) = - x(12)/tau_sb + AmplitudeB * S3(x(13) + INPUT_CURRENT_B);
    % P -> GABAb
    dx(13) = x(14) - x(13)/tau_mi;
    dx(14) = - x(14)/tau_si + AmplitudeEB * S2(x(1) + x(5) + x(9) + x(11) + INPUT_CURRENT_PY);
    % Parameters:
    dx(15) = 0; % alpha_i
    dx(16) = 0; % alpha_e
    dx(17) = 0; % alpha_re
    dx(18) = 0; % alpha_ri
    dx(19) = 0; % alpha_u
    dx(20) = 0; % alpha_b
    dx(21) = 0; % alpha_eb

end

function do_plot(x,t, Vm, LFP_I, f_e, f_i, f_b)    
    figure
    plot(t,x(:,[1 3 5 7 9 11 13]));
    legend({'x1' 'x3' 'x5' 'x7' 'x9' 'x11' 'x13'});
    ylabel('mV');
    xlabel('Time (s)');
    
    figure
    plot(t, x(:,1) + x(:,5) + x(:,9) + x(:,11));
    hold
    plot(t, x(:,3) + x(:,7));
    plot(t, x(:,13))
    legend({'V_{mp}' 'V_{mi}' 'V_{mGABA_{b}}'});
    ylabel('mV');
    xlabel('Time (s)');

    figure
    ax1 = subplot(2,1,1);
    plot(t, Vm, 'k');
    title('Vm_{Py}');
    ylabel('mV');
    xlabel('Time (s)');
    ax2 = subplot(2,1,2);
    plot(t, LFP_I, 'k');
    title('LFP (from current)');
    ylabel('pA?');
    xlabel('Time (s)');
    linkaxes([ax1 ax2],'x')
    
    figure
    plot(t, f_e); hold on;
    plot(t, f_i);
    plot(t, f_b);
    ylabel('Spike rate (Hz)');
    xlabel('Time (s)');
    legend({'Pyramidal', 'GABA_A', 'GABA_B'}, 'Location', 'best');
    
%     figure
%     plot(t, x(:,9));    
%     
%     figure
%     subplot(1,2,1)
%     plot(x(:,1), x(:,3));
%     subplot(1,2,2)
%     plot(x(:,1), [0; diff(x(:,1))]);
end
