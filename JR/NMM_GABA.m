% NMM_GABAb.m It is the NMM with 3 populations (1 pyramidal and 2
% inhibitory: fast (GABA_A) and slow (GABA_B)).
%
function [x, y, t, f_e, f_i, params, yy] = NMM_GABA(varargin)
    clear parameter
    if nargin >= 2
        parameter = cell(1,nargin/2);
        value = cell(1,nargin/2);
        for i = 2:2:nargin
            % In case the input is a regular parameter:
            parameter{round(i/2)} = varargin{i-1};
            value{round(i/2)} = varargin{i};
        end
    end
    
    N = 2000; % Number of samples: 1 sample = 1 milisecond
    u = 1;

%     params = set_parameters('seizure', u);
    params = set_parameters('gabab', u);
    params.time = N * params.dt;
    
    % Options  ------------------------------------------------------------
    params.options.ADD_NOISE = 1; % External input noise (0 = no noise, 1 = noise)
    params.options.CHANGE_U = 0; % 0: U doesn't change during simulation. Any other value of CHANGE_U: U changes.
    params.options.CHANGE_AGONIST = 0; % Agonist changes
    
    CURRENT = 0;% 50e-12; % 50
    if exist('injected_current','var'), CURRENT = injected_current; end % If 'CURRENT' was a varargin, ignore previous line
    params.options.CURRENT_TIME = 1490:1500; %1490:1500;
    params.options.INPUT_CURRENT_PY = 1000 * CURRENT / params.g_m_P; % 1000 for milivolts, then xe-12 A, where x is the amplitude in pA
    params.options.INPUT_CURRENT_IN = 1000 * CURRENT / params.g_m_I;   
    % --------------------------------------------------------- End Options
    
    % Parse inputs --------------------------------------------------------
    if exist('parameter','var')
        for i = 1:numel(parameter)
            try
                params.(parameter{i}) = params.(parameter{i}) * value{i};
            catch E
                if strcmp('MATLAB:nonExistentField', E.identifier) && strcmp('CURRENT', parameter{i})
                    % The parameter to modify is the injected current
                    params.options.INPUT_CURRENT_PY = 1000 * value{i} / params.g_m_P;
                    params.options.INPUT_CURRENT_IN = 1000 * value{i} / params.g_m_I;
                    disp('The injected current has been modified by an input argument to the function NMM_GABA.m')
                else
                    error(['Couldn''t assign value: ' num2str(value{i}) ' to the parameter: ' parameter{i}]);
                end
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
    
    e_0 = params.e0; 
    e_0i = params.e0i;
    
    % Synaptic functions
    S1 = @(x) gompertz_io(x, e_0i, ib, ic ,id);% sigmoid_io(x, e_0i, v0, r); %gaussian_io(x, params.gaussiani.a, params.gaussiani.b, params.gaussiani.c, params.gaussiani.d);% 
    S2 = @(x) gompertz_io(x, e_0, b, c, d);  % sigmoid_io(x, e_0, v0, r); % gaussian_io(x, params.gaussian.a, params.gaussian.b, params.gaussian.c, params.gaussian.d);% 
        
    dt = params.dt;
    t = 0:dt:(N-1)*dt;

    alpha_i = params.alpha_i;
    alpha_e = params.alpha_e; 
    alpha_re = params.alpha_re; 
    alpha_ri = params.alpha_ri; 
    alpha_u = params.alpha_u; 
    alpha_b = params.alpha_b; 
        
    eq =zeros(21,1);
        
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
        alpha_i;
        alpha_e;
        alpha_re;
        alpha_ri;
        alpha_u;
        alpha_b;
        0;
        0;
        ];
     
    
    [t,x,y] = ode45(@(t,x) ode(t,x,params,dt, S1, S2), t, x0); % use "t", instead of "[min(t) max(t)]" fix the output time vector
    %     [t,x,y] = ode23(@(t,x) ode(t,x,params,dt), [min(t) max(t)], x0);
    %     [t,x,y] = ode113(@(t,x) ode(t,x,params,dt), [min(t) max(t)], x0);
    
    % Create injected current vector
    I_py = zeros(size(x,1),1);
    I_in = zeros(size(x,1),1);
    I_py(params.options.CURRENT_TIME) = params.options.INPUT_CURRENT_PY;
    I_in(params.options.CURRENT_TIME) = params.options.INPUT_CURRENT_IN;
    
    yy = zeros(size(y));
    for i = 1:size(x,1)
        y(i) = x(i,1) + x(i,5) + x(i,7) + x(i,9) + I_py(i);
        yy(i) = (x(i,6) + x(i,10) - x(i,2) - x(i,8) + I_py(i))/params.g_m_P; % Current based LFP
    end
    
    % Calculate firing rate
    f_e = S2(x(:,1) + x(:,5) + x(:,7) + x(:,9) + I_py);
    f_i = S1(x(:,3) + x(:,11) + x(:,19) + I_in); 
    
    if nargin == 0        
        close all
        do_plot(x,t,y, yy, f_e, f_i);
    end
end

function dx = ode(t,x,params,dt, S1, S2)
    
    u = params.u;
    
    % Following lines are meant to change the input mid simulation, comment
    % them to run it with constant input.
    if isfield(params,'options') & isfield(params.options,'CHANGE_U') & params.options.CHANGE_U
        if t >= 3*params.time/4
            u = 2;
        elseif t >= 2*params.time/4
            u = 1;
        elseif t >= 1*params.time/4
            u = 0.5;
        end
        
    elseif isfield(params,'options') & isfield(params.options,'CHANGE_AGONIST') & params.options.CHANGE_AGONIST
        if t >= 1
            x(13) = 1.5 * params.alpha_i;
            x(16) = 1.5 * params.alpha_ri;
        end
    end
    
    if t >=  params.options.CURRENT_TIME(end) * 1e-3
        INPUT_CURRENT_PY = 0;
        INPUT_CURRENT_IN = 0;
    elseif t >= params.options.CURRENT_TIME(1) * 1e-3
        INPUT_CURRENT_PY = params.options.INPUT_CURRENT_PY;
        INPUT_CURRENT_IN = params.options.INPUT_CURRENT_IN;
    else
        INPUT_CURRENT_PY = 0;
        INPUT_CURRENT_IN = 0;
    end
        

    Tau_coeff = @(m, s) 1/(m*s);% Nicola Campbell
    
    c_constant = params.c_constant;
    c1 = 29.7217  * c_constant * params.P_inTOpy;   % Inhibitory synapses into pyramidal population
    c2 = 70.6482 * c_constant * params.P_pyTOin ;% Excitatory synapses into inhibitory population
    c3 = 141.1007 * c_constant * params.P_pyTOpy;   % Recursive excitation to pyramidal cells
    c4 = 9.5   * c_constant * params.P_inTOin;   % Recursive inhibition to inhibitory cells
    c5 = 17.6;                                      % External excitatory synapses into pyramidal population
    c6 = 170.3 * c_constant * params.P_inTOpy;      % Inhibitory synapses into pyramidal population % 166.4415
    
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
    AmplitudeI  = c1 * x(13) * Tau_coeff(tau_mp,  tau_sp);
    AmplitudeE  = c2 * x(14) * Tau_coeff(tau_mi,  tau_si);
    AmplitudeRE = c3 * x(15) * Tau_coeff(tau_mrp, tau_srp);
    AmplitudeRI = c4 * x(16) * Tau_coeff(tau_mri, tau_sri); % This is the same for GABAa->GABAa and GABAa->GABAb
    AmplitudeU  = c5 * x(17) * Tau_coeff(tau_mrp, tau_srp);
    AmplitudeU_interneuorns  = 0.3 * AmplitudeU;
    AmplitudeB  = c6 * x(18) * Tau_coeff(tau_mp,  tau_sb);
    
    %% Diff equations ------------------------------------------------------
    dx = zeros(18,1);

    % Double exponential from Nicola-Campbell (2013):
    % GABAa->P
    dx(1) = x(2) - x(1)/tau_mp;
    dx(2) = - x(2)/tau_sp + AmplitudeI *  S1(x(3) + x(11) + x(19) + INPUT_CURRENT_IN);
    
    % P->I
    dx(3) = x(4) - x(3)/tau_mi;
    dx(4) = - x(4)/tau_si + AmplitudeE *  S2(x(1) + x(5) + x(7) + x(9) + INPUT_CURRENT_PY);
    
    % Recurrent Pyramidal P->P
    dx(5) = x(6) - x(5)/tau_mrp;
    dx(6) = - x(6)/tau_srp + AmplitudeRE *  S2(x(1) + x(5) + x(7)+ x(9) + INPUT_CURRENT_PY );
    
    % External input u->P
    dx(7) = x(8) - x(7)/tau_mrp;
    dx(8) = - x(8)/tau_srp + AmplitudeU  *  (u + (params.options.ADD_NOISE *(sqrt(u).*randn(1,1))));
    
    % GABAb -> P
    dx(9) = x(10) - x(9)/tau_mp;
    dx(10) = - x(10)/tau_sb + AmplitudeB *  S1(x(3) + x(11) + x(19) + INPUT_CURRENT_IN); % S3(x(13) + x(7) + INPUT_CURRENT_B);
    
    % Recurrent Inhibition GABAa -> I
    dx(11) = x(12) - x(11)/tau_mri;
	dx(12) = - x(12)/tau_sri + AmplitudeRI * S1(x(3) + x(11) + x(19) + INPUT_CURRENT_IN);
    
    % External input u->I
    dx(19) = x(20) - x(19)/tau_mi;
    dx(20) = - x(20)/tau_si + AmplitudeU_interneuorns  * (u + (params.options.ADD_NOISE *(sqrt(u).*randn(1,1))));
    
    % Parameters:
    dx(13) = 0; % alpha_i
    dx(14) = 0; % alpha_e
    dx(15) = 0; % alpha_re
    dx(16) = 0; % alpha_ri
    dx(17) = 0; % alpha_u
    dx(18) = 0; % alpha_b

end

function do_plot(x,t, Vm, LFP_I, f_e, f_i)    
    figure
    plot(t,x(:,[1 3 5 7 9 11 19]));
    legend({'x1' 'x3' 'x5' 'x7' 'x9' 'x11'});
    ylabel('mV');
    xlabel('Time (s)');
    
    figure
    plot(t, x(:,1) + x(:,5) + x(:,7) + x(:,9));
    hold
    plot(t, x(:,3) + x(:,11) +  + x(:,19));
    legend({'V_{mp}' 'V_{mi}'});
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
    ylabel('Spike rate (Hz)');
    xlabel('Time (s)');
    legend({'Pyramidal', 'Interneurons'}, 'Location', 'best');
    
%     figure
%     plot(t, x(:,7));    
%     
%     figure
%     subplot(1,2,1)
%     plot(x(:,1), x(:,3));
%     subplot(1,2,2)
%     plot(x(:,1), [0; diff(x(:,1))]);
end
