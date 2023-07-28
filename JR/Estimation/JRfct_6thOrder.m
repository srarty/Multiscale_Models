function [out, f_p, f_e, f_i] = JRfct_6thOrder(parameters,x)

global dq log_sig

p = x(1:dq,:);                      % these are the sigma points for the parameters
xnl = x(dq+1:size(x,1),:);      % these are the sigma points for the states

dt = parameters{1};         % time step
mu = parameters{2};        % this is just the mean when this function is used in estimation

% get current parameters
e_0 = parameters{3};        % maximal firing parameter

if dq== 3
    v_0 = parameters{4};        % firing threshold
    r = parameters{5};              % sigmoid slope

    A = parameters{6};% p(2,:); %parameters{6};          % excitatory gain
    a = parameters{7};          % excitatory time constant
    B = p(1,:);%parameters{8};         % inhibitory gain
    b = p(2,:);%parameters{9};         % inhibitory time constant
    mu = p(3,:);                % mean input
    
    nrM = parameters{14};
    nra = parameters{15};
    nrb = parameters{16};
    nrs = parameters{17}; 
    
    ga = parameters{18}; 
    gb = parameters{19}; 
    gc = parameters{20}; 
    gd = parameters{21}; 
elseif dq==5
    v_0 = parameters{4};        % firing threshold
    r = parameters{5};              % sigmoid slope

    A = p(4,:); %parameters{6};          % excitatory gain
    a = p(5,:);          % excitatory time constant
    B = p(1,:);%parameters{8};         % inhibitory gain
    b = p(2,:);%parameters{9};         % inhibitory time constant
    mu = p(3,:);                % mean input
    
    nrM = parameters{14};
    nra = parameters{15};
    nrb = parameters{16};
    nrs = parameters{17}; 
    
    ga = parameters{18}; 
    gb = parameters{19}; 
    gc = parameters{20}; 
    gd = parameters{21}; 
end

C1 = parameters{10};        % connectivity parameter - excitatory feedback PSP
C2 = parameters{11};        % connectivity parameter - excitatory feedback firing
C3 = parameters{12};        % connectivity parameter - inhibitory feedback PSP
C4 = parameters{13};        % connectivity parameter - excitatory feedback firing

% get current states
v_f = xnl(1,:);
z_f = xnl(2,:);
v_p1 = xnl(3,:);
z_p1 = xnl(4,:);
v_p2 = xnl(5,:);
z_p2 = xnl(6,:);

v_p = v_p1 - v_p2;                   % the minus is where the inhibition is introduced

if (log_sig == 1)
    f_p = 2*e_0 ./ (1 + exp(r.*(v_0 - v_p)));           % pyramidal firing
    f_e = 2*e_0 ./ (1 + exp(r.*(v_0 - C1.*v_f)));       % excitatory feedback firing
    f_i = 2*e_0 ./ (1 + exp(r.*(v_0 - C3.*v_f)));       % inhibitory feedback firing
elseif (log_sig == 0)
    f_p = e_0.*(erf((v_p-v_0)./(sqrt(2).*r)) + 1);           % pyramidal firing
    f_e = e_0.*(erf((C1.*v_f-v_0)./(sqrt(2).*r)) + 1);       % excitatory feedback firing
    f_i = e_0.*(erf((C3.*v_f-v_0)./(sqrt(2).*r)) + 1);       % inhibitory feedback firing   
elseif log_sig == 2
    % Running the transfer function:
    f_p = naka_rushton_io(v_p, nrM, nra, nrb, nrs); % pyramidal firing
    f_e = naka_rushton_io(C1.*v_f, nrM, nra, nrb, nrs); % excitatory feedback firing
    f_i = naka_rushton_io(C3.*v_f, nrM, nra, nrb, nrs); % inhibitory feedback firing   
elseif log_sig == 3
    % Parameters (from the fit)
%     ga = 5.01;
%     gb = 2.439;
%     gc = 2.547;
%     gd = 0.3911;
    % Running the transfer function:
    f_p = gompertz_io(v_p, ga, gb, gc, gd); % pyramidal firing
    f_e = gompertz_io(C1.*v_f, ga, gb, gc, gd); % excitatory feedback firing
    f_i = gompertz_io(C3.*v_f, ga, gb, gc, gd); % inhibitory feedback firing   
end

% state update equation - propagate states
xnl(1,:) = v_f + dt*z_f;
xnl(2,:) = z_f + dt*(A.*a.*f_p - 2.*a.*z_f - a.^2.*v_f);
xnl(3,:) = v_p1 + dt*z_p1;
xnl(4,:) = z_p1 + dt*(A.*a.*(mu + C2.*f_e) - 2.*a.*z_p1 - a.^2.*v_p1);
xnl(5,:) = v_p2 + dt*z_p2;
xnl(6,:) = z_p2 + dt*(B.*b.*C4.*f_i - 2.*b.*z_p2 - b.^2.*v_p2);

% augmented state output
out = [x(1:dq,:)
    xnl];