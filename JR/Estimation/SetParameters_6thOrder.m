

% set parameters for the 6th order JR model

mu = 220;%370;%520;%220        % this is just the mean when this function is used in estimation
sigma = 5.74;


% get current parameters

e_0 = 2.5;        % maximal firing parameter
v_0 = 6;        % firing threshold
if (log_sig == 1)
    % logistic sigmoid
    r = 0.56;              % sigmoid slope
else
    % erf sigmoid
    r = 3.0285;
end


A = 3*3.25;%3.25;          % excitatory gain
a = 100;%;250;%100;%105;          % excitatory time constant
B = 3*22; %22;        % inhibitory gain
b = 50;%50%50;%85;         % inhibitory time constant

% connectivity parameters from J and R model
C = 135;%135;        % 135 for alpha activity
C1 = C;        % connectivity parameter - excitatory feedback PSP
C2 = 0.8*C;        % connectivity parameter - excitatory feedback firing
C3 = 0.25*C;        % connectivity parameter - inhibitory feedback PSP
C4 = 0.25*C;        % connectivity parameter - excitatory feedback firing

nrM = 5;%38.48; Asymmetric
nra = 3.528;
nrb = 0-10;
nrs = 23.53;
% nrM = 5.008; Symmetric
% nra = 9.007;
% nrb = -10;
% nrs = 15.88;

ga = 5.01;
gb = 2.439;
gc = 2.547;
gd = 0.3911;