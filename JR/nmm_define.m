% NMM_DEFINE Neural Mass Model based on Jensen and Rit. Defines the state based representation of the model.
%
% Inputs:
%   x0 - initial states
%   P0 - initial covariance matrix
%   params - parameters defined by Jensen and Rit
%   C_CONSTANT - (optional) Default is 135
%
% Outputs:
%   nmm - (struct) the neural mass model
%
function nmm = nmm_define(x0,P0,params, varargin)
if nargin > 3, C_CONSTANT = varargin{1}; else, C_CONSTANT = 135; end

% Indexes
v_idx = [1 3];
z_idx = [2 4];
u_idx = 5;
alpha_idx = [6 7];

% the parameters
dt          = params.dt;
e_0         = params.e0;
r           = params.r;	% varsigma
v0          = params.v0; % Threshold
decay_e     = params.decay_e; % inverse time constants 
decay_i     = params.decay_i; % 
alpha_ei    = params.alpha_ei; % synaptic gains (excitatory)
alpha_ie    = params.alpha_ie; % (inhibitory)
u           = params.u;	% mean input firing rate.
scale       = params.scale; % Scale to fix mismatch in state amplitudes. Not to be confused with the scale in analytic_kalman_filter_2

c_constant = C_CONSTANT; %135;%675; % This should reflect Brunel's 'N' variable
% c1 = 0.25*c_constant;%1*c_constant;	% number of synapses
% c2 = 0.25*c_constant;%0.8*c_constant;
% c1 = 4*c_constant;	% number of synapses (matches brunel?)
% c2 = 1*c_constant; % (matches brunel?)
c1 = 1*c_constant;
c2 = 0.8*c_constant;
% c1 = 4 * c_constant * params.P_pyTOin; % Excitatory synapses into inhibitory population
% c2 = 1 * c_constant * params.P_inTOpy; % Inhibitory synapses into pyramidal population

% Number of augmented states
xlen = length(x0);

% Linear component of model
A =     [1,                  dt*scale,      0,              0,          0,  0,  0; ...
  -decay_i^2*dt/scale,	1-2*decay_i*dt,     0,              0,          0,  0,  0; ...
         0,                     0,          1,            dt*scale,     0,  0,  0; ...
         0,                     0  -decay_e^2*dt/scale,	1-2*decay_e*dt,	0,  0,  0; ...
         0,                     0,          0,              0,          1,  0,  0; ...
         0,                     0,          0,              0,          0,  1,  0; ...
         0,                     0,          0,              0,          0,  0,  1];
     
% B Matrix (Augmented parameters)                                                                
%      
% B =     [0,	0,	0,	0,	0,	0,	0; ...
%          0,	0,	0,	0,	0,	1,	0; ...
%          0,	0,	0,	0,	0,	0,	0; ...
%          0,	0,	0,	0,	0,	0,	1; ...
%          0,	0,	0,	0,	0,	0,	0; ...
%          0,	0,	0,	0,	0,	0,	0; ...
%          0,	0,	0,	0,	0,  0,  0];
%
B = zeros(xlen);
B(z_idx, alpha_idx) = diag(ones(size(z_idx)));
% B(4,7) = 1; % ?
% B(2,6) = 1; % ?

% C Matrix (Augmented)
% C =     [0,	0,	0,	0,	0,	0,	0; ... 
%          0,	0,	1,	0,	0,	0,	0; ... % phi(x3)
%          0,	0,	0,	0,	0,	0,	0; ...
%          1,	0,	0,	0,	1,	0,	0; ... % phi(x5+x1) % x1 is negative thanks to alpha_i
%          0,	0,	0,	0,	0,	0,	0; ...
%          0,	0,	0,	0,	0,	0,	0; ...
%          0,	0,	0,	0,	0,  0,  0];
%
C = zeros(xlen);
% Basic connections
C(2,3) = 1; % phi(x3)
C(4,1) = 1; % x1
C(4,u_idx) = 1; % x5
% Extras
% C(2,u_idx) = 1; % input -> inhibitory % Only uncomment this to test external excitatory inputs to the inhibitory population
% C(4,3) = 1; % inhibitory -> inhibitory (recurrent inhibition) Need to add the synapses
% C(2,1) = 1; % exc -> exc
C = C./scale;

alpha_i = alpha_ie * c2 * 2 * e_0 * dt * decay_i; % lumped constant (inhibitory input to pyramidal cells)
alpha_e = alpha_ei * c1 * 2 * e_0 * dt * decay_e; % lumped constant (excitatory input to inhibitory interneurons)

delump(1) = 1/(c2 * 2 * e_0 * dt * decay_i);
delump(2) = 1/(c1 * 2 * e_0 * dt * decay_e);

% SCALE 1 - this is to avoid large differences between states upsetting the filter 
% (magnitude of the membrane potentials and their derivatives)
input = scale*u(1);
% SCALE 2 - this converts a constant input to its effect on the pyramidal
% membrane potential by taking the steady state limit of the synaptic kernel
% (assumption that the input varies much slower than the state variables).
% input = input * alpha_ei*decay_e /decay_e^2;
%       ~~~~~   ~~~~~~~~~~~~~~   ~~~~~~~~~
%       input   synaptic gain    integral of kernel

x0(5) = input(1);
x0(6) = alpha_i; % Inhibitory synapses, input to pyramidal cells
x0(7) = alpha_e; % Excitatory synapses, input to inhibitory interneurons

nmm = struct;
nmm.A = A;
nmm.B = B;
nmm.C = C;
nmm.x0 = x0;
nmm.P0 = P0;
nmm.delump = delump;
nmm.params = params;
nmm.options = struct;
nmm.iteration = 1;

end % end function - nmm_define
