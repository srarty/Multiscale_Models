%% Calculate the PSC and PSP of the three populations
% 
% Finding the parameters for the double exponential Wendling model
%
% Refs:
%   + Destexhe et al., 1996
%   + Wendling et al., 2002
%
% Artemio - November 2022
%% Params
% Synaptic time constants (seconds)
tauAMPA = 0.01172; %0.0024;
tauGABAa = 0.0011; %0.00525;
tauGABAb = 0.0229; %20*0.00525;

% Membrane time constants (seconds)
taui = 0.008457; % (interneurons membrane potential)
taup = 0.01788; % (Pyramidal membrane potential)

% Current gains (Amperes)
aAMPA = 300e-12;
aGABAa = 60e-12;
aGABAb = 60e-12;

alpha_ip = 10;
alpha_pp = -5.5;
alpha_pa = -1.19; % To ensure same AUC as the alpha function from wendling
alpha_pb = 90;
alpha_ab = 31.5;

%% Functions
% PSC
c = @(a,tau_syn,t) a*exp(-t/tau_syn); % Single exponential decay

% PSP
% TODO: Make the PSPs look like the wendling
v = @(a,tau_syn, tau_m, t) a*(exp(-t/tau_syn) - exp(-t/tau_m)); % Souble exponential function (incorporates the PSC)
alpha = @(a, tau, t) (a/tau).*t.*exp(-t/tau);

%% Run
T = 0:0.0001:0.5;

cAMPA = c(aAMPA, tauAMPA, T);
cGABAa = c(aGABAa, tauGABAa, T);
cGABAb = c(aGABAb, tauGABAb, T);

vip = v(alpha_ip, tauAMPA, taui, T);
vpp = v(alpha_pp, tauAMPA, taup, T);
vpa = v(alpha_pa, tauGABAa, taup, T);
vpb = v(alpha_pb, tauGABAb, taup, T);
vab = v(alpha_ab, tauGABAb, taui, T);

% Alpha functions from wendling 
% vp = alpha(3.25, 1/100, T);
% va = alpha(10, 1/500, T);
% vb = alpha(22, 1/50, T); 

%% Plot
% Currents
f = figure;
f.Position = [300 555 1252 420];
subplot(1,2,1)
plot([-0.05 0 T], [0 0 cAMPA*1e12], '-.');
hold;
plot([-0.05 0 T], [0 0 cGABAa*1e12], '-');
plot([-0.05 0 T], [0 0 cGABAb*1e12], '--');
legend({'AMPA', 'GABA_A', 'GABA_B'});
ylabel('PSC (pA)');
xlabel('Time (s)');

% Voltages
subplot(1,2,2)
plot(T, vip);
hold;
plot(T, vpa);
plot(T, vpb);
plot(T, vpp, '--');
plot(T, vab, '--');
legend({'AMPA', 'GABA_A', 'GABA_B'});
ylabel('PSP (mV)');
xlabel('Time (s)');