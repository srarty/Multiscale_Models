%% Nullclines
% External input
u_bkg = 1;
u_ext = 0;

u = u_bkg + u_ext;

% Nonlinearity parameters (constant)
M = 41.18;
a = 2.919;
s = 24.8;
Mi = 69.29;
ai = 2.734;
si = 23.94;

% Connectivity constants
c1 = 12.78;
c2 = 30.21;
c3 = 34.02;
c4 = 4.419;
c5 = 17.6;
c6 = 28.13;
c7 = 8.8;

% Gain parameters
alpha_i  = -0.705;
alpha_e  = 2.472;
alpha_re = 1.255;
alpha_ri = -3.316;
alpha_u  = 1.637;
alpha_b  = -0.32;
alpha_ui = 2.591;

tau_mp=0.02;
tau_mi=0.01;
tau_sp=0.00525;
tau_si=0.0012;
tau_srp = 0.0024;
tau_sri = 0.00525;
tau_sb = 0.012;

% Naka-Rushton firing-rate function
S1 = @(x) heaviside(x) .* (Mi .* x.^ai ./ (si.^ai + x.^ai));
S2 = @(x) heaviside(x) .* (M .* x.^a ./ (s.^a + x.^a));

% Tau coefficient
Tau_coeff = @(m, s) 1/(m*s);

% Calculate gain amplitude values
AmplitudeI  = c1 * alpha_i  * Tau_coeff(tau_mp,  tau_sp);
AmplitudeE  = c2 * alpha_e  * Tau_coeff(tau_mi,  tau_si);
AmplitudeRE = c3 * alpha_re * Tau_coeff(tau_mp, tau_srp);
AmplitudeRI = c4 * alpha_ri * Tau_coeff(tau_mi, tau_sri);
AmplitudeU  = c5 * alpha_u  * Tau_coeff(tau_mp, tau_srp);
AmplitudeB  = c6 * alpha_b  * Tau_coeff(tau_mp,  tau_sb);
AmplitudeUi = c7 * alpha_ui * Tau_coeff(tau_mi, tau_si);

% Equilibrium point:
x1 = -15.055811;
x2 = -752.79225;
x3 = 7.8693759;
x4 = 787.53048;
x5 = 4.4994306;
x6 = 224.97413;
x7 = 28.8112;
x8 = 1440.56;
x9 = -15.041952;
x10 = -752.09746;
x11 = -24.486244;
x12 = -2448.633;
x13 = 22.801494;
x14 = 2279.5709;

%% Define nullclines
nullclineX1 = @(x2) -x2 * tau_mp; 
nullclineX2 = @(x3, x11, x13) tau_sp * AmplitudeI * S1(x3 + x11 + x13);
nullclineX3 = @(x4) -x4 * tau_mi; 

%% Plot
%% X1 vs X2
x1 = [-100:0.1:100];
x2 = [-1000:0.1:100];

figure
plot(x1 , nullclineX2(x3, x11, x13)*ones(size(x1)) );hold
plot(nullclineX1(x2) , x2);
xlabel('x1');
ylabel('x2');

%% X1 vs X3
x1 = [-100:0.1:100];
x3 = [-100:0.1:100];
x2 = -752.79225;

figure
plot(x1 , nullclineX3(x4)*ones(size(x1)) );hold
plot(nullclineX1(x2)*ones(size(x3)) , x3);
xlabel('x1');
ylabel('x3');
