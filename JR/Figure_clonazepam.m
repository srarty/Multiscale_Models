% Figures to show the effect of GABA_A agonists
%   + barbiturates (e.g. ??) - increasing the duration of chloride ion
%       channel opening at the GABA_A receptor, i.e.increases the time 
%       constant of GABA_A inhibitory synapses.
%       NOTE: increases the toxicity in overdose because it directly
%       controls the opening of GABA_A Cloride ion channel.
%   + benzodiazepines (e.g. midazolam, diazepam) -  increase the frequency 
%       of the chloride ion channel opening at the GABAA receptor, i.e. 
%       increases the potency of GABA). 
%       NOTE: indirect action on the cloride channel.
%       REFS: Rovira et al. 1993 (Fig 2); 
%
% Beta activity increase: Nunez, page 435. , Niedermeyers and Lopes da
% Silva 1999, 
%
% TODO: Check which concentration or dose corresponds to which increase in
%       GABA_A gain (Rovira 1993)
% TODO: Bar plot of Beta-band energy with different concentrations of the
%       agonist

%% Options
% normalize = @(x) x./rms(x);
% normalize = @(x) x./max(x) - min(x);
% normalize = @(x) x;
normalize = @(x) x - mean(x);


%% NMM - Simulate
% Simulates NMM, one with default parameters, one with augmented GABA_A like diazepam
[~, y_d]=NMM_GABAb('alpha_i', 1);
% [~, y_a, t]=NMM_GABAb('tau_sp',2, 'tau_sri', 2);% Barbiturate (GABAa agonist, incereases Cl- channel opening)
% [~, y_a, t]=NMM_GABAb('alpha_i', 3); % Benzodiazepine (GABAa agonist, increases gain)
% [~, y_a, t]=NMM_GABAb('alpha_b', 3); % Baclofen (Increases inhibition in GABAb receptors, decreases release of excitatory neurotransimtters in presynaptic neurons)

range = 1000:length(t);
Fs = 1/diff(t(1:2));

t = t(range);
y_default = y_d(range);
y_agonist = y_a(range);

%% Beta band
b_default =bandpass(normalize(y_default)', [13 20], Fs);
b_agonist =bandpass(normalize(y_agonist)', [13 20], Fs);

E_default = sum(b_default.^2)/sum(y_default.^2);
E_agonist = sum(b_agonist.^2)/sum(y_agonist.^2);

%% Plot NMM
figure;
% Raw
subplot(2,1,1)
plot(t, normalize(y_default)');
hold
plot(t, normalize(y_agonist)');
ylabel('LFP');
box off
title('NMM')
l = legend({'Default', 'GABA_A agonist'});
l.Location = 'best';

% Beta
subplot(2,1,2)
plot(t, b_default);
hold
plot(t, b_agonist);
ylabel('\beta{}-band');
xlabel('Time (s)');
box off

% Spectrum
[fig, Xs_default] = fft_plot(normalize(y_default)',t');
[~, Xs_agonist, Freqs] = fft_plot(normalize(y_agonist)',t', fig);
xline(13, '--', 'LineWidth', 1.5, 'Color', [0.5 0.5 0.5]);
xline(20, '--', 'LineWidth', 1.5, 'Color', [0.5 0.5 0.5]);
box off
l = legend({'Default', 'GABA_A agonist', '\beta band'});
l.Location = 'best';
xlim([0 60]);
title('NMM');
ylabel('|X|');
xlabel('Frequency (Hz)');


return

%% LIF
file_default = 'C:/Users/artemios/Documents/Multiscale_Models_Data/lfp_77.mat'; % 1 x j_pi
% file_default = 'C:/Users/artemios/Documents/Multiscale_Models_Data/lfp_82.mat'; % 1 x j_pi (non-random th and refractory period)

% file_agonist = 'C:/Users/artemios/Documents/Multiscale_Models_Data/lfp_78.mat'; % 1.5 x j_pi
% file_agonist = 'C:/Users/artemios/Documents/Multiscale_Models_Data/lfp_81.mat'; % 1.5 x j_pi (non-random th and refractory period)
% file_agonist = 'C:/Users/artemios/Documents/Multiscale_Models_Data/lfp_83.mat'; % 2 x j_pi and j_bi (non-random th and refractory period)
file_agonist = 'C:/Users/artemios/Documents/Multiscale_Models_Data/lfp_84.mat'; % 2 x j_pi and j_bi (random th and refractory period)
% file_agonist = 'C:/Users/artemios/Documents/Multiscale_Models_Data/lfp_85.mat'; % 1.5 x j_pi and j_bi (random th and refractory period)


% Load
agonist = load(file_agonist);
default = load(file_default);

t2 = agonist.lfp_dt:agonist.lfp_dt:agonist.lfp_dt * length(agonist.LFP_V);
range2 = 6000:length(t2);
Fs2 = 1/diff(t2(1:2));

lfp_agonist = agonist.LFP_V(range2);
lfp_default = default.LFP_V(range2);

%% Beta band
beta_default = bandpass(normalize(lfp_default)', [13 20], Fs2);
beta_agonist = bandpass(normalize(lfp_agonist)', [13 20], Fs2);

Elif_default = sum(beta_default.^2)/sum(y_default.^2);
Elif_agonist = sum(beta_agonist.^2)/sum(y_agonist.^2);

%% Plot LIF
% Raw
figure;
subplot(2,1,1)
plot(t2(range2), normalize(lfp_default));
hold
plot(t2(range2), normalize(lfp_agonist));
ylabel('LFP');
box off
title('LIF')
l = legend({'Default', 'GABA_A agonist'});
l.Location = 'best';

% Beta band
subplot(2,1,2)
plot(t2(range2), beta_default);
hold
plot(t2(range2), beta_agonist);
box off
ylabel('\beta{}-band');
xlabel('Time (s)');

% Spectrum
[fig2, Xx] = fft_plot(normalize(lfp_default), t2(range2));
[~, Xx_, F_] = fft_plot(normalize(lfp_agonist), t2(range2), fig2);
xline(13, '--', 'LineWidth', 1.5, 'Color', [0.5 0.5 0.5]);
xline(20, '--', 'LineWidth', 1.5, 'Color', [0.5 0.5 0.5]);
box off
l = legend({'Default', 'GABA_A agonist', '\beta band'});
l.Location = 'best';
xlim([0 60]);
title('LIF');
ylabel('|X|');
xlabel('Frequency (Hz)');

%}