% Plots NMM's LFP and FFT for diff combinations of ri and i gains

cmap =[0 0 0; 0 0 1; 1 0 1; 1 0.8 0; 1 0 0];
f_ = @(F_,X_,colour_) plot(F_,X_/mean(X_), 'Color', cmap(colour_,:), 'LineWidth', 1);

%% -NMM--------------------------------------------------------------------

f1 = figure; hold
f2 = figure; hold

%% Down
colour_ = 1;
[x, ~, t, f_e, f_i, params, y] = NMM_GABA('u', 0, 'alpha_ri', 0.5, 'alpha_i', 1); 
figure(f1)
 plot(t(500:end),(y(500:end) - y(500))*1e3,'Color', cmap(colour_,:), 'LineWidth', 1);

[~, X_, F_] = fft_plot(y(500:end)-mean(y(500:end)), t(500:end),[],false);
figure(f2)
f_(F_,X_,colour_);

%% Normal
colour_ = 2;
[x, ~, t, f_e, f_i, params, y] = NMM_GABA('u', 0, 'alpha_ri', 1, 'alpha_i', 1); 
figure(f1)
 plot(t(500:end),(y(500:end) - y(500))*1e3,'Color', cmap(colour_,:), 'LineWidth', 1);
ylabel('LFP (mV)');
xlabel('Time (s)')

[~, X_, F_] = fft_plot(y(500:end)-mean(y(500:end)), t(500:end),[],false);
figure(f2)
f_(F_,X_,colour_);

%% Oscillation, can be called high amplitude/synchrony oscill
colour_ = 3;
[x, ~, t, f_e, f_i, params, y] = NMM_GABA('u', 0, 'alpha_ri', 2, 'alpha_i', 1.5); 
figure(f1)
 plot(t(500:end),(y(500:end) + 45e-3)*1e3,'Color', cmap(colour_,:), 'LineWidth', 1);

[~, X_, F_] = fft_plot(y(500:end)-mean(y(500:end)), t(500:end),[],false);
figure(f2)
f_(F_,X_,colour_);

%% HFO, can be called low amplitude oscill
colour_ = 4;
[x, ~, t, f_e, f_i, params, y] = NMM_GABA('u', 0, 'alpha_ri', 1.9, 'alpha_i', 1.8); 
figure(f1)
 plot(t(500:end),(y(500:end) - y(500))*1e3,'Color', cmap(colour_,:), 'LineWidth', 1);

[~, X_, F_] = fft_plot(y(500:end)-mean(y(500:end)), t(500:end),[],false);
figure(f2)
f_(F_,X_,colour_);

%% Saturation
colour_ = 5;
[x, ~, t, f_e, f_i, params, y] = NMM_GABA('u', 0, 'alpha_ri', 2, 'alpha_i', 0.5); 
figure(f1)
ylabel('LFP (mV)');
xlabel('Time (s)');
xlim([0.5 1]);
ylim([-70 350]);
ax = gca;
ax.FontSize = 12;

[~, X_, F_] = fft_plot(y(500:end)-mean(y(500:end)), t(500:end),[],false);
figure(f2)
f_(F_,X_,colour_);
xlim([0 100]);
xlabel('Frequency (Hz)');
ylabel('Normalized |X| (a.u.)');
ax = gca;
ax.YTick = [];
ax.FontSize = 12;


%% -CUBN-------------------------------------------------------------------
folder = 'C:\Users\artemios\Documents\Multiscale_Models_Data\2023\e_vs_i_fano\'; 
type_of_LIF = 'CUBN';
params = set_parameters('gabab');
cmap =[0 0 0; 0 0 1; 1 0 1; 1 0.8 0; 1 0 0];
f_ = @(F_,X_,colour_) plot(F_,X_/mean(X_), 'Color', cmap(colour_,:), 'LineWidth', 1);
f3 = figure; hold
f4 = figure; hold

% %% Down
colour_ = 1;

% load
d = dir([folder '*_e*.mat']); % Load all files with _e in the name
lif = load([folder '\' d(235).name]);

% get LFP
lif.i_pe = lif.i_pe(2501:10000);
lif.i_pi = lif.i_pi(2501:10000);
lif.i_ie = lif.i_ie(2501:10000);
lif.i_ii = lif.i_ii(2501:10000);

y = -(lif.i_pe - lif.i_pi)/params.g_m_P; % LFP
t = 2501*lif.lfp_dt : lif.lfp_dt : (2500 + numel(y)) * lif.lfp_dt;

% plot
figure(f3)
plot(t,(y - y(1))*1e3,'Color', cmap(colour_,:), 'LineWidth', 1);
ylabel('LFP (mV)');
xlabel('Time (s)');
xlim([0.5 1]);
ylim([0 500]);
ax = gca;
ax.FontSize = 12;

[~, X_, F_] = fft_plot(y-mean(y), t,[],false);
figure(f4)
f_(F_,X_,colour_);
xlim([0 100]);
xlabel('Frequency (Hz)');
ylabel('Normalized |X| (a.u.)');
ax = gca;
ax.YTick = [];
ax.FontSize = 12;



% %% Normal
colour_ = 2;

% load
d = dir([folder '*_e*.mat']); % Load all files with _e in the name
lif = load([folder '\' d(205).name]);

% get LFP
lif.i_pe = lif.i_pe(2501:10000);
lif.i_pi = lif.i_pi(2501:10000);
lif.i_ie = lif.i_ie(2501:10000);
lif.i_ii = lif.i_ii(2501:10000);

y = -(lif.i_pe - lif.i_pi)/params.g_m_P; % LFP
t = 2501*lif.lfp_dt : lif.lfp_dt : (2500 + numel(y)) * lif.lfp_dt;

% plot
figure(f3)
plot(t,(y - y(1))*1e3,'Color', cmap(colour_,:), 'LineWidth', 1);
ylabel('LFP (mV)');
xlabel('Time (s)');
xlim([0.5 1]);
ylim([0 500]);
ax = gca;
ax.FontSize = 12;

[~, X_, F_] = fft_plot(y-mean(y), t,[],false);
figure(f4)
f_(F_,X_,colour_);
xlim([0 100]);
xlabel('Frequency (Hz)');
ylabel('Normalized |X| (a.u.)');
ax = gca;
ax.YTick = [];
ax.FontSize = 12;


% %% LFO
colour_ = 3;

% load
d = dir([folder '*_e*.mat']); % Load all files with _e in the name
lif = load([folder '\' d(200).name]);

% get LFP
lif.i_pe = lif.i_pe(2501:10000);
lif.i_pi = lif.i_pi(2501:10000);
lif.i_ie = lif.i_ie(2501:10000);
lif.i_ii = lif.i_ii(2501:10000);

y = -(lif.i_pe - lif.i_pi)/params.g_m_P; % LFP
t = 2501*lif.lfp_dt : lif.lfp_dt : (2500 + numel(y)) * lif.lfp_dt;

% plor
figure(f3)
plot(t,(y - y(1))*1e3,'Color', cmap(colour_,:), 'LineWidth', 1);
ylabel('LFP (mV)');
xlabel('Time (s)');
% xlim([0.5 2]);
% ylim([-70 350]);
ax = gca;
ax.FontSize = 12;

[~, X_, F_] = fft_plot(y-mean(y), t,[],false);
figure(f4)
f_(F_,X_,colour_);
xlim([0 100]);
xlabel('Frequency (Hz)');
ylabel('Normalized |X| (a.u.)');
ax = gca;
ax.YTick = [];
ax.FontSize = 12;


% %% HFO
colour_ = 4;

% load
d = dir([folder '*_e*.mat']); % Load all files with _e in the name
lif = load([folder '\' d(277).name]);

% get LFP
lif.i_pe = lif.i_pe(2501:10000);
lif.i_pi = lif.i_pi(2501:10000);
lif.i_ie = lif.i_ie(2501:10000);
lif.i_ii = lif.i_ii(2501:10000);

y = -(lif.i_pe - lif.i_pi)/params.g_m_P; % LFP
t = 2501*lif.lfp_dt : lif.lfp_dt : (2500 + numel(y)) * lif.lfp_dt;

% plor
figure(f3)
plot(t,(y - y(1))*1e3,'Color', cmap(colour_,:), 'LineWidth', 1);
ylabel('LFP (mV)');
xlabel('Time (s)');
% xlim([0.5 2]);
ylim([-70 350]);
ax = gca;
ax.FontSize = 12;

[~, X_, F_] = fft_plot(y-mean(y), t,[],false);
figure(f4)
f_(F_,X_,colour_);
xlim([0 100]);
xlabel('Frequency (Hz)');
ylabel('Normalized |X| (a.u.)');
ax = gca;
ax.YTick = [];
ax.FontSize = 12;


