% Analyze excitability
%
% For NMM data, just call analyze excitability with the results of the NMM
%
% For LIF data, use:
%   data_file = 'C:/Users/artemios/Documents/Multiscale_Models_Data/lfp_79.mat'; % impulse response (500 pA), j_pi = 37, alpha_i = -0.5xxx
%   a = load(data_file)
%   y_lif = a.LFP_V;
%   t_lif = a.lfp_dt:a.lfp_dt:length(a.LFP_V)*a.lfp_dt;
%   recovery = analyze_excitability(y_lif', t_lif', 4899, -4e14)
%
function [recovery, y_, t_] = analyze_excitability(y,t,varargin)

if nargin >= 3, idx_stim = varargin{1}; else, idx_stim = 489; end
if nargin >= 4, min_peak = varargin{2}; else, min_peak = 0.5; end
   

% Normalize signal
y_ = -y/(max(-y));
% Find peak. The time between idx_stim and the peak (plus 20 ms) will be 
% the approximate recovery time given that the decay time constant from 
% peak to recovery is the same in all recordings and it is almost identical
% to the membrane time constant of the pyramidal population.
[~, idx] = findpeaks(y_, 'MinPeakHeight', min_peak);
idx(idx <= idx_stim+1) = []; % Remove peaks found before the impulse
idx = min(idx); % only leave smallest idx

y_ = y_(idx:end);
t_ = t(1:end-idx+1);

if ~isempty(y_)
    % If y_ is not empty, idx should have found a peak, here we're
    % measuring the decay time constant by fitting the decay to a single
    % exponential. This helps us estimate the total recovery time at the
    % end.
    ft = fittype( 'a*exp(-t/tau) - b', 'independent', 't', 'dependent', 'y');
    opts = fitoptions(ft);
    opts.StartPoint = [1 0 0.1];
    opts.Lower = [0 -1 0];
    opts.Upper = [300 1 1];
    opts.Robust = 'Off';
    fitresult = fit(t_, y_, ft, opts); % With options
    % fitresult = fit(t_, y_, ft); % No options
    
    fit_time_constant = fitresult.tau * 1000;
else
    % If y_ is empty, we assume the model saturated, so there's no recovery
    % and we set the recovefry time to infinity.
    fit_time_constant = Inf;
end

% 
% figure; 
% plot(t_, y_);
% hold;
% plot(fitresult);

% Recovery time:
recovery = fit_time_constant + max(idx) - idx_stim;
% integral = sum()% todo

% If large negative number, it means (most likely) that the populations are
% saturated after the probe
if isempty(recovery) || (recovery < 0), recovery = nan; end