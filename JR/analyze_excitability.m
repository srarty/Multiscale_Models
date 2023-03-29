% Analyze excitability
%
% For NMM data, just call analyze excitability with the results of the NMM
%
% Syntaxis: 
%   [t_recovery, y_, t_] = analyze_excitability(y, t);
%   [t_recovery, y_, t_] = analyze_excitability(y, t, idx_stim, min_peak, trim, PLOT, NEG);
% 
%
% Inputs:
%   y           : input recording (LFP or Membrane potential).
%   t           : time vector corresponding to the 'y'.
%   idx_stim    : (optional) Index at which the stimulus pulse started, 489
%                   is the default.
%   min_peak    : (optional) Minimum peak height to find for 'y_' after the
%                   stimulus pulse. The exponential decay is fitted from 
%                   the peak found above min_peak.
%   trim        : (optional) trims 'y' and 't' by these initial number of 
%                   samples prior to the analysis.
%   PLOT        : (optional) Boolean, plots the results if true (default).
%   NEG         : (optional) Multiplies 'y' by -1 when normalizing if true.
%   
%
% Outputs:
%   t_recovery  : recovery time in ms
%   y_          : Normalized (if NEG) and trimmed recording
%   t_          : trimmed time vector
%
%
% For LIF data, use:
%   data_file = 'C:/Users/artemios/Documents/Multiscale_Models_Data/lfp_86.mat'; % 500pA current pulse, u=0
%   data_file = 'C:/Users/artemios/Documents/Multiscale_Models_Data/lfp_89.mat'; % 50pA current pulse, u=0
%   data_file = 'C:\Users\artemios\Documents\Multiscale_Models_Data\2023\excitability\lfp_8';
%   a = load(data_file)
%   y_lif = a.LFP_V;
%   t_lif = a.lfp_dt:a.lfp_dt:length(a.LFP_V)*a.lfp_dt;
%   recovery = analyze_excitability(y_lif', t_lif', 4899, 0)
%
function [recovery, y_, t_] = analyze_excitability(y,t,varargin)

if nargin >= 3, idx_stim = varargin{1}; else, idx_stim = 489; end
if nargin >= 4, min_peak = varargin{2}; else, min_peak = 0.5; end
if nargin >= 5, trim = varargin{3}; else, trim = 1; end
if nargin >= 6, PLOT = varargin{4}; else, PLOT = true; end
if nargin >= 7, NEG = varargin{5}; else, NEG = true; end

% Trim the first part of the simulation
y = y(trim:end);
idx_stim = idx_stim - trim;

% Normalize signal
y_ = -y/abs(max(-y));

% Find peak. The time between idx_stim and the peak (plus 20 ms) will be 
% the approximate recovery time given that the decay time constant from 
% peak to recovery is the same in all recordings and it is almost identical
% to the membrane time constant of the pyramidal population.
try
    [~, idx] = findpeaks(y_, 'MinPeakHeight', min_peak);
catch E
    disp('MinPeakHeight higher than number of samples')
end
idx(idx <= idx_stim+1) = []; % Remove peaks found before the impulse
idx = min(idx); % only leave smallest idx

y_ = y_(idx:end);
try
    t_ = t(trim:end-idx+1);
    t_ = t_ - t_(1);
catch E
    disp('No peaks found');
    f = figure; plot(t(trim:end), y);
%     response = questdlg('No peaks found, what happened?', 'Manual input', 'Saturation', 'Oscillation', 'Oscillation');
    response = 'Oscillation';
    switch response, case 'Oscillation', failure_value = 0; case 'Saturation', failure_value = Inf; end
    close(f);
end
    
if ~isempty(y_)
    % If y_ is not empty, idx should have found a peak, here we're
    % measuring the decay time constant by fitting the decay to a single
    % exponential. This helps us estimate the total recovery time at the
    % end.
    ft = fittype( 'a*exp(-t/tau) - b', 'independent', 't', 'dependent', 'y');
    opts = fitoptions(ft);
    opts.StartPoint = [1 0 0.1];
    opts.Lower = [0 -5 0];
    opts.Upper = [2000 5 10];
    opts.Robust = 'Off';
    fitresult = fit(t_, y_, ft, opts); % With options
    % fitresult = fit(t_, y_, ft); % No options
    
    fit_time_constant = fitresult.tau; % In seconds
else
    % If y_ is empty, we assume the model saturated, so there's no recovery
    % and we set the recovefry time to infinity.
    try
        fit_time_constant = failure_value;
    catch E
        disp('''failure_value'' should have been set in a dialog question. If the question didn''t appear, you need to debug.');
    end
end

if PLOT
    figure; 
    plot(t_, y_);
    hold;
    plot(fitresult);
end

% Recovery time:
dt = t(2)-t(1);
recovery = (6 * fit_time_constant) + dt*(max(idx) - idx_stim); % Multiply by 6 because we calculated the time constant and now we want the actual decay time

% If large negative number, it means (most likely) that the populations are
% saturated after the probe
if isempty(recovery) || (recovery < 0), recovery = nan; end