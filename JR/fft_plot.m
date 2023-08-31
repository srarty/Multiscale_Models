%% FFT quickly plots fft of signal
%
function varargout = fft_plot(x, t, varargin)

if nargin > 2, fig = varargin{1}; else, fig = []; end
if nargin > 3, PLOT = varargin{2}; else, PLOT = false; end

% If column vector, transform to row.
if size(x,1) > size(x,2), x = transpose(x); end

Fs = 1/diff(t(1:2));
L = length(t);
n = 2^nextpow2(L);
X = fft(x,n);
P2 = abs(X/L);
P1 = P2(:,1:n/2+1);
P1(:,2:end-1)=2*P1(:,2:end-1);

if PLOT
    if ~isempty(fig), figure(fig); else, fig = figure; end
    try
        if ~ishold(fig.Children), hold; end 
    catch
        disp('Couldn''t hold fig');
    end
    plot(0:(Fs/n):(Fs/2-Fs/n), P1(1:n/2));
%     plot(0:(Fs/(0.5*n)):(Fs), P1); % This line works for NMM, tested it counting the actual number of peaks.
    xlim([0 200]);
end

% Need to make sure the following lines make sense for LIF. They have been
% tested on NMM
X_ = P1(1:n/2);
F_ = 0:(Fs/n):(Fs/2-Fs/n);
% X_ = P1; % Tested for NMM
% F_ = 0:(Fs/(0.5*n)):(Fs); % Tested for NMM
varargout = {fig, X_, F_};

end