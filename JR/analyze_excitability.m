% Analyze excitability of NMM

y_ = -y/(max(-y));
[a, idx] = findpeaks(y_, 'MinPeakHeight', 0.5);

y_ = y_(idx:end);
t_ = t(1:end-idx+1);

ft = fittype( 'a*exp(-t/tau) - b', 'independent', 't', 'dependent', 'y');
opts = fitoptions(ft);
opts.StartPoint = [1 0 0.1];
opts.Lower = [0 -1 0];
opts.Upper = [300 1 1];
opts.Robust = 'Off';
fitresult = fit(t_, y_, ft, opts) % With options
% fitresult = fit(t_, y_, ft) % No options


figure; 
plot(t_, y_);
hold
plot(fitresult)