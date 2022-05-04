% EXPLORE_VALUES runs NMM_diff_equations_DblExp.m with different values
% of amplitude and/or time constants to find inflexion points, etc
%
% Artemio - May 2022

% range = [-1e7:0.5e6:-0.1e6]; % Range for AmplitudeI
% range = [5e6:0.2e6:5e7]; % Range for AmplitudeE
% range2 = [0.5e7:0.5e6:2e7]; % Range for AmplitudeE
range = [5:30];
range2 = 1;
freqs = []; amps = [];
for i = 1:length(range)
    for ii = 1:length(range2)
%     [x, y, t] = NMM_diff_equations_DblExp('AmplitudeE', range(i));
%         [x, y, t] = NMM_diff_equations_DblExp('AmplitudeI', range(i), 'AmplitudeE', 1.42e7);
%         [x, y, t] = NMM_diff_equations_DblExp('AmplitudeI', range(i), 'AmplitudeE', range2(ii));
        [x, y, t] = NMM_diff_equations_DblExp('u', range(i));
        freqs(i,ii) = spectrum(x,y,t, false);        
        disp([num2str(i) ' , ' num2str(ii) '/' num2str(length(range2))]);
    end
end