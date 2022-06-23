% EXPLORE_VALUES runs NMM_diff_equations_DblExp.m with different values
% of amplitude and/or time constants to find inflexion points, etc
%
% Artemio - May 2022

%% Run NMM and store values
% range = [-1e7:0.5e6:-0.1e6]; % Range for AmplitudeI
% range = [5e6:0.2e6:5e7]; % Range for AmplitudeE
% range2 = [0.5e7:0.5e6:2e7]; % Range for AmplitudeE
% range = [0:5:300]; % u
range = [15:25]; % u
range2 = 1;
% range = [0.002:0.001:0.05]; % tau_m_e
% range2 = [0.005:0.001:0.055]; % tau_m_i
freqs = [];
for i = 1:length(range)
    for ii = 1:length(range2)
%     [x, y, t] = NMM_diff_equations_DblExp('AmplitudeE', range(i));
%         [x, y, t] = NMM_diff_equations_DblExp('AmplitudeI', range(i), 'AmplitudeE', 1.42e7);
%         [x, y, t] = NMM_diff_equations_DblExp('AmplitudeI', range(i), 'AmplitudeE', range2(ii));
        [x, y, t] = NMM_diff_equations_DblExp('u', range(i));
%         [x, y, t] = NMM_diff_equations_DblExp('tau_m_e', range(i), 'tau_m_i', range2(ii));
        freqs(i,ii) = spectrum(x,y,t, false);
        disp([num2str(i) '/' num2str(length(range)) ' , ' num2str(ii) '/' num2str(length(range2))]);
    end
end

%% Plot 3d Mesh
figure;
mesh(range2, range, freqs, 'FaceColor', 'flat', 'EdgeColor', 'black')
xlabel('\tau_{m_{i}}');
ylabel('\tau_{m_{e}}');
zlabel('Frequency (Hz)');
hold
% plot3(0.01638,0.008115,25.6339,'rx','LineWidth',3)
% title('NMM | u = 15 spikes/ms/cell');
title('NMM');
c = colorbar;
c.Label.String = 'Frequency (Hz)';


%% Load LIF results in a loop and store values
location = 'C:/Users/artemios/Documents/Multiscale_Models_Data/spartan/';
d = dir([location '*.mat']);

freqs = []; u = [];
for i = 1: length(d)
    [freqs(i) u(i)] = spectrum(x,y,t,false,[d(i).folder '\' d(i).name]);
    disp([num2str(i) '/' num2str(length(d))])
end

%%
figure;
stem(u,freqs)
xlabel('Input spike rate')
ylabel('Frequency (Hz)')