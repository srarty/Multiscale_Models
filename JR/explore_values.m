% EXPLORE_VALUES runs NMM_diff_equations_DblExp.m with different values
% of amplitude and/or time constants to find inflexion points, etc
%
% Artemio - May 2022

%% Run NMM and store values
% range = [-1e7:0.5e6:-0.1e6]; % Range for AmplitudeI
% range = [5e6:0.2e6:5e7]; % Range for AmplitudeE
% range2 = [0.5e7:0.5e6:2e7]; % Range for AmplitudeE
% range = [15:25]; % u
% range = [0.002:0.001:0.05]; % tau_m_e
% range2 = [0.005:0.001:0.055]; % tau_m_i

value = 'c_constant'; %'e0';%'pII';
range = logspace(-2,3,100)%[0.1:1:100]; % P[II] or P[PP]
value2 = 'u';%'e0i';
range2 = 1;%1*[10:5:100];%[0:0.5:5]; %[0:10]; % P[II] or P[PP]

freqs = [];
freqs_py = [];
freqs_in = [];
for i = 1:length(range)
    for ii = 1:length(range2)
%         [x, y, t] = NMM_diff_equations_DblExp_recursive(value, range(i));
        [x, y, t, f_e, f_i] = NMM_diff_equations_DblExp_recursive(value, range(i), value2, range2(ii));

%         freqs(ii,i) = spectrum(x,y,t, false); % Oscillations
        freqs_py(ii,i) = mean(f_e(250:end)); % spike rate Py
        freqs_in(ii,i) = mean(f_i(250:end)); % spike rate Py
        disp([num2str(i) '/' num2str(length(range)) ' , ' num2str(ii) '/' num2str(length(range2))]);
    end
end

%% Store values
results = struct;
results.value = value;
results.range = range;
results.value2 = value2;
results.range2 = range2;
results.freqs = freqs;
results.freqs_py = freqs_py;
results.freqs_in = freqs_in;

folder = 'C:\Users\artemios\Dropbox\University of Melbourne\Epilepsy\Resources for meetings\2022 07 14\';
name = [value ' vs ' value2];
if isempty(dir([folder name]))
    save([folder name '.mat'], 'results');
    disp('results saved');
else
    disp('Results not saved, file exists');
end

%% Plot 3d Mesh for oscillations
try
    f_handle = figure;
    mesh(range, range2, freqs, 'FaceColor', 'flat', 'EdgeColor', 'black')
    % xlabel('\tau_{m_{i}}');
    % ylabel('\tau_{m_{e}}');
    xlabel(value);%('Input rate');
    ylabel(value2);
    zlabel('Frequency (Hz)');
    hold
    % plot3(0.01638,0.008115,25.6339,'rx','LineWidth',3)
    title('NMM | u = 9');
    % title('NMM');
    colormap hsv
    c = colorbar;
    c.Label.String = 'Frequency (Hz)';
    caxis([25 100]);
    c.Limits = [25 65];
    zlim([25 100]);
catch
    close(f_handle);
    disp('Couldn''t plot oscillation frequencies graph.');
end

%% Plot spike rates mesh
if results.range2(end)<0, angle=[0 -90]; else, angle=[0 90]; end

figure('Position', [325 404 1112 364]);
ax = subplot(1,2,1);

mesh(results.range, results.range2, results.freqs_py, 'FaceColor', 'flat', 'EdgeColor', 'none')
xlabel(results.value);%('Input rate');
ylabel(results.value2);
zlabel('Firing rate (Py)');
c = colorbar;

c.Label.String = 'Mean firing rate (Hz)';
caxis([0 1.5]);
c.Limits = [0 1.5];
% zlim([0 1.5]);
title('Pyramidal')
%xlim([0 5])
ax.View = (angle);

ax = subplot(1,2,2);
mesh(results.range, results.range2, results.freqs_in, 'FaceColor', 'flat', 'EdgeColor', 'none')
xlabel(results.value);%('Input rate');
ylabel(results.value2);
zlabel('Firing rate (In)');
title('Inhibitory');
colormap jet
c = colorbar;
c.Label.String = 'Mean firing rate (Hz)';
caxis([0 3.5]);
c.Limits = [0 3.5];
% zlim([0 1.5]);
% xlim([0 5])
ax.View = (angle);

%%
pause(0.1);
%{
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
% stem(u,freqs)
% xlabel('Input spike rate')
stem(range,freqs)
xlabel('Probability of recurrent inhibitory synapses')
ylabel('Frequency (Hz)')
%}