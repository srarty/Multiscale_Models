% EXPLORE_VALUES_3POPS runs NMM_GABAb.m with different values
% of amplitude and/or time constants to find inflexion points, etc
%
% Artemio - Nov 2022

%% Run NMM

values = 1:0.1:2;
recovery_times = zeros(3, length(values));
w = zeros(3, length(values));
rates = zeros(3, length(values));
for i = 1:length(values)
    disp(['Loop ' num2str(i) ' out of ' num2str(length(values))]);
    
    [x, y, t, f_e, f_i, f_b]=NMM_GABAb('alpha_i', values(i), 'alpha_ri', values(i)); % Benzodiazepine
    rates(1,i) = mean(f_e(ceil(end/2):end));
    rates(2,i) = mean(f_i(ceil(end/2):end));
    rates(3,i) = mean(f_b(ceil(end/2):end));
end

%% Plot firing rates
figure; plot(values, 100 * rates(1,:)/ rates(1,1));
hold
plot(values,  100 * rates(2,:)/ rates(2,1));
plot(values,  100 * rates(3,:)/ rates(3,1));
legend({'Pyramidal' 'GABA_A', 'GABA_B'})
ylabel('Firing rate (%)')
xlabel('GABA_A boost')

%% Plot 3d Mesh for oscillations
% z_axis = freqs; z_label = 'Frequency (Hz)'; limits = [25 100];
z_axis = w_nmm; z_label = 'Autocorrelation Width'; limits = [25 100];
% z_axis = recovery; z_label = 't_{recovery}'; limits = [15 20];
try
    f_handle = figure;
    mesh(range, range2, z_axis, 'FaceColor', 'flat', 'EdgeColor', 'none')
    % xlabel('\tau_{m_{i}}');
    % ylabel('\tau_{m_{e}}');
    xlabel(value);%('Input rate');
    ylabel(value2);
    zlabel(z_label);
    hold
    % plot3(0.01638,0.008115,25.6339,'rx','LineWidth',3)
%     title('NMM | u = 9');
    title('NMM');
    colormap cool
    c = colorbar;
    c.Label.String = z_label;
    caxis(limits);
    c.Limits = limits;
    zlim(limits);
catch ME
    if strcmp('MATLAB:surfchk:NonMatrixData', ME.identifier)
        plot(range,z_axis)
    else
        close(f_handle);
        disp(['Couldn''t plot oscillation frequencies graph: ' ME.message]);
    end
end

%% Plot value vs recovery time
figure; stem(range, recovery);
xlabel(['\' results.value]);
ylabel('Recovery time (ms)');

return

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
caxis([0 2]);
c.Limits = [0 2];
% zlim([0 1]);
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
caxis([0 2]);
c.Limits = [0 2];
% zlim([0 2]);
% xlim([0 5])
ax.View = (angle);

%%
%{
pause(0.1);
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