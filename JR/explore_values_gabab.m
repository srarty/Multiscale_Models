% EXPLORE_VALUES_GABAB runs NMM_GABA.m with different values
% of parameters
%
% Plots two varying parameters in a surface plot.
%
% Artemio - May 2022

%% Run NMM and store values
% range = [-1e7:0.5e6:-0.1e6]; % Range for AmplitudeI
% range = [5e6:0.2e6:5e7]; % Range for AmplitudeE
% range2 = [0.5e7:0.5e6:2e7]; % Range for AmplitudeE
% range = [15:25]; % u
% range = [0.002:0.001:0.05]; % tau_m_e
% range2 = [0.005:0.001:0.055]; % tau_m_i

value = 'alpha_i';
range = 0:0.1:2;%0:0.1:2;%-1*[0:0.05:3];
value2 = 'u';%'u';
range2 = 0:0.1:2;%0:0.1:2;%[0:0.05:3];%[0:0.05:1];%

freqs = [];
freqs_py = [];
freqs_in = [];
w_nmm = [];
w_lif = [];
recovery = [];%size(range);
for i = 1:length(range)
    for ii = 1:length(range2)
%         [x, y, t] = NMM_diff_equations_DblExp_recursive(value, range(i));
        [x, y, t, f_e, f_i] = NMM_GABA(value, range(i), value2, range2(ii));

%         freqs(ii,i) = spectrum(x,y,t, false); % Oscillations
%         try
%             [w_nmm(ii,i), w_lif(ii,i)] = spectrum(x,y,t, false); % Width of Autocorrelation function
%         catch E
%             if strcmp('MATLAB:subsassigndimmismatch', E.identifier)
%                 w_nmm = 2000;
%             else
%                 rethrow(E);
%             end
%         end


        freqs_py(ii,i) = mean(f_e(250:end)); % spike rate Py
        freqs_in(ii,i) = mean(f_i(250:end)); % spike rate Py

        disp([num2str(i) '/' num2str(length(range)) ' , ' num2str(ii) '/' num2str(length(range2))]);

%         recovery(i,ii) = analyze_excitability(y,t);
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

folder = 'C:\Users\artemios\Documents\Multiscale_Models_Data\Spartan\firing_rates_nmm';
name = [value ' vs ' value2];
if isempty(dir([folder name]))
    save([folder name '.mat'], 'results');
    disp('results saved');
else
    disp('Results not saved, file exists');
end



% %% Plot 3d Mesh for oscillations
% % z_axis = freqs; z_label = 'Frequency (Hz)'; limits = [25 100];
% z_axis = w_nmm; z_label = 'Autocorrelation Width'; limits = [25 100];
% % z_axis = recovery; z_label = 't_{recovery}'; limits = [15 20];
% try
%     f_handle = figure;
%     mesh(range, range2, z_axis, 'FaceColor', 'flat', 'EdgeColor', 'none')
%     % xlabel('\tau_{m_{i}}');
%     % ylabel('\tau_{m_{e}}');
%     xlabel(value);%('Input rate');
%     ylabel(value2);
%     zlabel(z_label);
%     hold
%     % plot3(0.01638,0.008115,25.6339,'rx','LineWidth',3)
% %     title('NMM | u = 9');
%     title('NMM');
%     colormap cool
%     c = colorbar;
%     c.Label.String = z_label;
%     caxis(limits);
%     c.Limits = limits;
%     zlim(limits);
% catch ME
%     if strcmp('MATLAB:surfchk:NonMatrixData', ME.identifier)
%         plot(range,z_axis)
%     else
%         close(f_handle);
%         disp(['Couldn''t plot oscillation frequencies graph: ' ME.message]);
%     end
% end
% 
% %% Plot value vs recovery time
% figure; stem(range, recovery);
% xlabel(['\' results.value]);
% ylabel('Recovery time (ms)');
% 
% return

%% Plot spike rates mesh
% if results.range2(end)<0, angle=[0 -90]; else, angle=[0 90]; end
angle = [0 -90];

f = figure('Position', [325 404 1112 364]);
f.Position = [589 411 758 228];
ax = subplot(1,2,1);

imagesc(results.range, results.range2, results.freqs_py);% , 'FaceColor', 'flat', 'EdgeColor', 'none')
xlabel(['\' results.value]);%('Input rate');
% ylabel(['\' results.value2]);
ylabel(results.value2);
zlabel('Firing rate (Py)');

cmap = colormap(jet(256)); % use jet colormap as a starting point
cmap(end,:) = [0 0 0]; % set values greater than 10 to black

% apply custom colormap
colormap(cmap);
c = colorbar;
c.Label.String = 'Mean firing rate (Hz)';
caxis([0 0.6]);
c.Limits = [0 0.6];
% zlim([0 1]);
title('Pyramidal')
%xlim([0 5])
ax.View = (angle);
ax.FontSize = 12;

ax = subplot(1,2,2);
imagesc(results.range, results.range2, results.freqs_in);% , 'FaceColor', 'flat', 'EdgeColor', 'none')
xlabel(['\' results.value]);%('Input rate');
% ylabel(['\' results.value2]);
ylabel(results.value2);
zlabel('Firing rate (In)');
title('Inhibitory');
c = colorbar;
c.Label.String = 'Mean firing rate (Hz)';
caxis([0 6]);
c.Limits = [0 6];
% zlim([0 2]);
% xlim([0 5])
ax.View = (angle);
ax.FontSize = 12;

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