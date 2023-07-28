% Plots nonlinearities and results of the fit. Also plots results of fitting EPSP and IPSP
%

%% Plot nonlinearities

% Pyramidal, Naka-Rushton
[x, nonlinearity_post, max_firing_rate, ~, ~, ~, ~, potential_integral, firing_rates] = analyze_brunel_injected_current('Py', 'N');

fig = figure;
s_py = scatter(potential_integral, firing_rates, 5, 'x', 'MarkerEdgeColor', [0.0745 0.6235 1]);
hold
l_py = plot(x, nonlinearity_post, 'r', 'LineWidth', 3, 'Color', [0 0.4471 0.7412]); xlabel('Membrane potential (mV)');

% In, Naka-Rushton
[x, nonlinearity_post, max_firing_rate, ~, ~, ~, ~, potential_integral, firing_rates] = analyze_brunel_injected_current('In', 'N');

figure(fig)
s_in = scatter(potential_integral, firing_rates, 5, 'x', 'MarkerEdgeColor', [1 0.4118 0.1608]);
l_in = plot(x, nonlinearity_post, 'r', 'LineWidth', 3, 'Color', [0.8510 0.3255 0.0980]); xlabel('Membrane potential (mV)');

xlabel('Membrane potential (mV)');
ylabel('Firing rate (spikes/s)');
ax = gca;
ax.FontSize = 12;
l = legend([l_py l_in], {'Py' 'In'});
l.Location = 'best';
l.FontSize = 12;
xlim([-20 80])

%% Plot rmse comparison
fig2 = figure;
hold

[~, ~, ~, bestrmse, ~, gof] = analyze_brunel_injected_current('Py', 'S'); figure(fig2); bar(1, gof{bestrmse}.rmse,'EdgeColor', 'none', 'FaceColor', [0.72 0.27 1]);
[~, ~, ~, bestrmse, ~, gof] = analyze_brunel_injected_current('Py', 'L'); figure(fig2); bar(2, gof{bestrmse}.rmse,'EdgeColor', 'none', 'FaceColor', [0.30,0.75,0.93]);
[~, ~, ~, bestrmse, ~, gof] = analyze_brunel_injected_current('Py', 'G'); figure(fig2); bar(3, gof{bestrmse}.rmse,'EdgeColor', 'none', 'FaceColor', [0.47 0.67 0.19]);
[~, ~, ~, bestrmse, ~, gof] = analyze_brunel_injected_current('Py', 'N'); figure(fig2); bar(4, gof{bestrmse}.rmse,'EdgeColor', 'none', 'FaceColor', [0.93 0.69 0.13]);

[~, ~, ~, bestrmse, ~, gof] = analyze_brunel_injected_current('In', 'S'); figure(fig2); bar(6, gof{bestrmse}.rmse,'EdgeColor', 'none', 'FaceColor', [0.72 0.27 1]);
[~, ~, ~, bestrmse, ~, gof] = analyze_brunel_injected_current('In', 'L'); figure(fig2); bar(7, gof{bestrmse}.rmse,'EdgeColor', 'none', 'FaceColor', [0.30,0.75,0.93]);
[~, ~, ~, bestrmse, ~, gof] = analyze_brunel_injected_current('In', 'G'); figure(fig2); bar(8, gof{bestrmse}.rmse,'EdgeColor', 'none', 'FaceColor', [0.47 0.67 0.19]);
[~, ~, ~, bestrmse, ~, gof] = analyze_brunel_injected_current('In', 'N'); figure(fig2); bar(9, gof{bestrmse}.rmse,'EdgeColor', 'none', 'FaceColor', [0.93 0.69 0.13]);

l = legend({'Erf', 'Logistic', 'Gompertz', 'Naka-Rushton'});
l.Location = 'northwest';
l.FontSize = 12;
ax = gca;
ax.FontSize = 12;
ylabel('RMSE');
xticks([2 6]);
xticklabels({'Pyramidal' 'Inhibitory'});

%% EPSP and IPSP
synapses = {'pp', 'ip', 'pu', 'iu', 'pi', 'pb', 'ii'};

psp_fig = figure('Position', [918 417 805 420]);
hold
for i = 1:numel(synapses)
    [fitresult, T, psp] = fit_alpha(synapses{i});
    
    figure(psp_fig);
    handle_fit(i) = plot(fitresult);
    handle_fit(i).LineWidth = 1;
    handle(i) = plot(T,psp, '--', 'LineWidth', 2);%,'b--')
end
l = legend([handle], {'EPSP^{P}_{AMPA}' 'EPSP^{I}_{AMPA}' 'EPSP^{P}_{AMPA_{u}}' 'EPSP^{I}_{AMPA_u}' 'IPSP^{P}_{GABA_A}' 'IPSP^{P}_{GABA_B}' 'IPSP^{I}_{GABA_A}'});
l.FontSize = 12;
l.Location = 'eastoutside';
xlim([0 0.2]);
xlabel('Time (s)');
ylabel('Unitary PSP (mV)');
grid
ax = gca;
ax.FontSize = 12;