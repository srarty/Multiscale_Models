
% [x,y,t,f_e,f_i, params, yy]=NMM_GABA();


close all


milivolts_scale = 1e-3; % The results are in mv, this scale is to express them in volts

% nmm_i_pi = (x(:,2) + x(:,10))*milivolts_scale * params.C_P;
% nmm_i_pe = (x(:,6) + x(:,8))*milivolts_scale * params.C_P;
% nmm_i_ie = (x(:,4) + x(:,20))*milivolts_scale * params.C_I;
% nmm_i_ii = (x(:,12))* milivolts_scale * params.C_I;
nmm_i_pi = (x(:,1) + x(:,9)) * milivolts_scale * params.g_m_P;
nmm_i_pe = (x(:,5) + x(:,7)) * milivolts_scale * params.g_m_P;
nmm_i_ie = (x(:,3) + x(:,19)) * milivolts_scale * params.g_m_I;
nmm_i_ii = (x(:,11)) * milivolts_scale * params.g_m_I;

%%
% load('C:\Users\artemios\Documents\Multiscale_Models_Data\2023\lfp_last.mat')
% Los datos de la derecha:
% load('C:\Users\artemios\Documents\Multiscale_Models_Data\Spartan\firing_rates\lfp_py_j_GABAb_12_u0.6000000000000001.mat');
load('C:\Users\artemios\Documents\Multiscale_Models_Data\Spartan\firing_rates\lfp_py_j_GABAb_1_u2.0.mat')


% PLOT
figure; 
plot([0 5], [0 0], '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 2); 
hold on

% LIF
l1 = plot(1, mean(i_pe)*1e9, 'b^', 'MarkerSize', 8, 'LineWidth', 2, 'MarkerFaceColor', 'none'); 
l2 = plot(1, mean(i_pi)*1e9, 'r^', 'MarkerSize', 8, 'LineWidth', 2, 'MarkerFaceColor', 'none');
l3 = plot(1, (mean(i_pi) + mean(i_pe))*1e9, 'k^', 'MarkerSize', 8, 'LineWidth', 2, 'MarkerFaceColor', 'none');
plot(1, mean(i_ie)*1e9, 'bo', 'MarkerSize', 8, 'LineWidth', 2, 'MarkerFaceColor', 'none');
plot(1, mean(i_ii)*1e9, 'ro', 'MarkerSize', 8, 'LineWidth', 2, 'MarkerFaceColor', 'none');
plot(1, (mean(i_ii) + mean(i_ie))*1e9, 'ko', 'MarkerSize', 8, 'LineWidth', 2, 'MarkerFaceColor', 'none');
%%
% NMM
L = numel(nmm_i_pe);
scale = 1e9;
plot(2,-mean(nmm_i_pe(round(L/2):end))*scale, 'b^', 'MarkerSize', 8, 'LineWidth', 2, 'MarkerFaceColor', 'none');
plot(2,-mean(nmm_i_pi(round(L/2):end))*scale, 'r^', 'MarkerSize', 8, 'LineWidth', 2, 'MarkerFaceColor', 'none');
plot(2,-(mean(nmm_i_pi(round(L/2):end)) + mean(nmm_i_pe(round(L/2):end)))*scale, 'k^', 'MarkerSize', 8, 'LineWidth', 2, 'MarkerFaceColor', 'none');
l4 = plot(2,-mean(nmm_i_ie(round(L/2):end))*scale, 'bo', 'MarkerSize', 8, 'LineWidth', 2, 'MarkerFaceColor', 'none');
l5 = plot(2,-mean(nmm_i_ii(round(L/2):end))*scale, 'ro', 'MarkerSize', 8, 'LineWidth', 2, 'MarkerFaceColor', 'none');
l6 = plot(2,-(mean(nmm_i_ii(round(L/2):end)) + mean(nmm_i_ie(round(L/2):end)))*scale, 'ko', 'MarkerSize', 8, 'LineWidth', 2, 'MarkerFaceColor', 'none');

xlim([0.5 2.5])
l=legend([l1 l2 l3 l4 l5 l6],{'E_p' 'I_p' '(E + I)_p' 'E_i' 'I_i' '(E + I)_i'});
l.Location = 'eastoutside';
xticks([1 2]);
% xticklabels({'Py_{LIF}' 'In_{LIF}' 'Py_{NMM}' 'In_{NMM}'});
xticklabels({'LIF' 'NMM'});
ylabel('input current (nA)');
ax = gca;
ax.FontSize = 12;
