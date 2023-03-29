load('C:\Users\artemios\Documents\Multiscale_Models_Data\2023\lfp_last.mat')
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

% PLOT
figure; 
plot([0 5], [0 0], '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 2); 
hold on

% LIF
l1 = plot(1, mean(i_pe)*1e9, 'bo', 'MarkerSize', 8, 'LineWidth', 3, 'MarkerFaceColor', 'b'); 
l2 = plot(1, mean(i_pi)*1e9, 'ro', 'MarkerSize', 8, 'LineWidth', 3, 'MarkerFaceColor', 'r');
l3 = plot(1, (mean(i_pi) + mean(i_pe))*1e9, 'ko', 'MarkerSize', 8, 'LineWidth', 3, 'MarkerFaceColor', 'k');
errorbar(2, mean(i_ie)*1e9, 0*std(i_ie)*1e9, 'bo', 'MarkerSize', 8, 'LineWidth', 3, 'MarkerFaceColor', 'auto');
errorbar(2, mean(i_ii)*1e9, 0*std(i_ii)*1e9, 'ro', 'MarkerSize', 8, 'LineWidth', 3, 'MarkerFaceColor', 'auto');
errorbar(2, (mean(i_ii) + mean(i_ie))*1e9, (0*std(i_ii + i_ie))*1e9, 'ko', 'MarkerSize', 8, 'LineWidth', 3, 'MarkerFaceColor', 'auto');

% NMM
L = numel(nmm_i_pe);
scale = 1e9;
errorbar(3,-mean(nmm_i_pe(round(L/2):end))*scale ,-0*std(nmm_i_pe(round(L/2):end))*scale, 'bo', 'MarkerSize', 8, 'LineWidth', 3, 'MarkerFaceColor', 'auto');
errorbar(3,-mean(nmm_i_pi(round(L/2):end))*scale ,-0*std(nmm_i_pi(round(L/2):end))*scale, 'ro', 'MarkerSize', 8, 'LineWidth', 3, 'MarkerFaceColor', 'auto');
errorbar(3,-(mean(nmm_i_pi(round(L/2):end)) + mean(nmm_i_pe(round(L/2):end)))*scale ,-(0*std(nmm_i_pi(round(L/2):end) + nmm_i_pe(round(L/2):end)))*scale, 'ko', 'MarkerSize', 8, 'LineWidth', 3, 'MarkerFaceColor', 'auto');
errorbar(4,-mean(nmm_i_ie(round(L/2):end))*scale ,-0*std(nmm_i_ie(round(L/2):end))*scale, 'bo', 'MarkerSize', 8, 'LineWidth', 3, 'MarkerFaceColor', 'auto');
errorbar(4,-mean(nmm_i_ii(round(L/2):end))*scale, -0*std(nmm_i_ii(round(L/2):end))*scale, 'ro', 'MarkerSize', 8, 'LineWidth', 3, 'MarkerFaceColor', 'auto');
errorbar(4,-(mean(nmm_i_ii(round(L/2):end)) + mean(nmm_i_ie(round(L/2):end)))*scale ,-(0*std(nmm_i_ii(round(L/2):end) + nmm_i_ie(round(L/2):end)))*scale, 'ko', 'MarkerSize', 8, 'LineWidth', 3, 'MarkerFaceColor', 'auto');

xlim([0.5 4.5])
l=legend([l1 l2 l3],{'E' 'I' 'E + I'});
l.Location = 'best';
xticks([1 2 3 4]);
xticklabels({'Py_{LIF}' 'In_{LIF}' 'Py_{NMM}' 'In_{NMM}'});
ylabel('input current (nA)');
ax = gca;
ax.FontSize = 12;
