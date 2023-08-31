
% [x,y,t,f_e,f_i, params, yy]=NMM_GABA();


%close all

function Figure_balance(x, params, lif, varargin)

nanoamps_scale = 1e9;
milivolts_scale = 1e-3; % The results are in mv, this scale is to express them in volts

% nmm_i_pi = (x(:,2) + x(:,10))*milivolts_scale * params.C_P;
% nmm_i_pe = (x(:,6) + x(:,8))*milivolts_scale * params.C_P;
% nmm_i_ie = (x(:,4) + x(:,20))*milivolts_scale * params.C_I;
% nmm_i_ii = (x(:,12))* milivolts_scale * params.C_I;
nmm_i_pi = (x(:,1) + x(:,9)) * milivolts_scale * params.g_m_P;
nmm_i_pe = (x(:,5) + x(:,7)) * milivolts_scale * params.g_m_P;
nmm_i_ie = (x(:,3) + x(:,19)) * milivolts_scale * params.g_m_I;
nmm_i_ii = (x(:,11)) * milivolts_scale * params.g_m_I;

L = 2;%numel(nmm_i_pe);

%%
% load('C:\Users\artemios\Documents\Multiscale_Models_Data\2023\lfp_last.mat')
% Los datos de la derecha:
% load('C:\Users\artemios\Documents\Multiscale_Models_Data\Spartan\firing_rates\lfp_py_j_GABAb_12_u0.6000000000000001.mat');
% load('C:\Users\artemios\Documents\Multiscale_Models_Data\Spartan\firing_rates\lfp_py_j_GABAb_1_u2.0.mat')


% PLOT
figure; 
plot([0 5], [0 0], '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 2); 
hold on

% Pyramidal
l1 = plot(1, mean(lif.i_pe)*1e9, 'bx', 'MarkerSize', 14, 'LineWidth', 2, 'MarkerFaceColor', 'none'); 
l2 = plot(1, mean(lif.i_pi)*1e9, 'rx', 'MarkerSize', 14, 'LineWidth', 2, 'MarkerFaceColor', 'none');
l3 = plot(1, (mean(lif.i_pi) + mean(lif.i_pe))*1e9, 'kx', 'MarkerSize', 14, 'LineWidth', 2, 'MarkerFaceColor', 'none');
if nargin == 5
    lif_co = varargin{2};
    plot(1, mean(lif_co.i_pe)*1e9, 'b+', 'MarkerSize', 14, 'LineWidth', 2, 'MarkerFaceColor', 'none'); 
    plot(1, mean(lif_co.i_pi)*1e9, 'r+', 'MarkerSize', 14, 'LineWidth', 2, 'MarkerFaceColor', 'none');
    l7=plot(1, (mean(lif_co.i_pi) + mean(lif_co.i_pe))*1e9, 'k+', 'MarkerSize', 14, 'LineWidth', 2, 'MarkerFaceColor', 'none');
end
plot(1,-mean(nmm_i_pe(round(L/2):end))*nanoamps_scale, 'bo', 'MarkerSize', 10, 'LineWidth', 2, 'MarkerFaceColor', 'none');
plot(1,-mean(nmm_i_pi(round(L/2):end))*nanoamps_scale, 'ro', 'MarkerSize', 10, 'LineWidth', 2, 'MarkerFaceColor', 'none');
plot(1,-(mean(nmm_i_pi(round(L/2):end)) + mean(nmm_i_pe(round(L/2):end)))*nanoamps_scale, 'ko', 'MarkerSize', 10, 'LineWidth', 2, 'MarkerFaceColor', 'none');

% Interneurons
plot(2, mean(lif.i_ie)*1e9, 'bx', 'MarkerSize', 14, 'LineWidth', 2, 'MarkerFaceColor', 'none');
plot(2, mean(lif.i_ii)*1e9, 'rx', 'MarkerSize', 14, 'LineWidth', 2, 'MarkerFaceColor', 'none');
plot(2, (mean(lif.i_ii) + mean(lif.i_ie))*1e9, 'kx', 'MarkerSize', 14, 'LineWidth', 2, 'MarkerFaceColor', 'none');
if nargin == 5
    plot(2, mean(lif_co.i_ie)*1e9, 'b+', 'MarkerSize', 14, 'LineWidth', 2, 'MarkerFaceColor', 'none');
    plot(2, mean(lif_co.i_ii)*1e9, 'r+', 'MarkerSize', 14, 'LineWidth', 2, 'MarkerFaceColor', 'none');
    plot(2, (mean(lif_co.i_ii) + mean(lif_co.i_ie))*1e9, 'k+', 'MarkerSize', 14, 'LineWidth', 2, 'MarkerFaceColor', 'none');
end
l4 = plot(2,-mean(nmm_i_ie(round(L/2):end))*nanoamps_scale, 'bo', 'MarkerSize', 10, 'LineWidth', 2, 'MarkerFaceColor', 'none');
l5 = plot(2,-mean(nmm_i_ii(round(L/2):end))*nanoamps_scale, 'ro', 'MarkerSize', 10, 'LineWidth', 2, 'MarkerFaceColor', 'none');
l6 = plot(2,-(mean(nmm_i_ii(round(L/2):end)) + mean(nmm_i_ie(round(L/2):end)))*nanoamps_scale, 'ko', 'MarkerSize', 10, 'LineWidth', 2, 'MarkerFaceColor', 'none');

xlim([0.5 2.5])
% l=legend([l1 l2 l3 l4 l5 l6],{'E_p' 'I_p' '(E + I)_p' 'E_i' 'I_i' '(E + I)_i'});
if nargin == 5, l=legend([l3 l7 l6],{'CUBN' 'COBN' 'NMM'}); else, l=legend([l3 l6],{'LIF' 'NMM'}); end
l.Location = 'eastoutside';
xticks([1 2]);
% xticklabels({'Py_{LIF}' 'In_{LIF}' 'Py_{NMM}' 'In_{NMM}'});
xticklabels({'Pyramidal' 'Interneurons'});
ylabel('input current (nA)');
ax = gca;
ax.FontSize = 12;

if nargin > 3, title(varargin{1}); end
