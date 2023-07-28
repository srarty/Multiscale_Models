% Plots a comparison of LFP and Firing rates of LIF and NMM. Also plots LFP histogram.
%
% Called by Figure_sample_recordings.m

function plot_lif_and_nmm(lif, t, f_e, f_i, yy, x, params, title_str)

    yy = yy(500:1000) - yy(1); % Adjust to start at zero like the LIF
    t = t(500:1000);
    f_e = f_e(500:1000);
    f_i = f_i(500:1000);
    lif.i_pe = lif.i_pe(5000:10000);
    lif.i_pi = lif.i_pi(5000:10000);
    lif.R_py = lif.R_py(5000:10000);
    lif.R_in = lif.R_in(5000:10000);

    f = figure; 
    f.Position = [657 252 560 680];
    subplot(211); 
    hold; 
    plot(linspace(t(1),t(end),length(lif.i_pe)), (-(lif.i_pe - lif.i_pi)/params.g_m_P)*1e3, '--', 'Color', [0.5 0.5 0.5]); % LIF
    plot(t, yy*1e3, 'k'); % NMM
    xlim([0.5 1]);
    ylabel('LFP (mV)');
    title(title_str);
    l = legend({'LIF', 'NMM'});
    l.Location = 'best';
    ax = gca;
    ax.FontSize = 12;
    
    subplot(212); l1=plot(t, f_e); hold; l2=plot(t, f_i);
    plot(linspace(t(1),t(end),length(lif.R_py)), lif.R_py, '--', 'Color', l1.Color); plot(linspace(t(1),t(end),length(lif.R_in)), lif.R_in, '--', 'Color', l2.Color);
    xlim([0.5 1]);
    xlabel('Time (s)');
    ylabel('Firing rates (Hz)');
    l = legend({'Py_{NMM}', 'In_{NMM}', 'Py_{LIF}', 'In_{LIF}'});
    l.Location = 'best';
    ax = gca;
    ax.FontSize = 12;
    
    f = figure;
    f.Position = [1219 253 420 300];
    % subplot(313); 
    yy_lif = downsample(-(lif.i_pe - lif.i_pi)/params.g_m_P,10);
    h = histogram(yy_lif*1e3, 'FaceColor', [0.5 0.5 0.5]);%, 'BinWidth', 0.01);
    hold
    histogram(yy*1e3, 'FaceColor', [0 0 0], 'BinWidth', h.BinWidth);
%     title(title_str);
    legend({['LIF (k = ' num2str(kurtosis(yy_lif,0)) ')'] ['NMM (k = ' num2str(kurtosis(yy,0)) ')'] });
    ax = gca;
    ax.FontSize = 12;
    xlabel('milivolts');
    ylabel('Distribution');
    
    
    %Balance
    % To plot excitation-inhibition balance, run Figure_balance.m 
    Figure_balance(x, params, lif, title_str);
    f = gcf;
    f.Position = [1219 645 420 280];

    drawnow

end