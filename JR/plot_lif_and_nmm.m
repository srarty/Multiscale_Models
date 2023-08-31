% Plots a comparison of LFP and Firing rates of LIF and NMM. Also plots LFP histogram.
%
% Called by Figure_sample_recordings.m

function plot_lif_and_nmm(lif, t, f_e, f_i, yy, x, params, title_str, varargin)
    
    yy = yy(500:1000) - yy(1); % Adjust to start at zero like the LIF
    t = t(500:1000);
    f_e = f_e(500:1000);
    f_i = f_i(500:1000);
    lif.i_pe = lif.i_pe(5000:10000);
    lif.i_pi = lif.i_pi(5000:10000);
    lif.R_py = lif.R_py(5000:10000);
    lif.R_in = lif.R_in(5000:10000);
    if nargin == 9
        lif_co = varargin{1};
        lif_co.i_pe = lif_co.i_pe(5000:10000);
        lif_co.i_pi = lif_co.i_pi(5000:10000);
        lif_co.R_py = lif_co.R_py(5000:10000);
        lif_co.R_in = lif_co.R_in(5000:10000);
    end
    
    f = figure; 
    f.Position = [657 252 560 680];
    subplot(211); 
    hold; 
    plot(linspace(t(1),t(end),length(lif.i_pe)), (-(lif.i_pe - lif.i_pi)/params.g_m_P)*1e3, '-', 'Color', [0.8 0.8 0.8]); % LIF
    if nargin == 9, plot(linspace(t(1),t(end),length(lif_co.i_pe)), (-(lif_co.i_pe - lif_co.i_pi)/params.g_m_P)*1e3, '-', 'Color', [0.5 0.5 0.5]); end% LIF
    plot(t, yy*1e3, 'k'); % NMM
    xlim([0.5 0.6]);
    ylabel('LFP (mV)');
    title(title_str);
    if nargin == 9, l = legend({'CUBN', 'COBN', 'NMM'}); else, l = legend({'LIF', 'NMM'});end
    l.Location = 'eastoutside';
    ax = gca;
    ax.FontSize = 12;
    
    subplot(212); l1=plot(t, f_e); hold; l2=plot(t, f_i);
    l3 = plot(linspace(t(1),t(end),length(lif.R_py)), lif.R_py, '--', 'Color', l1.Color); l4 = plot(linspace(t(1),t(end),length(lif.R_in)), lif.R_in, '--', 'Color', l2.Color);
    if nargin == 9, l5 = plot(linspace(t(1),t(end),length(lif_co.R_py)), lif_co.R_py, '-.', 'Color', l1.Color); l6 = plot(linspace(t(1),t(end),length(lif_co.R_in)), lif_co.R_in, '-.', 'Color', l2.Color); end
    xlim([0.5 0.6]);
    xlabel('Time (s)');
    ylabel('Firing rates (Hz)');
    if nargin == 9
        l = legend({'NMM_{Py}', 'NMM_{In}', 'CUBN_{Py}', 'CUBN_{In}', 'COBN_{Py}', 'COBN_{In}'});
    else
        l = legend({'Py_{NMM}', 'In_{NMM}', 'Py_{LIF}', 'In_{LIF}'});
    end
    l.Location = 'eastoutside';
    ax = gca;
    ax.FontSize = 12;
    box off;
    
    f = figure;
    f.Position = [1219 253 420 300];
    yy_lif = downsample(-(lif.i_pe - lif.i_pi)/params.g_m_P,10);
    h = histogram(yy_lif*1e3, 'FaceColor', [0.8 0.8 0.8]);%, 'BinWidth', 0.01);
    hold
    if nargin == 9
        yy_lif_co = downsample(-(lif_co.i_pe - lif_co.i_pi)/params.g_m_P,10);
        histogram(yy_lif_co*1e3, 'FaceColor', [0.5 0.5 0.5]); 
    end
    histogram(yy*1e3, 'FaceColor', [0 0 0], 'BinWidth', h.BinWidth);
    if nargin == 9
        legend({['CUBN (k = ' num2str(kurtosis(yy_lif,0)) ')'] ['COBN (k = ' num2str(kurtosis(yy_lif_co,0)) ')'] ['NMM (k = ' num2str(kurtosis(yy,0)) ')'] });
    else
        legend({['LIF (k = ' num2str(kurtosis(yy_lif,0)) ')'] ['NMM (k = ' num2str(kurtosis(yy,0)) ')'] });
    end
    ax = gca;
    ax.FontSize = 12;
    xlabel('milivolts');
    ylabel('Distribution');
    
    
    %Balance
    % To plot excitation-inhibition balance, run Figure_balance.m 
    if nargin == 9
        Figure_balance(x, params, lif, title_str, lif_co);
    else
        Figure_balance(x, params, lif, title_str);
    end
    f = gcf;
    f.Position = [1219 645 420 280];

    drawnow

end