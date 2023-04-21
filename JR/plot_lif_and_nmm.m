% Plots a comparison of LFP and Firing rates of LIF and NMM. Also plots LFP histogram.
%
% Called by Figure_sample_recordings.m

function plot_lif_and_nmm(lif, t, f_e, f_i, yy, params, title_str)

    yy = yy(200:end) - yy(1); % Adjust to start at zero like the LIF
    t = t(200:end);
    f_e = f_e(200:end);
    f_i = f_i(200:end);
    lif.i_pe = lif.i_pe(2000:end);
    lif.i_pi = lif.i_pi(2000:end);
    lif.R_py = lif.R_py(2000:end);
    lif.R_in = lif.R_in(2000:end);

    figure; subplot(211); plot(t, yy, 'k');
    hold; plot(linspace(t(1),t(end),length(lif.i_pe)), -(lif.i_pe - lif.i_pi)/params.g_m_P, '--k')
    xlim([0.2 2]);
    ylabel('LFP (V)');
    title(title_str);
    l = legend({'NMM', 'LIF'});
    l.Location = 'best';
    subplot(212); l1=plot(t, f_e); hold; l2=plot(t, f_i);
    plot(linspace(t(1),t(end),length(lif.R_py)), lif.R_py, '--', 'Color', l1.Color); plot(linspace(t(1),t(end),length(lif.R_in)), lif.R_in, '--', 'Color', l2.Color);
    xlim([0.2 2]);
    xlabel('Time (s)');
    ylabel('Firing rates (Hz)');

    figure
    histogram(downsample(-(lif.i_pe - lif.i_pi)/params.g_m_P,10));%, 'BinWidth', 0.01);
    hold
    histogram(yy);%, 'BinWidth', 0.01);
    title(title_str);

    drawnow

end