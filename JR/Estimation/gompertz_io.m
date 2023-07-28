% Gompertz function - nonlinearity
function out = gompertz_io(x, a, b, c, d)
%     a = 5;
    out = a * exp(-b*exp(-d*(x-c)));
end

%% Plot the Gompertz function:
%{
    p = set_parameters('default');

    x = -20:0.1:80;
    y = gompertz_io(x, p.e0, p.gompertz.b, p.gompertz.c, p.gompertz.d);
    th_idx = find(y>=0.5*p.e0, 1);

    figure
    plot(x, y, 'LineWidth', 2);
    hold;
    plot([min(x) max(x)],[y(th_idx) y(th_idx)],'--k');
    plot([x(th_idx) x(th_idx)], [0 30], '--k'); %[0 1]*p.e0,'--k');
    xlim([-20 40]);
    ylim([0 30]);

    box off
    grid on
    ylabel('Output firing rate (spikes/s)');
    xlabel('Input membrane potential (mV)');
    title('Pyramidal');
%}

%{
    p = set_parameters('default');

    x = -20:0.1:80;
    y = gompertz_io(x, p.e0i, p.gompertzi.b, p.gompertzi.c, p.gompertzi.d);
    th_idx = find(y>=0.5*p.e0i, 1);

    figure
    plot(x, y, 'LineWidth', 2, 'color', [1 0.4118 0.1608]);
    hold;
    plot([min(x) max(x)],[y(th_idx) y(th_idx)],'--k');
    plot([x(th_idx) x(th_idx)], [0 60], '--k'); %[0 1]*p.e0,'--k');
    xlim([-20 40]);
    ylim([0 60]);

    box off
    grid on
    ylabel('Output firing rate (spikes/s)');
    xlabel('Input membrane potential (mV)');
    title('Interneurons');
%}