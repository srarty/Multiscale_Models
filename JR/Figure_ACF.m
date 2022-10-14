% Autocorrelation Figure generator for different parameters
%
%
free_param = 'alpha_i';

legends = {};
values = 0.1:0.1:2;
color = [0 0.4470 0.7410];
for i = 1:length(values)
    [x, y, t, f_e, f_i] = NMM_diff_equations_DblExp_recursive(free_param, values(i));
    x_nmm = y';%(250:end)';
    [auto_nmm, lag_nmm] = autocorr(x_nmm, 'NumLags', 100);

    figure(99)
    if ~ishold, hold on; end
    
    % Compute width of the half section
    w_nmm = 2 * (find(auto_nmm <= 0.5, 1) - 1); % -1 because the first sample is at Lag=0

    plty = [fliplr(auto_nmm) auto_nmm];
    pltx = [fliplr(-lag_nmm) lag_nmm];
    l_nmm(i) = plot(pltx', plty', 'Color', values(i) * color * 0.5);
    % plot([-w_nmm/2 w_nmm/2], [0.5 0.5], '--', 'Color', l_nmm.Color);
    legends{i,1} = [num2str(values(i)) '\' free_param];
end

xlabel('Lag (@10^{-3}samples/s)');
ylabel('Autocorrelation');
grid on
box off
legend(legends);
