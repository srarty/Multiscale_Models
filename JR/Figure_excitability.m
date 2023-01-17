% Varies each alpha parameter and plots the trend of ACFW and t_recovery
%

DO_RECOVERY_TIME = true;
DO_ACFW = false;

%% alpha_i
values = 0.1:0.1:2; 
w_alpha_i = zeros(size(values));
tr_alpha_i = zeros(size(values));
for i = 1:length(values)
    [x, y, t, f_e, f_i] = NMM_diff_equations_DblExp_recursive('alpha_i', values(i));
    if DO_ACFW, w_alpha_i(i) = spectrum(x, y, t, false); end
    if DO_RECOVERY_TIME, tr_alpha_i(i) = analyze_excitability(y,t,1489, -3, 1000); end
end

%% alpha_e
values = 0.1:0.1:2; 
w_alpha_e = zeros(size(values));
tr_alpha_e = zeros(size(values));
for i = 1:length(values)
    [x, y, t, f_e, f_i] = NMM_diff_equations_DblExp_recursive('alpha_e', values(i));
    if DO_ACFW, w_alpha_e(i) = spectrum(x, y, t, false); end
    if DO_RECOVERY_TIME, tr_alpha_e(i) = analyze_excitability(y,t,1489, -3, 1000); end
end

%% alpha_ri
values = 0.1:0.1:2; 
w_alpha_ri = zeros(size(values));
tr_alpha_ri = zeros(size(values));
for i = 1:length(values)
    [x, y, t, f_e, f_i] = NMM_diff_equations_DblExp_recursive('alpha_ri', values(i));
    if DO_ACFW, w_alpha_ri(i) = spectrum(x, y, t, false); end
    if DO_RECOVERY_TIME, tr_alpha_ri(i) = analyze_excitability(y,t,1489, -3, 1000); end
end

%% alpha_re
values = 0.1:0.1:2; 
w_alpha_re = zeros(size(values));
tr_alpha_re = zeros(size(values));
for i = 1:length(values)
    [x, y, t, f_e, f_i] = NMM_diff_equations_DblExp_recursive('alpha_re', values(i));
    if DO_ACFW, w_alpha_re(i) = spectrum(x, y, t, false); end
    if DO_RECOVERY_TIME, tr_alpha_re(i) = analyze_excitability(y,t,1489, -3, 1000); end
end

%% Plot
% ACFW
figure
if DO_ACFW
    subplot(2,1,1)
    plot(values, w_alpha_i,'-o', 'MarkerSize', 2.5);
    hold
    plot(values, w_alpha_e,'-o', 'MarkerSize', 2.5);
    plot(values, w_alpha_ri,'-o', 'MarkerSize', 2.5);
    plot(values, w_alpha_re,'-o', 'MarkerSize', 2.5);
    plot([0.1 2], [200 200], '--k')
    box off
    grid
    l = legend({'\alpha_i', '\alpha_e', '\alpha_{ii}', '\alpha_{ee}'});
    l.Location = 'best';
    l.EdgeColor = 'none';
    ylabel('ACFW (samples)');
end

%% Recovery time
subplot(2,1,2)
plot(values, tr_alpha_i,'-o', 'MarkerSize', 2.5);
hold
plot(values, tr_alpha_e,'-o', 'MarkerSize', 2.5);
plot(values, tr_alpha_ri,'-o', 'MarkerSize', 2.5);
plot(values, tr_alpha_re,'-o', 'MarkerSize', 2.5);
plot([0.1 2], [30 30], '--k')
box off
grid
l = legend({'\alpha_i', '\alpha_e', '\alpha_{ii}', '\alpha_{ee}'});
l.Location = 'best';
l.EdgeColor = 'none';
ylabel('t_{recovery} (ms)');
xlabel('Multplier');

%% Recovery time (again on its own figure)
figure
plot(values, tr_alpha_i,'-o', 'MarkerSize', 5);
hold
plot(values, tr_alpha_e,'-o', 'MarkerSize', 5);
plot(values, tr_alpha_ri,'-o', 'MarkerSize', 5);
plot(values, tr_alpha_re,'-o', 'MarkerSize', 5);
plot([0.1 2], [30 30], '--k')
box off
grid
l = legend({'\alpha_i', '\alpha_e', '\alpha_{ii}', '\alpha_{ee}'});
l.Location = 'northeastoutside';
l.EdgeColor = 'black';
ylabel('t_{recovery} (ms)');
xlabel('Multplier');

%% ACFW * t_recovery / 6000
if DO_ACFW
    figure
    plot(values, w_alpha_i .* tr_alpha_i ./ 6e3,'-o', 'MarkerSize', 2.5);
    hold
    plot(values, w_alpha_e .* tr_alpha_e ./ 6e3,'-o', 'MarkerSize', 2.5);
    plot(values, w_alpha_ri .* tr_alpha_ri ./ 6e3,'-o', 'MarkerSize', 2.5);
    plot(values, w_alpha_re .* tr_alpha_re ./ 6e3,'-o', 'MarkerSize', 2.5);
    plot([0.1 2], [1 1], '--k')
    box off
    grid
    l = legend({'\alpha_i', '\alpha_e', '\alpha_{ii}', '\alpha_{ee}'});
    l.Location = 'best';
    l.EdgeColor = 'none';
    ylabel('Excitability (a.u.)');
    xlabel('Multplier');
    title('Product')


    %% Weighted sum:
    % ( 2*(t_recovery/30) + ACFW/200 )/3
    weighted_sum = @(w, tr) ( 2*(tr/30) + w/200 )/3; % Evaluation function

    figure
    plot(values, weighted_sum(w_alpha_i, tr_alpha_i),'-o', 'MarkerSize', 2.5);
    hold
    plot(values, weighted_sum(w_alpha_e, tr_alpha_e),'-o', 'MarkerSize', 2.5);
    plot(values, weighted_sum(w_alpha_ri, tr_alpha_ri),'-o', 'MarkerSize', 2.5);
    plot(values, weighted_sum(w_alpha_re, tr_alpha_re),'-o', 'MarkerSize', 2.5);
    plot([0.1 2], [1 1], '--k')
    box off
    grid
    l = legend({'\alpha_i', '\alpha_e', '\alpha_{ii}', '\alpha_{ee}'});
    l.Location = 'best';
    l.EdgeColor = 'none';
    ylabel('Excitability (a.u.)');
    xlabel('Multplier');
    title('Weighted sum')
end


%% Normalized recovery time
figure
plot(values, tr_alpha_i/30,'-o', 'MarkerSize', 2.5);
hold
plot(values, tr_alpha_e/30,'-o', 'MarkerSize', 2.5);
plot(values, tr_alpha_ri/30,'-o', 'MarkerSize', 2.5);
plot(values, tr_alpha_re/30,'-o', 'MarkerSize', 2.5);
plot([0 2], [1 1], '--k')
box off
grid
l = legend({'\alpha_i', '\alpha_e', '\alpha_{ii}', '\alpha_{ee}'});
l.Location = 'best';
l.EdgeColor = 'none';
ylabel('t_{recovery} (a.u.)');
xlabel('Multplier');