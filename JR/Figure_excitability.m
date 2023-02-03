% Varies each alpha parameter and plots the trend of ACFW and t_recovery
% 
% It is only for NMM, to get the LIF figure run Figure_LIF_excitabilty.m

DO_RECOVERY_TIME = true;
DO_ACFW = false;

%% u (external input rate)
disp('Varying u...')
values = 0:0.25:3; 
w_alpha_i = zeros(size(values));
tr_u = zeros(size(values));
fr_in_u = zeros(size(values));
fr_py_u = zeros(size(values));
for i = 1:length(values)
%     [x, y_{i}, t, f_e, f_i] = NMM_diff_equations_DblExp_recursive('alpha_i', values(i));
    [x, y_{i}, t, f_e, f_i] = NMM_GABA('u', values(i));
    if DO_ACFW, w_alpha_i(i) = spectrum(x, y_{i}, t, false); end
    if DO_RECOVERY_TIME, tr_u(i) = analyze_excitability(y_{i},t,1489, -3, 1000, false); end
    fr_in_u(i) = mean(f_i(400:1400));
    fr_py_u(i) = mean(f_e(400:1400));
end
% Save
disp('Saving results for ''u''...');
save 'parameter_sweeps\u.mat' y_ tr_u fr_in_u fr_py_u values
clear y_ tr_u fr_in_u fr_py_u

%% alpha_i
disp('Varying alpha_i...')
values = 0.1:0.1:2; 
w_alpha_i = zeros(size(values));
tr_alpha_i = zeros(size(values));
fr_in_alpha_i = zeros(size(values));
fr_py_alpha_i = zeros(size(values));
for i = 1:length(values)
%     [x, y_{i}, t, f_e, f_i] = NMM_diff_equations_DblExp_recursive('alpha_i', values(i));
    [x, y_{i}, t, f_e, f_i] = NMM_GABA('alpha_i', values(i));
    if DO_ACFW, w_alpha_i(i) = spectrum(x, y_{i}, t, false); end
    if DO_RECOVERY_TIME, tr_alpha_i(i) = analyze_excitability(y_{i},t,1489, -3, 1000, false); end
    fr_in_alpha_i(i) = mean(f_i(400:1400));
    fr_py_alpha_i(i) = mean(f_e(400:1400));
end
% Save
disp('Saving results for ''alpha_i''...');
save 'parameter_sweeps\alpha_i.mat' y_ tr_alpha_i fr_in_alpha_i fr_py_alpha_i values
clear y_ tr_alpha_i fr_in_alpha_i fr_py_alpha_i

%% alpha_e
disp('Varying alpha_e...')
values = 0.1:0.1:2; 
w_alpha_e = zeros(size(values));
tr_alpha_e = zeros(size(values));
fr_in_alpha_e = zeros(size(values));
fr_py_alpha_e = zeros(size(values));
for i = 1:length(values)
%     [x, y_{i}, t, f_e, f_i] = NMM_diff_equations_DblExp_recursive('alpha_e', values(i));
    [x, y_{i}, t, f_e, f_i] = NMM_GABA('alpha_e', values(i));
    if DO_ACFW, w_alpha_e(i) = spectrum(x, y_{i}, t, false); end
    if DO_RECOVERY_TIME, tr_alpha_e(i) = analyze_excitability(y_{i},t,1489, -3, 1000, false); end
    fr_in_alpha_e(i) = mean(f_i(400:1400));
    fr_py_alpha_e(i) = mean(f_e(400:1400));
end
% Save
disp('Saving results for ''alpha_e''...');
save 'parameter_sweeps\alpha_e.mat' y_ tr_alpha_e fr_in_alpha_e fr_py_alpha_e values
clear y_ tr_alpha_e fr_in_alpha_e fr_py_alpha_e

%% alpha_ri
disp('Varying alpha_ri...')
values = 0.1:0.1:2; 
w_alpha_ri = zeros(size(values));
tr_alpha_ri = zeros(size(values));
fr_in_alpha_ri = zeros(size(values));
fr_py_alpha_ri = zeros(size(values));
for i = 1:length(values)
%     [x, y_{i}, t, f_e, f_i] = NMM_diff_equations_DblExp_recursive('alpha_ri', values(i));
    [x, y_{i}, t, f_e, f_i] = NMM_GABA('alpha_ri', values(i));
    if DO_ACFW, w_alpha_ri(i) = spectrum(x, y_{i}, t, false); end
    if DO_RECOVERY_TIME, tr_alpha_ri(i) = analyze_excitability(y_{i},t,1489, -3, 1000, false); end
    fr_in_alpha_ri(i) = mean(f_i(400:1400));
    fr_py_alpha_ri(i) = mean(f_e(400:1400));
end
% Save
disp('Varying alpha_ri...')
disp('Saving results for ''u''...');
save 'parameter_sweeps\alpha_ri.mat' y_ tr_alpha_ri fr_in_alpha_ri fr_py_alpha_ri values
clear y_ tr_alpha_ri fr_in_alpha_ri fr_py_alpha_ri


%% alpha_re
disp('Varying alpha_re...')
values = 0.1:0.1:2; 
w_alpha_re = zeros(size(values));
tr_alpha_re = zeros(size(values));
fr_in_alpha_re = zeros(size(values));
fr_py_alpha_re = zeros(size(values));
for i = 1:length(values)
%     [x, y_{i}, t, f_e, f_i] = NMM_diff_equations_DblExp_recursive('alpha_re', values(i));
    [x, y_{i}, t, f_e, f_i] = NMM_GABA('alpha_re', values(i));
    if DO_ACFW, w_alpha_re(i) = spectrum(x, y_{i}, t, false); end
    if DO_RECOVERY_TIME, tr_alpha_re(i) = analyze_excitability(y_{i},t,1489, -3, 1000, false); end
    fr_in_alpha_re(i) = mean(f_i(400:1400));
    fr_py_alpha_re(i) = mean(f_e(400:1400));
end
% Save
disp('Saving results for ''alpha_re''...');
save 'parameter_sweeps\alpha_re.mat' y_ tr_alpha_re fr_in_alpha_re fr_py_alpha_re values
clear y_ tr_alpha_re fr_in_alpha_re fr_py_alpha_re

%% alpha_rb
disp('Varying alpha_rb...')
values = 0.1:0.1:2; 
w_alpha_b = zeros(size(values));
tr_alpha_b = zeros(size(values));
fr_in_alpha_b = zeros(size(values));
fr_py_alpha_b = zeros(size(values));
y_ = cell(size(values));
for i = 1:length(values)
%     [x, y_{i}, t, f_e, f_i] = NMM_diff_equations_DblExp_recursive('alpha_re', values(i));
    [x, y_{i}, t, f_e, f_i] = NMM_GABA('alpha_b', values(i));
    if DO_ACFW, w_alpha_b(i) = spectrum(x, y_{i}, t, false); end
    if DO_RECOVERY_TIME, tr_alpha_b(i) = analyze_excitability(y_{i},t,1489, -3, 1000, false); end
    fr_in_alpha_b(i) = mean(f_i(400:1400));
    fr_py_alpha_b(i) = mean(f_e(400:1400));
end
% Save
disp('Saving results for ''alpha_rb''...');
save 'parameter_sweeps\alpha_b.mat' y_ tr_alpha_b fr_in_alpha_b fr_py_alpha_b values
clear y_ tr_alpha_b fr_in_alpha_b fr_py_alpha_b

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
    plot(values, w_alpha_b,'-o', 'MarkerSize', 2.5);
    plot([0.1 2], [200 200], '--k')
    box off
    grid
    l = legend({'\alpha_i', '\alpha_e', '\alpha_{ii}', '\alpha_{ee}', '\alpha_{b}'});
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
plot(values, tr_alpha_b,'-o', 'MarkerSize', 2.5);
plot([0.1 2], [30 30], '--k')
box off
grid
l = legend({'\alpha_i', '\alpha_e', '\alpha_{ii}', '\alpha_{ee}', '\alpha_{b}'});
l.Location = 'best';
l.EdgeColor = 'none';
ylabel('t_{recovery} (ms)');
xlabel('Multplier');

%% Load values
load parameter_sweeps/u.mat
load parameter_sweeps/alpha_e.mat
load parameter_sweeps/alpha_i.mat
load parameter_sweeps/alpha_re.mat
load parameter_sweeps/alpha_ri.mat
load parameter_sweeps/alpha_b.mat

%% Recovery time (again on its own figure)
figure
plot(values, tr_alpha_i,'-o', 'MarkerSize', 5);
hold
plot(values, tr_alpha_e,'-o', 'MarkerSize', 5);
plot(values, tr_alpha_ri,'-o', 'MarkerSize', 5);
plot(values, tr_alpha_re,'-o', 'MarkerSize', 5);
plot(values, tr_alpha_b,'-o', 'MarkerSize', 5);
plot([0.1 2], [30 30], '--k')
box off
grid
l = legend({'\alpha_i', '\alpha_e', '\alpha_{ii}', '\alpha_{ee}', '\alpha_{b}'});
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
    plot(values, w_alpha_b .* tr_alpha_b ./ 6e3,'-o', 'MarkerSize', 2.5);
    plot([0.1 2], [1 1], '--k')
    box off
    grid
    l = legend({'\alpha_i', '\alpha_e', '\alpha_{ii}', '\alpha_{ee}', '\alpha_{b}'});
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
    plot(values, weighted_sum(w_alpha_b, tr_alpha_b),'-o', 'MarkerSize', 2.5);
    plot([0.1 2], [1 1], '--k')
    box off
    grid
    l = legend({'\alpha_i', '\alpha_e', '\alpha_{ii}', '\alpha_{ee}', '\alpha_{b}'});
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
plot(values, tr_alpha_b/30,'-o', 'MarkerSize', 2.5);
plot([0 2], [1 1], '--k')
box off
grid
l = legend({'\alpha_i', '\alpha_e', '\alpha_{ii}', '\alpha_{ee}', '\alpha_{b}'});
l.Location = 'best';
l.EdgeColor = 'none';
ylabel('t_{recovery} (a.u.)');
xlabel('Multplier');