% Varies each alpha parameter and plots the trend of balance, t_recovery and ACFW
% 
% It is only for NMM, to get the LIF figure run Figure_LIF_excitabilty.m

DO_RECOVERY_TIME = true;
DO_BALANCE = true;
DO_STATE = true;
DO_ACFW = false;

f1 = figure();
f2 = figure();

% var_vec = {'alpha_e', 'alpha_i', 'alpha_re', 'alpha_ri', 'alpha_b'};
var_vec = {'alpha_e', 'alpha_i', 'alpha_b'};
% var_vec = {'alpha_re', 'alpha_ri'};
% var_vec = {'alpha_e'};
values = 0.0:0.1:2;
%%

for j = 1:length(var_vec)
    current_var = var_vec{j};
    
    disp(['Varying ' current_var '...']);

    w = zeros(length(var_vec), length(values));
    tr = zeros(length(var_vec), length(values));
    fr_in = zeros(length(var_vec), length(values));
    fr_py = zeros(length(var_vec), length(values));
    bal_p = zeros(length(var_vec), length(values));
    bal_i = zeros(length(var_vec), length(values));
    state = zeros(length(var_vec), length(values));
    value_ = zeros(length(var_vec), length(values));

    for i = 1:length(values)
    %     [x, y_{i}, t, f_e, f_i] = NMM_diff_equations_DblExp_recursive('alpha_i', values(i));
    %     [x, y_{i}, t, f_e, f_i] = NMM_GABA('alpha_i', values(i), 'u', 0);
    
        % Run the model
        if DO_RECOVERY_TIME, [~, ~, t_excitability, ~, ~, ~, y_excitability{j,i}] = NMM_GABA(current_var, values(i), 'u', 0, 'CURRENT', 50e-12); end
        [x, ~, t, f_e, f_i, params, y_{j,i}] = NMM_GABA(current_var, values(i), 'u', 0);
        
        % Run the tests:
        if DO_ACFW, w(j,i) = spectrum(x, y_{j,i}, t, false); end %#ok<UNRCH>
        if DO_RECOVERY_TIME, tr(j,i) = analyze_excitability(y_excitability{j,i}, t_excitability, 1489, -3, 1000, false); end
        if DO_BALANCE, [bal_p(j,i), bal_i(j,i)] = calculate_balance(x, params); end
        
        if DO_STATE
            % Calculate fft to estimate oscillatory activity
            [~, X_] = fft_plot(y_{j,i}(1000:end)-mean(y_{j,i}(1000:end)), t(1000:end),[],false);
            if (mean(f_i(500:end)) > 60) || (mean(f_e(500:end)) > 35)
                state(j,i) = 2; % Saturation
            elseif max(X_) > 1e8
                state(j,i) = 1; % Oscillation
            else
                state(j,i) = 0; % Normal
            end
        end
        
        % Save firing rate
        fr_in(j,i) = mean(f_i(400:1400));
        fr_py(j,i) = mean(f_e(400:1400));
        
        % Store the current value
        value_(j,i) = values(i);
    end
    
    figure(f1)
%     subplot(211)
    if ~ishold(gca), hold; end
    l = plot(value_(j,:), bal_p(j,:), '-^');
    plot(value_(j,:), bal_i(j,:), '-x', 'Color', l.Color);
    xlabel('Gain factor');
    ylabel('Balance (nA)');

    
%     subplot(212)
    figure(f2)
    if ~ishold(gca), hold; end
    plot(value_(j,:), 1e3 * tr(j,:), '-o', 'Color', l.Color);
    ylabel('t_{r} (ms)');
    xlabel('Gain factor');
    legend(var_vec)

end
return


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