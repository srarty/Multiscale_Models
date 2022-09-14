% Generates figure for comparison between alpha's and j's

%% LIF parameters (pA)
j_pi = 37;
j_pp = -73.5;
j_pu = -11.275; % -1.375 * 8.2 = -11.275 (check lif_model.py)
j_ip = -165;
j_ii = 17.55;

%% Alpha parameters
p = set_parameters('recursive');

ratio_pi = j_pi / p.alpha_i;
ratio_pp = j_pp / p.alpha_re;
ratio_pu = j_pu / p.alpha_u;
ratio_ip = j_ip / p.alpha_e;
ratio_ii = j_ii / p.alpha_ri;

%% Plotting
f = @(a,ratio) a*ratio;

figure

% P -> P
alphas = [0:5];
js = f([alphas], ratio_pp);
subplot(2,2,1);
plot(js, alphas, '-ok', 'MarkerFaceColor', 'k', 'LineWidth', 1);
hold
plot([j_pp j_pp], [min(alphas) max(alphas)], '--k');
plot([min(js) max(js)], [p.alpha_re p.alpha_re], '--k');
ylim([min(alphas) max(alphas)]);
xlim([min(js) max(js)]);
box off
grid on
ylabel('\alpha_{pp}');
xlabel('j_{P_{AMPA}}');

% I -> P
alphas = -1*[0:5];
js = f([alphas], ratio_pi);
subplot(2,2,2);
plot(js, alphas, '-ok', 'MarkerFaceColor', 'k', 'LineWidth', 1);
hold
plot([j_pi j_pi], [min(alphas) max(alphas)], '--k');
plot([min(js) max(js)], [p.alpha_i p.alpha_i], '--k');
ylim([min(alphas) max(alphas)]);
xlim([min(js) max(js)]);
box off
grid on
ylabel('\alpha_{pi}');
xlabel('j_{P_{GABA}}');

% P -> I
alphas = [0:5];
js = f([alphas], ratio_ip);
subplot(2,2,3);
plot(js, alphas, '-ok', 'MarkerFaceColor', 'k', 'LineWidth', 1);
hold
plot([j_ip j_ip], [min(alphas) max(alphas)], '--k');
plot([min(js) max(js)], [p.alpha_e p.alpha_e], '--k');
ylim([min(alphas) max(alphas)]);
xlim([min(js) max(js)]);
box off
grid on
ylabel('\alpha_{ip}');
xlabel('j_{I_{AMPA}}');

% I -> I
alphas = -1*[0:5];
js = f([alphas], ratio_ii);
subplot(2,2,4);
plot(js, alphas, '-ok', 'MarkerFaceColor', 'k', 'LineWidth', 1);
hold
plot([j_ii j_ii], [min(alphas) max(alphas)], '--k');
plot([min(js) max(js)], [p.alpha_ri p.alpha_ri], '--k');
ylim([min(alphas) max(alphas)]);
xlim([min(js) max(js)]);
box off
grid on
ylabel('\alpha_{ii}');
xlabel('j_{I_{GABA}}');