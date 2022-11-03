% Generates figure:  Excitability vs GABA_A agonist
%
% TODO: Vary the gains to see trend

values = 0.5:0.1:2.5;
recovery_times = zeros(3, length(values));
w = zeros(3, length(values));

for i = 1:length(values)
    
    [x, y, t]=NMM_GABAb('alpha_i', values(i), 'alpha_ri', values(i)); % Benzodiazepine
    recovery_times(1, i) = analyze_excitability(y,t, 1489, 0.05, 1000, false);
    
    [auto_nmm, lag_nmm] = autocorr(y, 'NumLags', 100);
    w(1,i) = 2 * (find(auto_nmm <= 0.5, 1) - 1);
    
    [~, y, t]=NMM_GABAb('alpha_i', values(i)); % Inhibitory -> Pyramidal
    recovery_times(2, i) = analyze_excitability(y,t, 1489, 0.05, 1000, false);
    [auto_nmm, lag_nmm] = autocorr(y, 'NumLags', 100);
    w(2,i) = 2 * (find(auto_nmm <= 0.5, 1) - 1);
    
    [~, y, t]=NMM_GABAb('tau_sp',values(i), 'tau_sri', values(i), 'alpha_i', values(i), 'alpha_ri', values(i));% Barbiturate 
    recovery_times(3, i) = analyze_excitability(y,t, 1489, 0.05, 1000, false);
    [auto_nmm, lag_nmm] = autocorr(y, 'NumLags', 100);
    w(3,i) = 2 * (find(auto_nmm <= 0.5, 1) - 1);
end

%% Plot
% Recovery time
figure; 
plot(values, recovery_times(1,:), 'o');
hold
plot(values, recovery_times(2,:), 'x');
plot(values, recovery_times(3,:), '^');
legend({'\alpha_{GABA_{A}}', '\alpha_{pi}', '\tau_{GABA_{A}}'});

% ACFW
figure; 
plot(values, w(1,:), 'o');
hold
plot(values, w(2,:), 'x');
plot(values, w(3,:), '^');
legend({'\alpha_{GABA_{A}}', '\alpha_{pi}', '\tau_{GABA_{A}}'});
