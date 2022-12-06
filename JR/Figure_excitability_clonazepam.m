% Generates figure:  Excitability vs GABA_A agonist
%
%   + barbiturates (e.g. ??) - increasing the duration of chloride ion
%       channel opening at the GABA_A receptor, i.e.increases the time 
%       constant of GABA_A inhibitory synapses.
%       NOTE: increases the toxicity in overdose because it directly
%       controls the opening of GABA_A Cloride ion channel.
%   + benzodiazepines (e.g. midazolam, diazepam) -  increase the frequency 
%       of the chloride ion channel opening at the GABAA receptor, i.e. 
%       increases the potency of GABA). 
%       NOTE: indirect action on the cloride channel.
%       REFS: Rovira et al. 1993 (Fig 2); 
%
% Beta activity increase: Nunez, page 435. , Niedermeyers and Lopes da
% Silva 1999, 

values = 0.5:0.1:2.5;

%% Recovery time (Response to pulse)
recovery_times = zeros(3, length(values));
for i = 1:length(values)
    [x, y, t]=NMM_GABAb('alpha_i', values(i), 'alpha_ri', values(i), 'u', 1); % Benzodiazepine
%     [x, y, t]=NMM_GABAb('tau_sp', values(i), 'tau_sri', values(i), 'u', 0); % Benzodiazepine
%     [x, y, t]=NMM_GABAb('tau_sri', values(i), 'u', 0); % Benzodiazepine
    recovery_times(1, i) = analyze_excitability(y,t, 1489, -3, 1000, false);
    
    [~, y, t]=NMM_GABAb('alpha_i', values(i), 'u', 1); % Inhibitory -> Pyramidal
    recovery_times(2, i) = analyze_excitability(y,t, 1489, -3, 1000, false);
        
    [~, y, t]=NMM_GABAb('tau_sp',values(i), 'tau_sri', values(i), 'alpha_i', values(i), 'alpha_ri', values(i), 'u', 1);% Barbiturate 
    recovery_times(3, i) = analyze_excitability(y,t, 1489, -3, 1000, false);
end    
    
%% ACFW (Background activity)
w = zeros(3, length(values));
for i = 1:length(values)    
    [x, y, t]=NMM_GABAb('alpha_i', values(i), 'alpha_ri', values(i), 'CURRENT', 0, 'u', 1.5); % Benzodiazepine
    [auto_nmm, lag_nmm] = autocorr(y(1000:end), 'NumLags', 100);
    w(1,i) = 2 * (find(auto_nmm <= 0.5, 1) - 1);
    
    [~, y, t]=NMM_GABAb('alpha_i', values(i), 'CURRENT', 0, 'u', 1.5); % Inhibitory -> Pyramidal
    [auto_nmm, lag_nmm] = autocorr(y(1000:end), 'NumLags', 100);
    w(2,i) = 2 * (find(auto_nmm <= 0.5, 1) - 1);
    
    [~, y, t]=NMM_GABAb('tau_sp',values(i), 'tau_sri', values(i), 'alpha_i', values(i), 'alpha_ri', values(i), 'CURRENT', 0, 'u', 1.5);% Barbiturate    
    [auto_nmm, lag_nmm] = autocorr(y(1000:end), 'NumLags', 100);
    w(3,i) = 2 * (find(auto_nmm <= 0.5, 1) - 1);
end

%% Plot Recovery time
figure; 
plot(values, recovery_times(1,:), 'o', 'LineWidth', 2);
hold
plot(values, recovery_times(2,:), 'x', 'LineWidth', 2);
plot(values, recovery_times(3,:), '^', 'LineWidth', 2);
legend({'Benzodiazepine', '\alpha_{pi}', 'Barbiturate'});
xlabel('Concentration (a.u.)')
ylabel('t_{recovery} (ms)')
% xlim([1 2.5]);
title('t_r | u=0')

%% Plot ACFW
figure; 
plot(values, w(1,:), 'o', 'LineWidth', 2);
hold
plot(values, w(2,:), 'x', 'LineWidth', 2);
plot(values, w(3,:), '^', 'LineWidth', 2);
legend({'\alpha_{GABA_{A}}', '\alpha_{pi}', '\tau_{GABA_{A}}'});
title('ACFW | u=1.5')

%% Plot autocorrelation function for single runs:
%{
figure
plty = [fliplr(auto_nmm) auto_nmm];
pltx = [fliplr(-lag_nmm) lag_nmm];
l = plot(pltx, plty, 'k');
%}