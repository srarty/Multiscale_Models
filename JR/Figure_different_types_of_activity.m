% Figure with different types of EEG signal (like in wendling2002)
%
% November - 2022

figure;
parameter = 'alpha_ri';
values = [0.1:0.1:0.9 1:5];
N = length(values);
for i = 1:N
%     [x,y,t,f_e,f_i,f_b, p]=NMM_GABAb(parameter, values(i), 'u', 5); model = 'gabab';
%     [x,y,t,f_e,f_i,f_b]=NMM_wendling('u', values(i)); model = 'wendling';
    [x,y,t,f_e,f_i,p] = NMM_diff_equations_DblExp_recursive
    
    subplot(N, 1, i);
    plot(t(1000:end), y(1000:end));
    ylabel([parameter ' = ' num2str(values(i))]);
end
subplot(N,1,1);
title_str = [model ' | u=' num2str(p.u)];
title(title_str);
subplot(N,1,N);
xlabel('Time');


%% Plot a single value
u = 1;
alpha = 5;
syn = 'ri';

[x,y,t,f_e,f_i,f_b]=NMM_GABAb('u', u, ['alpha_' syn], alpha);
figure; plot(t,f_e); hold;
plot(t,f_i)
plot(t,f_b)
legend({'Pyramidal', 'GABA_A', 'GABA_B'}, 'Location', 'best');
xlim([1 3])
title(['\alpha_{' syn '}X=' num2str(alpha) ' | u=' num2str(u) ' |model:gabab'], 'FontWeight', 'normal')
xlabel('Time (s)')
ylabel('Firing rate')
figure;
plot(t(1000:end), y(1000:end));
title(['\alpha_{' syn '}X=' num2str(alpha) ' | u=' num2str(u) ' |model:gabab'], 'FontWeight', 'normal')
xlabel('Time (s)')
ylabel('V_m [mV]')
