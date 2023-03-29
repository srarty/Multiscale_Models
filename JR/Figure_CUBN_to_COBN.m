% Generates the supplementary material figure with the results of the CUBN
% to COBN conversion

%% Load the results
% load('C:\Users\artemios\Dropbox\University of Melbourne\Epilepsy\Multiscale_Models\JR\parameter_sweeps\CUBNtoCOBNresults.mat');
load('C:\Users\artemios\Dropbox\University of Melbourne\Epilepsy\Multiscale_Models\JR\parameter_sweeps\CUBNtoCOBNresults2.mat');

%% Plot
figure
plot(Py_diff*1e3);
hold
plot(In_diff*1e3);
ylim([-1 10]);
ylim([-0.25 10]);
box off
plot((Py_diff + In_diff)*1e3, 'LineWidth', 2);
plot([0 numel(k_AMPA_I)],[0.01 0.01],'--k');
legend({'Error_{P}' 'Error_{I}' 'Error_{P} + Error_{I}'});
xlabel('iteration');
ylabel('Error (mV)');

figure
plot(k_GABA_P);
hold
plot(k_GABAb_P);
plot(k_AMPA_P);
plot(k_GABA_I);
plot(k_AMPA_I);
xlabel('iteration');
ylabel('k_{syn_a}');
box off
grid
legend({'k_{GABA_A_P}' 'k_{GABA_B_P}' 'k_{AMPA_P}' 'k_{GABA_A_I}' 'k_{AMPA_I}'});