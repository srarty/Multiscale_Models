%% Runs u vs gain and saves the values
% Run this code before Figure_colormaps
%
% Artemio - May 2023

var_vec = {'alpha_i', 'alpha_e', 'alpha_ri', 'alpha_re', 'alpha_b'};
% var_vec = {'alpha_b'};

for j = 1:length(var_vec)

    value = var_vec{j}; %'alpha_i';
    range = 0.4:0.1:4;%0:0.1:2;
    value2 = 'u';
    range2 = 0:0.1:2;%0.4:0.1:4;%

    freqs = [];
    freqs_py = [];
    freqs_in = [];
    w_nmm = [];
    w_lif = [];
    recovery = [];%size(range);
    for i = 1:length(range)
        for ii = 1:length(range2)
    %         [x, y, t] = NMM_diff_equations_DblExp_recursive(value, range(i));
            [x, y, t, f_e, f_i] = NMM_GABA(value, range(i), value2, range2(ii));

    %         freqs(ii,i) = spectrum(x,y,t, false); % Oscillations
    %         try
    %             [w_nmm(ii,i), w_lif(ii,i)] = spectrum(x,y,t, false); % Width of Autocorrelation function
    %         catch E
    %             if strcmp('MATLAB:subsassigndimmismatch', E.identifier)
    %                 w_nmm = 2000;
    %             else
    %                 rethrow(E);
    %             end
    %         end


            freqs_py(ii,i) = mean(f_e(250:end)); % spike rate Py
            freqs_in(ii,i) = mean(f_i(250:end)); % spike rate Py

            disp([num2str(i) '/' num2str(length(range)) ' , ' num2str(ii) '/' num2str(length(range2))]);

    %         recovery(i,ii) = analyze_excitability(y,t);
        end
    end


    %% Store values

    results = struct;
    results.value = value;
    results.range = range;
    results.value2 = value2;
    results.range2 = range2;
    results.freqs = freqs;
    results.freqs_py = freqs_py;
    results.freqs_in = freqs_in;

    folder = 'C:\Users\artemios\Documents\Multiscale_Models_Data\Spartan\firing_rates_nmm3\';
    name = [value ' vs ' value2];
    if isempty(dir([folder name]))
        save([folder name '.mat'], 'results');
        disp('results saved');
    else
        error(['Results not saved, file exists | j = ' num2str(j)]);
    end
    
end
