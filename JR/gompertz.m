%% GOMPERTZ Computes the nonlinear gompertz sigmoid.
% Input-output function for pyramidal cells. Experimental results (on the
% Brunel) give a different nonlinearity for each neuron type.
%
% Note: Plot the non-linearity with the following code:
%{
        params = set_parameters('recursive');       % Chose params.u from a constant value in set_params
        nonlinearity = [];
        count = 0;
        x = -20:0.1:80;
        for i = x
            count = count +1;
            nonlinearity(count) = gompertz(i, params.e0i, params.gompertzi.b, params.gompertzi.c, params.gompertzi.d);
        end
        figure
        %plot(x, nonlinearity, 'LineWidth', 2, 'color', [0.0745 0.6235 1]);
        plot(x, nonlinearity, 'LineWidth', 2, 'color', [1 0.4118 0.1608]);
        box off
        grid on
        ylabel('Output firing rate (spikes/s)');
        xlabel('Input membrane potential (mV)');
        hold;
        maximum = max(nonlinearity);
        idx = nonlinearity();
        plot([min(x) max(x)],[0.5 0.5]*maximum,'--k');
        plot([7.2 7.2], [0 1]*maximum,'--k');
        title('Interneurons population');
%}
% Artemio - February 2022

function out = gompertz(x, a, b, c, d)    
    out = a * exp(-b*exp(-d*(x-c)));