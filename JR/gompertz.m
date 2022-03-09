%% GOMPERTZ Computes the nonlinear gompertz sigmoid.
% Input-output function for pyramidal cells. Experimental results (on the
% Brunel) give a different nonlinearity for each neuron type.
%
% Note: Plot the non-linearity with the following code:
%{
        params = set_parameters('allen');       % Chose params.u from a constant value in set_params
        nonlinearity = [];
        count = 0;
        x = -20:0.1:80;
        for i = x
            count = count +1;
            nonlinearity(count) = gompertz(i, params);
        end
        figure
        plot(x, nonlinearity * params.e0, 'LineWidth', 2);
        box off
        grid on
        ylabel('Output firing rate');
        xlabel('Input membrane potential');
        hold;
        maximum = max(nonlinearity * params.e0);
        idx = nonlinearity();
        plot([min(x) max(x)],[0.5 0.5]*maximum,'--k');
        plot([10.7673 10.7673], [0 1]*params.e0,'--k');
%}
% Artemio - February 2022

function out = gompertz(x, params)
    a = params.gompertz.a;
    b = params.gompertz.b;
    c = params.gompertz.c;
    d = params.gompertz.d;
    
    out = a*exp(-b*exp(-d*(x-c)));