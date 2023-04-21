% CALCULATE_BALANCE computes the excitation inhibition balance of the NMM
%
% Input:
%   x       - NMM's state vector.
%   params  - parameters of the NMM.
%
% Output:
%   balance_p - absolute sum of excitatory and inhibitory input currents
%               onto the pyramidal population (nA).
%   balance_i - absolute sum of excitatory and inhibitory input currents
%               onto the inhibitory population (nA).
%
% Artemio - April 2023

function [balance_p, balance_i] = calculate_balance(x, params)

    milivolts_scale = 1e-3; % The results are in mv, this scale is to express them in volts
    nanoamp_scale = 1e9; % To express the output in nA

    nmm_i_pi = (x(:,1) + x(:,9)) * milivolts_scale * params.g_m_P;
    nmm_i_pe = (x(:,5) + x(:,7)) * milivolts_scale * params.g_m_P;
    nmm_i_ie = (x(:,3) + x(:,19)) * milivolts_scale * params.g_m_I;
    nmm_i_ii = (x(:,11)) * milivolts_scale * params.g_m_I;

    L = numel(nmm_i_pe);

    balance_p = -(mean(nmm_i_pi(round(L/2):end)) + mean(nmm_i_pe(round(L/2):end))) * nanoamp_scale;
    balance_i = -(mean(nmm_i_ii(round(L/2):end)) + mean(nmm_i_ie(round(L/2):end))) * nanoamp_scale;

end