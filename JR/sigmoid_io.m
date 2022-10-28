% Error function
% Nonlinearity
function out = sigmoid_io(x, e0, v0, r)
    out = e0 * (0.5 * erf((x - v0) / (sqrt(2)*r)) + 0.5);
end