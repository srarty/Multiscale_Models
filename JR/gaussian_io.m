% Gaussian function - Nonlinearity
function out = gaussian_io(x, a, b, c, d)
    out = a*exp(-((x-b)/c).^2);
end