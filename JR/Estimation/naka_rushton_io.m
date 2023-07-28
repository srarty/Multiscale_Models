% Naka-Rushton function - nonlinearity
function out = naka_rushton_io(x, M, a, b, s)
%     M = 5;
    out = heaviside(x - b) .* (M .* (x - b).^a ./ (s.^a + (x - b).^a));
end
