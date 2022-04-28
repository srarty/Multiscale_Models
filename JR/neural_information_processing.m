K = 5;
f = @(t)-K*exp(-t/0.001);
g = @(t)K*exp(-t/0.008);
t = 0:0.00001:0.02;
figure
plot(t,f(t))
hold
plot(t, g(t))
plot(t,g(t)+f(t))