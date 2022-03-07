
dt = 0.001;
T = 1;
t = linspace(0,T,1/dt);
L = length(t);

state = 1;
% x_ = x(state,:);
x_ = y;
% x_ = v_ip(1,:);
x_ = x_ - mean(x_);

Fs = 1/dt;
L = length(t);
n = 2^nextpow2(L);
X = fft(x_,n);
P2 = abs(X/L);
P1 = P2(:,1:n/2+1);
P1(:,2:end-1)=2*P1(:,2:end-1);

figure; 
plot(0:(Fs/n):(Fs/2-Fs/n), P1(1:n/2));