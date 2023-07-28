% UKF for Jansen and Rit

function [xhat, Pxx, K, Xa] = UKFJR(xhat,Pxx,y,fct,obsfct,parameters,R,Sigma)

global dx 


kappa = 3-dx;

N = 2*dx;       % number of sigma points
Pxx = (Pxx + Pxx')/2;   % make sure it is symmetric for numerical stability

xsigma = chol(dx*Pxx)'; % Pxx = chol(Pxx)'*chol(Pxx) and sigmas  = sqrt(dx)sqrt(P)

Xa = xhat*ones(1,N) + [xsigma, -xsigma];    % these are the sigma points!!

% now we propogate the sigma points through the system eqs.
X = feval(fct,parameters,Xa);       % note: fct = 'vossFNfct';

xtilde = sum(X,2)/N;        % mean of the propogated sigma points
X1 = X - xtilde*ones(1,size(X,2)); % substract the mean from the columns of X

% this is different to the Schiff / Voss example where I add Sigma to the
% covariance
% ~~~~~~~~~~~~~~~~~
Pxx = X1*X1' /N + Sigma;
% ~~~~~~~~~~~~~~~~~

Pxx = (Pxx + Pxx')/2;       % make sure we have a symetric matrix again

% ~~~~~~~~~~~~~~~~~~~~~~
Y=feval(obsfct,X);

ytilde = sum(Y,2)/N;            % the mean
Y1 = Y - ytilde*ones(1,size(Y,2)); 
Pyy = Y1*Y1' /N + R;

Pxy = X1*Y1' /N;

% calculate the Kalman gain and update
K = Pxy/Pyy;
xhat=xtilde+K*(y-ytilde);
Pxx=Pxx-K*Pxy';
Pxx = (Pxx + Pxx')/2;
