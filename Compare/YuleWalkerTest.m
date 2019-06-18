%%
%Part One: Develop and m AR model
%X[n] = a[1]*X[n-1] + a[2]*X[n-2] + ... + a[m]*X[n-m] + Z[n]

clear; 
close all;
clc

N = 100000;   % Number of points
m = 2;      % Desired order of the AR model
stdZ = 2;   % Desired std of the random gaussean variable
a = [0.5 0.2]';

if(m ~= length(a))
    disp('coefficient matrix must have m length')
    return
end

X = zeros(N,1);

for i = m+1:N
    X(i) = sum(X(i-1:-1:i-m).*a) + stdZ.*randn(1);
end

figure
plot(X)

%% 
% Part 2: Yule Walker Estimation
[rxx,lags] = xcov(X,m,'unbiased');
r_n = rxx(lags>0);
idx0 = find(lags==0);
R_n = zeros(m,m);

for i=1:m
    tmpidx = idx0 -(i-1);
    R_n(:,i) = rxx(tmpidx:tmpidx+m-1);
end

Theta_hat = R_n\r_n; % Predictor for AR coefficients
SigmaZ_hat = rxx(lags>=0)'*[1; -Theta_hat];
