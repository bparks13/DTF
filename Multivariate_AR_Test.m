CCC;

%% Create ground-truth models

N = 100000;                 % Number of samples
m = 2;                      % Model order
stdZ = 1;                   % Standard deviation of the error term
d = 2;                      % Number of time series
% A = [0.5,0.001;0.003,0.4];  % MAR coefficients. Should be [d x d x m]
A = zeros(d,d,m);
I = eye(d);

for i = 1:m
    randAR = rand(d);
    A(:,:,i) = I.*(randAR - 0.5) + (1 - I).*((randAR - 0.5) * 0.05);
end

X = zeros(N + m,d);             % Data, should be [N x d]. 
%   NOTE: Initiating with N + m so that the zeros at the beginning can be discarded

for i = m+1:N
    for j = 1:m
        X(i,:) = X(i,:) + X(i-j,:) * A(:,:,j) + stdZ*randn(1,d);
    end
end

X = X(m+1:end,:);

%% Estimate the AR coefficients

R_all=zeros(m*2+1,m^2);

% [R_xx,lags_xx] = xcov(X(:,1),m,'unbiased');
% [R_yy,lags_yy] = xcov(X(:,2),m,'unbiased');
% [R_xy,lags_xy] = xcov(X(:,1),X(:,2),m,'unbiased');
% [R_yx,lags_yx] = xcov(X(:,2),X(:,1),m,'unbiased');

for i=1:m^2
    a=ceil(i/2);
    b=mod(i-1,2)+1;
    [R_all(:,i),lags] = xcov(X(:,a),X(:,b),m,'unbiased');
end

r_n = zeros(m,m*d);

for i=1:m
    for j=1:m*d
        a=(i-1)*2 + mod(j-1,2)+1;
        b=-ceil(j/2);
        r_n(i,j)=R_all(lags==b,a);
    end
end

% r_n(1,1)=R_xx(lags_xx==-1);
% r_n(1,2)=R_xy(lags_xy==-1);
% r_n(1,3)=R_xx(lags_xx==-2);
% r_n(1,4)=R_xy(lags_xy==-2);
% r_n(2,1)=R_yx(lags_yx==-1);
% r_n(2,2)=R_yy(lags_yy==-1);
% r_n(2,3)=R_yx(lags_yx==-2);
% r_n(2,4)=R_yy(lags_yy==-2);

R_n = zeros(m*d);

for i=1:m*d
    for j=1:m*d
        a=mod((i-1)*2+mod(j-1,2),4)+1;
        b=-floor((j-1)/2)*mod(ceil(i/2),2)+mod(ceil(j/2),2)*mod(ceil(i/2)-1,2);
        R_n(i,j)=R_all(lags==b,a);
    end
end

% R_n(1,1)=R_xx(lags_xx==0);
% R_n(1,2)=R_xy(lags_xy==0);
% R_n(1,3)=R_xx(lags_xx==-1);
% R_n(1,4)=R_xy(lags_xy==-1);
% R_n(2,1)=R_yx(lags_yx==0);
% R_n(2,2)=R_yy(lags_yy==0);
% R_n(2,3)=R_yx(lags_yx==-1);
% R_n(2,4)=R_yy(lags_yy==-1);
% R_n(3,1)=R_xx(lags_xx==1);
% R_n(3,2)=R_xy(lags_xy==1);
% R_n(3,3)=R_xx(lags_xx==0);
% R_n(3,4)=R_xy(lags_xy==0);
% R_n(4,1)=R_yx(lags_yx==1);
% R_n(4,2)=R_yy(lags_yy==1);
% R_n(4,3)=R_yx(lags_yx==0);
% R_n(4,4)=R_yy(lags_yy==0);

Theta_hat = reshape(r_n / R_n,d,d,m);

A
Theta_hat
A - Theta_hat

