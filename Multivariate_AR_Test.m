CCC;

%% Create ground-truth models

N = 10000;                 % Number of samples
m = 30;                      % Model order
stdZ = 1;                   % Standard deviation of the error term
d = 8;                      % Number of time series
% A = [0.5,0.001;0.003,0.4];  % MAR coefficients. Should be [d x d x m]
A = zeros(d,d,m);
I = eye(d);

for i = 1:m
    randAR = rand(d);
    A(:,:,i) = I.*((randAR) - 0.5) / i + (1 - I).*((randAR - 0.5) * 0.05);
end

X = zeros(N + m,d);             % Data, should be [N x d]. 
%   NOTE: Initiating with N + m so that the zeros at the beginning can be discarded

for i = m+1:N+m
    for j = 1:m
        X(i,:) = X(i,:) + X(i-j,:) * A(:,:,j) + stdZ*randn(1,d);
    end
end

X = X(m+1:end,:);

%% Estimate the AR coefficients

R_all=zeros(m*2+1,d^2);

for i=1:d^2
    a=ceil(i/d);
    b=mod(i-1,d)+1;
    [R_all(:,i),lags] = xcov(X(:,a),X(:,b),m,'unbiased');
end

r_n = zeros(d,m*d);

for i=1:d
    for j=1:m*d
        a=(i-1)*d + mod(j-1,d)+1;
        b=-ceil(j/d);
%         fprintf('(%d,%d) a = %d, b = %d\n',i,j,a,b);
        r_n(i,j)=R_all(lags==b,a);
    end
end

R_n = zeros(m*d);

for i=1:m*d
    for j=1:m*d
        a=mod((i-1)*d+mod(j-1,d),d^2)+1;
        b=-floor((j-1)/d)+floor((i-1)/d);
%         fprintf('(%d,%d) a = %d, b = %d\n',i,j,a,b);
        R_n(i,j)=R_all(lags==b,a);
    end
end

Theta_hat = reshape(r_n / R_n,d,d,m);

return

%%

A %#ok<*UNRCH>
Theta_hat
diff=A - Theta_hat;

