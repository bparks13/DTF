function [X]=create_data(N,m,stdZ,a)
%% [X]=create_data(N,m,stdZ,a)
%
%  Create a vector of values that is based on a known AR model. Used for testing accuracy
%  of AR coefficient estimation. If no inputs are given, uses N = 100000, m = 2, stdZ = 2,
%  a = [0.5 0.2]'.
%
%   Inputs:
%    - N: Number of samples
%    - m: AR model order
%    - stdZ: Standard deviation of the random Gaussian error injected into the signal
%    - a: Known AR coefficients for testing purposes. Size is [m x 1]
%
%   Outputs:
%    - X: Vector of signals created by the AR model. Size is [N x 1]
%
%  See also: mvar, estimate_ar_coeffecients, dtf, estimate_residuals
%

if nargin == 0
    N = 100000;   % Number of points
    m = 2;      % Desired order of the AR model
    stdZ = 2;   % Desired std of the random gaussean variable
    a = [0.5 0.2]';
end

if(m ~= length(a))
    disp('coefficient matrix must have m length')
    return
end

X = zeros(N+m,1);

for i = m+1:N
    X(i) = sum(X(i-1:-1:i-m).*a) + stdZ.*randn(1);
end

X=X(m+1:end);


end