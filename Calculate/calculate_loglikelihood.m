function logL=calculate_loglikelihood(E,C,method)
%% logL=calculate_loglikelihood(E,C,method)
%
%  Calculate the log-likelihood of the model using the residuals and the covariance
%
%   Inputs:
%    -  E: Residuals between the data and the model (AR), size is [N - m x d], where the N
%       is the number of points in the original data, m is the model order, and d is the
%       number of series
%    - C: Covariance as calculated from the Yule-Walker Equations, in a [d x d] matrix
%    - method: Int defining whether to use the Matlab Log-Likelihood [1] or the Ding
%       Log-Likelihood [2]
%
%   Outputs:
%    - logL: Log-likelihood of the model
%
%  See also: infer, estimate_residuals, estimate_ar_coefficients
%

% Matlab log-likelihood comes from internal functions, \matlab\toolbox\stats\stats\private\statmvnrobj.m
% Ding log-likelihood comes from https://onlinelibrary.wiley.com/doi/abs/10.1002/9783527609970.ch17

if nargin < 3
    method=1;
end

if nargin==1 || isempty(C)
    C=(E' * E) / length(E);
end

if method == 1
    [N,K]=size(E);  % N is number of samples, K is number of series
    R = chol(C);    % Upper triangular Cholesky factor

    logL = -0.5 * sum(sum((E / R).^2,2)) - 0.5 * N * (K * log(2*pi) + 2 * log(det(R)));
elseif method == 2 || method == 3
    logL=log(det(C));
end

end