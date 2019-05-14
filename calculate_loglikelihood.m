function logL=calculate_loglikelihood(E,C)
%% logL=calculate_loglikelihood(E,C)
%
%  Calculate the log-likelihood of the model using the residuals and the covariance
%
%   Inputs:
%    -  E: Residuals between the data and the model (AR), size is [N - m x d], where the N
%       is the number of points in the original data, m is the model order, and d is the
%       number of series
%    - C: Covariance as calculated from the Yule-Walker Equations, in a [d x d] matrix
%
%   Outputs:
%    - logL: Log-likelihood of the model
%
%  See also: infer, estimate_residuals, estimate_ar_coefficients
%

[N,K]=size(E);  % N is number of samples, K is number of series
R = chol(C);    % Upper triangular Cholesky factor

logL = -0.5 * sum(sum((E / R).^2,2)) - 0.5 * N * (K * log(2*pi) + 2 * log(det(R)));

end