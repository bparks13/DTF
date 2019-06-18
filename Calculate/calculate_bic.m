function bic=calculate_bic(logL,modelOrder,N)
%% bic=calculate_bic(logL,modelOrder,N)
%
%  Calculate the Bayesian Information Criterion for the AR model
%
%   Inputs:
%    - logL: Log-likelihood of the model
%    - modelOrder: Model order of the model
%    - N: Number of samples in the original dataset
%
%   Outputs:
%    - bic: Bayesian Information Criterion
%
%  See also: calculate_loglikelihood, estimate_residuals, estimate_ar_coefficients
%

bic = -2 * logL + modelOrder * log(N - modelOrder);
% bic = 2 * logL + 2 * (modelOrder + 1) / N / log(N);

end