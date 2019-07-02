function aic=calculate_aic(logL,modelOrder,N)
%% aic=calculate_aic(logL,modelOrder,N)
%
%  Calculate the Akaike Information Criterion for the AR model
%
%   Inputs:
%    - logL: Log-likelihood of the model
%    - modelOrder: Model order of the model
%    - N: Number of samples in the original dataset
%
%   Outputs:
%    - bic: Bayesian Information Criterion
%
%  See also: calculate_loglikelihood, estimate_residuals, estimate_ar_coefficients,
%       calculate_bic
%

aic = -2 * logL + modelOrder * 2;

end