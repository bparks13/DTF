function aic=calculate_aic(logL,modelOrder,numChannels,N,method)
%% aic=calculate_aic(logL,modelOrder,N,method)
%
%  Calculate the Akaike Information Criterion for the AR model
%
%   Inputs:
%    - logL: Log-likelihood of the model
%    - modelOrder: Model order of the model
%    - numChannels: Number of channels used in the model
%    - N: Number of samples in the original dataset
%    - method: Int defining whether to use the Matlab Log-Likelihood [1] or the Ding
%       Log-Likelihood [2]
%
%   Outputs:
%    - aic: Akaike Information Criterion
%
%  See also: calculate_loglikelihood, estimate_residuals, estimate_ar_coefficients,
%       calculate_bic
%

numParams = numChannels^2 * (modelOrder);

if method == 1
    aic = -2 * logL + modelOrder * 2;
elseif method == 2
    aic = 2 * logL + ((2 * numParams) / (N));
end

end