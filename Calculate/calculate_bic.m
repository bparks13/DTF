function bic=calculate_bic(logL,modelOrder,numChannels,N,method)
%% bic=calculate_bic(logL,modelOrder,numChannels,N,method)
%
%  Calculate the Bayesian Information Criterion for the AR model
%
%   Inputs:
%    - logL: Log-likelihood of the model
%    - modelOrder: Model order of the model
%    - numChannels: Number of channels used in the model
%    - N: Number of samples in the original dataset
%    - method: Int defining whether to use the Matlab Log-Likelihood [1] or the Ding
%       Log-Likelihood [2], or the Awareness during anaesthesia paper method [3]
%
%   Outputs:
%    - bic: Bayesian Information Criterion
%
%  See also: calculate_loglikelihood, estimate_residuals, estimate_ar_coefficients
%

numParams=numChannels^2 * (modelOrder + 1);

if method == 1
    bic = -2 * logL + numParams * log(N - modelOrder);
elseif method == 2
    bic = logL + ((2 * numParams * log(N)) / N);
elseif method == 3
    bic = logL + numParams * log(N - modelOrder);
end

end