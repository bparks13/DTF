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

bic=-2*logL+modelOrder*log(N);

end