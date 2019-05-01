function [pass]=test_model(E,len)
%% [pass]=test_model(E,len)
%
% Given the residuals of the AR model, and either the length of the trial directly or the
% trials themselves to calculate the length, test the error for whiteness. In this case,
% for the model to 'pass' means that the null hypothesis is not rejected (pass = 0), where 
% the null hypothesis is that there is no autocorrelation between the residuals and the noise
% is therefore 'white'.
%
%   Inputs:
%    - E: Residuals, given from the output of estimate.
%    - len: Either an int that is the length of the trials use, or a vector/matrix of
%       signals used to create the model
%   
%   Outputs:
%    - pass: Boolean denoting whether all models do not reject the null hypothesis of the 
%       Portmanteau test for whiteness (pass = 0 means the null is not rejected)
%
% See also: varm, estimate, lbqtest, mvar
%

if length(len) > 1
    len=length(len);
end

lags=round(log(len));

numSeries=size(E,2);

h=ones(numSeries,1);

for i=1:numSeries
    h(i)=lbqtest(E(:,i),'lags',lags);
end

pass=any(h);

end