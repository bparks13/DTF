function [pass,h,pVal,cVal]=test_model(E,len,config)
%% [pass,h,pVal,cVal]=test_model(E,len,config)
%
% Given the residuals of the AR model, and either the length of the trial directly or the
% trials themselves to calculate the length, test the error for whiteness. In this case,
% for the model to 'pass' means that the null hypothesis is not rejected (pass = 1), where 
% the null hypothesis is that there is no autocorrelation between the residuals and the noise
% is therefore 'white'.
%
%   Inputs:
%    - E: Residuals, given from the output of estimate. Matrix of residuals, with size
%       [t x c], where t is the number of samples and c is the number of channels
%    - len: Either an int that is the length of the trials used, or a vector/matrix of
%       signals used to create the model
%    - config: Optional struct containing additional parameters on how to run the test
%       overrideLags: Boolean denoting whether to override the preset lags defined as the
%           natural log of the length of the data used. Default is false, where the log of
%           the length of the data is used as the number of lags
%       lags: Manual lags. Must be defined if overrideLags is true
%       changeDOF: Boolean defining if the degrees of freedom should be changed
%       numParameters: If changeDOF is true, this parameter must be defined. The new DOF
%           is then defined as the number of lags minus the number of parameters of the
%           model
%   
%   Outputs:
%    - pass: Boolean denoting whether all models do not reject the null hypothesis of the 
%       Portmanteau test for whiteness (pass = 1 means the null is not rejected)
%    - h: Hypothesis rejection values of each series. '0' means that the null hypothesis
%       is not rejected, while '1' means that the null hypothesis is rejected. Size is
%       [t x c], where t is the number of trials in that condition and c is the number of
%       channels
%    - pVal: P-value of the test for each series
%    - cVal: Optional output. Critical values of the hypothesis testing
%
% See also: mvar, estimate_ar_coefficients, estimate_residuals, calculate_loglikelihood,
%   calculate_bic
%

if length(len) > 1
    len=length(len);
end

lags=round(log(len));
dof=lags;

if nargin > 2
    if isstruct(config)
        if isfield(config,'overrideLags') && ~isempty(config.overrideLags)
            if config.overrideLags
                if isfield(config,'lags') && ~isempty(config.lags)
                    lags=config.lags;
                    dof=lags;
                end
            end
        end
        
        if isfield(config,'changeDOF') && ~isempty(config.changeDOF)
            if config.changeDOF
                if isfield(config,'numParameters') && ~isempty(config.numParameters)
                    dof=lags-config.numParameters;
                else
                    disp('WARNING: To implement a new DOF value, numParameters must be defined in config.')
                end
            end
        end
    end
end

numSeries=size(E,2);

h=ones(numSeries,1);
pVal=zeros(numSeries,1);

for i=1:numSeries
    [h(i),pVal(i),cVal]=portmanteau(E(:,i),lags,dof,0.05);
end

pass=~any(h);

if nargout == 3
    clear cVal
end

    function [h,p,c]=portmanteau(E,lags,dof,alpha)
    %% [h,p,c]=portmanteau(E,lags,dof,alpha)
    % 
    %  Internal function to use the Ljung-Box Test to test if the residuals exhibit
    %  autocorrelation or not, with the null hypothesis being that there is no
    %  autocorrelation for the set number of lags given. All calculations are pulled from
    %  function 'lbqtest', this is to work-around the licensing issue resulting from too
    %  many people trying to use the Econometrics toolbox.
    %
    %   Inputs:
    %    - E: Residuals of the model fit
    %    - lags: Number of lags to test. Default is min(20,log(length(E)))
    %    - dof: Degrees of freedom to use for the Chi-Square test. Default is lags, but is
    %       recommended to be decreased by the number of parameters used to fit the model
    %       (which is mostly impossible for my data, since the model order is typically
    %       equal to or higher than the log length of the trial)
    %    - alpha: Alpha value to test the hypothesis against
    %
    %   Outputs:
    %    - h: Hypothesis value, either 0 (null hypothesis) or 1 (alternative hypothesis)
    %    - p: P-value associated with the hypothesis value. Does not correct for multiple
    %       comparisons, and uses a default value of alpha = 0.05
    %    - c: Critical value for the chi-squared density function
    %
    %  See also: lbqtest
    %
    
    T=length(E);
    lag=min(20,round(log(T)));
    df=lag;
    
    if nargin == 2
        lag=lags;
        df=lags;
    end
    
    if nargin == 3 || nargin == 4
        lag=lags;
        df=dof;
    end
    
    [ACF,lagInd]=xcorr(E,lag,'coeff');
    ACF=ACF(lagInd>0);
    idx=(T-(1:lag))';
    stat=T*(T+2)*cumsum((ACF.^2)./idx);
    stat=stat(end);
    
    p=1-chi2cdf(stat,df);
    
    h=alpha >= p;
    
    c = chi2inv(1-alpha,df);
        
    end

end