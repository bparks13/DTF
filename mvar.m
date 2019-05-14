function [mdl,E,criterion]=mvar(x,config)
%% [mdl,E,criterion]=mvar(x,config)
%  
%  Given a matrix of signals, and optionally a configuration struct, construct a
%  multi-variate autoregressive model of the signals.
%
%   Inputs:
%    - x: Signals in a matrix format [n x s], where n is the number of samples and s is
%       the number of series to model
%    - config: Struct containing optional parameters
%       orderRange: Vector containing the model orders to consider. Default [1:30]
%       crit: String defining which information criterion to use, 'aic' or 'bic' [default]
%       output: Int defining level of verbosity for output. 0 (none), 1 (model number)
%           [default], 2 (model order and criterion tested)
%       method: String defining which method to use; Matlab's varm ('varm') or Yule-Walker
%           equations ('yule') [default].
%
%   Outputs:
%    - mdl: Struct containing the AR model fit for the data given
%       AR: Autoregressive coefficients as found by estimate_ar_coefficients
%       C: Covariance matrx as calculated by the Yule-Walker equations
%       logL: Log-likelihood of the model fit, used for calculating information criterion
%       order: Model order that is found to have the lowest information criterion
%       numSeries: Number of series in the model
%    - E: Residuals of the model fit, used for testing the whiteness of the model. Size is
%       [(n - o) x s], where n is the number of samples, o is the model order, and s is
%       the number of series
%    - criterion: Optional output, contains the information criterion chosen in config
%
% See also: estimate_ar_coefficients, estimate_residuals, calculate_loglikelihood,
% calculate_bic
%

% Defaults
orderRange=1:30;
crit='bic';
output=1;
method='yule';

if nargin > 1 && isstruct(config)
    if isfield(config,'orderRange')
        orderRange=config.orderRange;
    end
    
    if isfield(config,'crit')
        crit=config.crit;
    end
    
    if isfield(config,'output')
        output=config.crit;
    end
    
    if isfield(config,'method')
        method=config.method;
    end
end

numSeries=size(x,2);
numOrders=length(orderRange);
numSamples=length(x);

criterion=zeros(numOrders,1);

minCrit=inf;

mdl=struct('AR',[],'C',[],'logL',[],'order',[],'numSeries',numSeries);

if output ~= 0
    fprintf('Beginning order estimation:\n');
end

for i=1:numOrders
    if strcmp(method,'varm')
        mdl=varm(numSeries,i);
        [estMdl,~,logL,E]=estimate(mdl,x);
        results=summarize(estMdl);
        
        if strcmp(crit,'bic')
            criterion(i)=results.BIC;
            
            if results.BIC < minCrit
                minCrit=results.BIC;
                mdl.AR=zeros(numSeries,numSeries,orderRange(i));
                for j=1:orderRange(i)
                    mdl.AR(:,:,j)=reshape(estMdl.AR{j});
                end
                mdl.C=estMdl.Covariance;
                mdl.logL=logL;
                mdl.order=orderRange(i);
            end
        elseif strcmp(crit,'aic')
            criterion(i)=results.AIC;
            
            if results.AIC < minCrit
                minCrit=results.AIC;
            end
        end
    elseif strcmp(method,'yule')
        [AR,C]=estimate_ar_coefficients(x,orderRange(i));
        E=estimate_residuals(x,AR);
        logL=calculate_loglikelihood(E,C);
        criterion(i)=calculate_bic(logL,orderRange(i),numSamples-orderRange(i));

        if strcmp(crit,'bic')
            if criterion(i) < minCrit
                minCrit=criterion(i);
                mdl.AR=AR;
                mdl.C=C;
                mdl.logL=logL;
                mdl.order=orderRange(i);
            end
        elseif strcmp(crit,'aic')
            disp('WARNING: No AIC calculation implemented. No criterion tested.');
        end    
    end
    
    
    if output == 1
        fprintf('%d,',i);
    elseif output == 2
        fprintf('%d: Current %s = %.2f. Minimum %s = %.2f\n',i,crit,criterion(i),crit,minCrit);
    end
end

if output ~= 0
    fprintf('\nDone: Minimum %s found at model order %d\n',crit,mdl.order);
end

end

