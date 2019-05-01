function [estMdl,E,criterion]=mvar(x,config)
%% [estMdl,E,criterion]=mvar(x,config)
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
%       [default], 2 (model order and criterion tested)
%
%   Outputs:
%    - estMdl: Output from varm and estimate, containing the estimated AR model
%       coefficients, among other things
%    - E: Residuals of the model fit, used for testing the whiteness of the model
%    - criterion: Optional output, contains the information criterion chosen in config
%
% See also: varm, estimate, summarize, test_model
%

% Defaults
orderRange=1:30;
crit='bic';
output=1;

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
end

numSeries=size(x,2);

criterion=zeros(max(orderRange),1);

minCrit=inf;

if output ~= 0
    fprintf('Beginning order estimation:\n');
end

for i=orderRange
    mdl=varm(numSeries,i);

    [estMdl,~,~,~]=estimate(mdl,x);

    results=summarize(estMdl);
    
    if strcmp(crit,'bic')
        criterion(i)=results.BIC;
        if results.BIC < minCrit
            minCrit=results.BIC;
        end
    elseif strcmp(crit,'aic')
        criterion(i)=results.AIC;
        if results.AIC < minCrit
            minCrit=results.AIC;
        end
    end
    
    if output == 1
        fprintf('%d,',i);
    elseif output == 2
        fprintf('%d: Current %s = %.2f. Minimum %s = %.2f\n',i,crit,criterion(i),crit,minCrit);
    end
end

modelOrder=find(criterion==minCrit);

if output ~= 0
    fprintf('\nDone: Minimum %s found at model order %d\n',crit,modelOrder);
end

mdl=varm(numSeries,modelOrder);
[estMdl,~,~,E]=estimate(mdl,x);

end

