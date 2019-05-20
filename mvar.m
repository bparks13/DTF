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
%           equations ('yule') [default]
%       orderSelection: String defining how to algorithmically choose the model order.
%           'min' [default] uses the minimum information criterion found, while 'diff'
%           uses the first model order that the abs(difference) between successive ICs is
%           smaller than some epsilon, which is either 0.01 or user specified
%       epsilon: Can be specified if orderSelection is set to 'diff', where epsilon is the
%           threshold for defining when differences have decreased to a small enough
%           degree to select the model order
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
orderSelection='min';
epsilon=0.01;

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
    
    if isfield(config,'orderSelection')
        orderSelection=config.orderSelection;
        
        if isfield(config,'epsilon')
            epsilon=config.epsilon;
        end
    end
end

numSeries=size(x,2);
numOrders=length(orderRange);
numSamples=length(x);

criterion=zeros(numOrders,1);

minCrit=inf;
bool_minDiffFound=false;

mdl=struct('AR',[],'C',[],'logL',[],'order',[],'numSeries',numSeries);

if output ~= 0
    fprintf('Beginning order estimation:\n');
end

for i=1:numOrders
    if strcmp(method,'varm')
        mdl_varm=varm(numSeries,i);
        [estMdl,~,logL,tmp_E]=estimate(mdl_varm,x);
        results=summarize(estMdl);
        
        if strcmp(crit,'bic')
            criterion(i)=results.BIC;
            
            if strcmp(orderSelection,'min')
                result = results.BIC < minCrit;
            elseif strcmp(orderSelection,'diff')
                if i>1 && ~bool_minDiffFound
                    result = abs(criterion(i) - criterion(i-1)) < epsilon || ...
                        criterion(i) - criterion(i-1) > 0;
                    bool_minDiffFound=true;
                else
                    result=false;
                end
            end
            
            if result
                minCrit=results.BIC;
                mdl.AR=zeros(numSeries,numSeries,orderRange(i));
                for j=1:orderRange(i)
                    mdl.AR(:,:,j)=estMdl.AR{j};
                end
                mdl.C=estMdl.Covariance;
                mdl.logL=logL;
                mdl.order=orderRange(i);
                E=tmp_E;
                
                if strcmp(orderSelection,'diff')
                    bool_minDiffFound=true;
                end
            end
        elseif strcmp(crit,'aic')
            criterion(i)=results.AIC;
            
            if results.AIC < minCrit
                minCrit=results.AIC;
            end
        end
    elseif strcmp(method,'yule')
        [AR]=estimate_ar_coefficients(x,orderRange(i));
        [tmp_E,C]=estimate_residuals(x,AR);
        logL=calculate_loglikelihood(tmp_E,C);
        criterion(i)=calculate_bic(logL,orderRange(i),numSamples-orderRange(i));

        if strcmp(crit,'bic')
            if strcmp(orderSelection,'min')
                result = criterion(i) < minCrit;
            elseif strcmp(orderSelection,'diff')
                if i>1 && ~bool_minDiffFound
                    result = abs(criterion(i) - criterion(i-1)) < epsilon;
                else
                    result=false;
                end
            end
            
            if result
                minCrit=criterion(i);
                mdl.AR=AR;
                mdl.C=C;
                mdl.logL=logL;
                mdl.order=orderRange(i);
                E=tmp_E;
                
                if strcmp(orderSelection,'diff')
                    bool_minDiffFound=true;
                end
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

if isempty(mdl.order)
    fprintf('\nWARNING: NO MODEL FOUND WITH ORDER SELECTION ''%s''\nFINDING MINIMUM VALUE...\n',orderSelection);
    
    minCritInd=find(criterion==min(criterion));
    
    [mdl.AR]=estimate_ar_coefficients(x,minCritInd);
    [E,mdl.C]=estimate_residuals(x,mdl.AR);
    mdl.logL=calculate_loglikelihood(E,mdl.C);
    mdl.order=minCritInd;
end

if output ~= 0
    fprintf('\nDone: Minimum %s found at model order %d\n',crit,mdl.order);
end

end

