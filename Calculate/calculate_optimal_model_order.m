function [modelOrder,criterion]=calculate_optimal_model_order(x,config)
%% [modelOrder,criterion]=calculate_optimal_model_order(x,config)
%
%  Instead of calculating the optimal order inside the mvar.m function, find the order
%  here first, and apply that one model order to every realization within a single
%  condition (does not apply between conditions). 
%
%   Inputs:
%    - x: Signals in a matrix format [n x c x r], where n is the number of samples, c is
%       the number of channels to model, and r is the number of realizations
%    - config: Struct containing additional parameters
%    -- orderRange: Vector containing the model orders to consider. 
%    -- crit: String defining which information criterion to use, 'aic', 'bic',
%           'psd', or 'spectra'
%    -- method: String defining which method to use; Matlab's varm ('varm') or Yule-Walker
%           equations ('yule'). Additionally can use the Signal Processing for
%           Neuroscientists method ('arfit')
%    -- orderSelection: String defining how to algorithmically choose the model order.
%           'min' uses the minimum information criterion found, while 'diff1'
%           uses the first model order that the abs(difference) between successive ICs is
%           smaller than some epsilon, which is either 0.01 or user specified. 'diff2'
%           uses the same criteria as 'diff1', but chooses the model order (m-1) from when
%           the criteria is met; this is the difference betweeen choosing the first bend
%           in the criterion and the second point after the bend
%    -- epsilon: Can be specified if orderSelection is set to 'diff', where epsilon is the
%           threshold for defining when differences have decreased to a small enough
%           degree to select the model order, defined as a percentage (i.e. 0.01 == 1%).
%           Default value is 0.1% (= 0.001). Does not apply if using the 'arfit' method
%
%   Outputs:
%    - modelOrder: A scalar defining the optimal order defined by averaging the given
%       criterion over all realizations, and choosing based on the orderSelection given
%    - criterion: Returns all the criterion values found while searching for the optimal,
%       using the given criteria.
%
%  See also: pipeline.m, mvar.m, calculate_loglikelihood.m
%

if nargin == 1
    error('config is currently required; no defaults are set for this function');
end

numRealizations=size(x,3);

criterion=zeros(length(config.orderRange),numRealizations);

if strcmp(config.method,'arfit')
    if strcmp(config.crit,'bic')
        for i=1:numRealizations
            [~,~,~,criterion(:,i),~,~]=arfit(x(:,:,i),config.orderRange(1),config.orderRange(end)); % SBC (BIC)
        end
    elseif strcmp(config.crit,'aic')
        for i=1:numRealizations
            [~,~,~,~,criterion(:,i),~]=arfit(x(:,:,i),config.orderRange(1),config.orderRange(end)); % FPE (AIC)
        end
    end
end

crit=mean(criterion,2);

if strcmp(config.orderSelection,'min')
    [~,modelOrder]=min(crit);
elseif strcmp(config.orderSelection,'diff1')
    ind=find(abs((crit(2:end) - crit(1:end-1)) ./ crit(2:end)) < config.epsilon | crit(2:end) - crit(1:end-1) > 0,1)+1;
    modelOrder=config.orderRange(ind);
elseif strcmp(config.orderSelection,'diff2')
    ind=find(abs((crit(2:end) - crit(1:end-1)) ./ crit(2:end)) < config.epsilon | crit(2:end) - crit(1:end-1) > 0,1)+1;
    modelOrder=config.orderRange(ind-1);
end

end