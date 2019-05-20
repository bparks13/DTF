function plot_criterion(criterion,config)
%% plot_criterion(criterion,config)
%
%  Plot the criterion used to determine the model order. Assumed to be the criterion for a
%  single trial.
%
%   Inputs:
%    - criterion: The information criterion used (AIC/BIC/etc.) in a vector
%    - config: Optional struct containing additional plotting parameters
%       hFig: Handle to a figure to plot multiple criterion on the same plot
%       orderSelection: String defining which method was used to calculate the model
%           order; 'min' uses the minimum [default], and 'diff' uses the first time the
%           absolute difference is smaller than some epsilon, either 0.01 or user defined.
%       epsilon: Can be specified if orderSelection is set to 'diff', where epsilon is the
%           threshold for defining when differences have decreased to a small enough
%           degree to select the model order
%       modelOrder: The model that was picked no matter which method was used. If defined,
%           overrides manual selection of the model order based on criterion and method
%       figTitle: String containing a figure title, containing for example the condition
%           being tested, or the patient/date/run combo, or all of the above
%
%   Outputs:
%    Figure containing the plot of the criterion
%

bool_newFig=true;
epsilon=0.01;
figTitle='';

if nargin > 1 
    if isstruct(config)
        if isfield(config,'hFig')
            bool_newFig=false;
        end
        
        if isfield(config,'modelOrder')
            critInd=config.modelOrder;
        elseif isfield(config,'orderSelection')
            if isfield(config,'epsilon')
                epsilon=config.epsilon;
            end
            
            if strcmp(config.orderSelection,'min')
                critInd=find(criterion==min(criterion));
            elseif strcmp(config.orderSelection,'diff')
                critInd=find(abs(diff(criterion)) < epsilon,1);
            end
        else
            critInd=find(criterion==min(criterion));
        end
        
        if isfield(config,'figTitle')
            figTitle=config.figTitle;
        end
    end
end

if bool_newFig
    figure;
else
    figure(config.hFig);
end

plot(criterion,'k-o'); hold on;

plot(critInd,criterion(critInd),'r*');

if ~isempty(figTitle)
    tmp_hFig=gcf;
    tmp_hFig.Name=figTitle;
end

drawnow;

end