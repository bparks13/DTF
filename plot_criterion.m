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
%
%   Outputs:
%    Figure containing the plot of the criterion
%

bool_newFig=true;

if nargin > 1 
    if isstruct(config)
        if isfield(config,'hFig')
            bool_newFig=false;
        end
    end
end

if bool_newFig
    figure;
else
    figure(config.hFig);
end

plot(criterion,'k-o'); hold on;

minCritInd=find(criterion==min(criterion));
plot(minCritInd,criterion(minCritInd),'r*');

end