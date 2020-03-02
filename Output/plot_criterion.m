function plot_criterion(criterion,config)
%% plot_criterion(criterion,config)
%
%  Plot the criterion used to determine the model order.
%
%   Inputs:
%    - criterion: The information criterion used (AIC/BIC/etc.) in a vector. If criterion
%       is given as a struct containing criteria from all trials, then it will
%       automatically go through all trials and plot them in the same figure
%    - config: Optional struct containing additional plotting parameters
%    -- hFig: Handle to a figure to plot multiple criterion on the same plot
%    -- hAx: Handle to an existing axis, typically a subplot
%    -- orderSelection: String defining which method was used to calculate the model
%        order; 'min' uses the minimum [default], and 'diff' uses the first time the
%        absolute difference is smaller than some epsilon, either 0.001 or user defined.
%    -- epsilon: Can be specified if orderSelection is set to 'diff', where epsilon is the
%        threshold for defining when differences have decreased to a small enough
%        degree to select the model order
%    -- modelOrder: The model that was picked no matter which method was used. If defined,
%        overrides manual selection of the model order based on criterion and method
%    -- figTitle: String containing a figure title, containing for example the condition
%        being tested, or the patient/date/run combo, or all of the above
%    -- average: Boolean defining whether or not to only plot the average values. Default
%        is false.
%
%   Outputs:
%    Figure containing the plot of the criterion
%
%  See also: mvar, calculate_loglikelihood, calculate_aic, calculate_bic,
%   calculate_ar_spectra
%

bool_newFig=true;
bool_newAx=true;
bool_plotModelOrder=false;
bool_criteriaKnown=false;
yLabel='';
epsilon=0.001;
figTitle='';
average=false;

if nargin > 1 
    if isstruct(config)
        if isfield(config,'hFig') && ~isempty(config.hFig)
            bool_newFig=false;
        end
        
        if isfield(config,'hAx') && ~isempty(config.hAx)
            bool_newAx=false;
        end
        
        if isfield(config,'modelOrder') && ~isempty(config.modelOrder)
            critInd=config.modelOrder;
            bool_plotModelOrder=true;
        elseif isfield(config,'orderSelection') && ~isempty(config.orderSelection)
            if isfield(config,'epsilon') && ~isempty(config.epsilon)
                epsilon=config.epsilon;
            end
            
            if strcmp(config.orderSelection,'min')
                critInd=find(criterion==min(criterion));
                bool_plotModelOrder=true;
            elseif strcmp(config.orderSelection,'diff')
                critInd=find(abs(diff(criterion)) < epsilon,1);
                bool_plotModelOrder=true;
            end
        end
        
        if isfield(config,'figTitle') && ~isempty(config.figTitle)
            figTitle=config.figTitle;
        end
        
        if isfield(config,'crit') && ~isempty(config.crit)
            if strcmp(config.crit,'aic')
                yLabel='AIC';
                bool_criteriaKnown=true;
            elseif strcmp(config.crit,'bic')
                yLabel='BIC';
                bool_criteriaKnown=true;
            elseif strcmp(config.crit,'psd') || strcmp(config.crit,'spectra')
                yLabel='MSE';
                bool_criteriaKnown=true;
            end
        end
        
        if isfield(config,'average') && ~isempty(config.average)
            average=config.average;
        end
    end
end

if bool_newFig
    figure;
else
    figure(config.hFig);
end

if isstruct(criterion)
    field=fieldnames(criterion);
    
    if length(field) > 1
        warning('Unexpected number of fields given.');
    end
    
    if nargin == 1
        config.hFig=gcf;
    end
    
    if average
        tmp_crit=zeros(length(criterion(1).(field{1})),1);
        
        for i=1:length(criterion)
            tmp_crit=tmp_crit+criterion(i).(field{1});
        end
        
        tmp_crit=tmp_crit/length(criterion);
        
        config.figTitle=[figTitle,'_Average'];
        plot_criterion(tmp_crit,config);
    else
        for i=1:length(criterion)
            plot_criterion(criterion(i).(field{1}),config);
        end
    end
    
    return
else
    if average
        tmp_crit=zeros(size(criterion,1),1);
        
        for i=1:size(criterion,2)
            tmp_crit=tmp_crit+criterion(:,i);
        end
        
        tmp_crit=tmp_crit/size(criterion,2);
        criterion=tmp_crit;
        
        config.figTitle=[figTitle,'_Average'];
    end
end

if bool_newAx
    plot(criterion,'k-o'); hold on;

    if bool_plotModelOrder
        plot(critInd,criterion(critInd),'r*');
    end
else
    plot(config.hAx,criterion,'k-o'); hold on;

    if bool_plotModelOrder
        plot(config.hAx,critInd,criterion(critInd),'r*');
    end
end

xlabel('Model Order'); 

if bool_criteriaKnown
    ylabel(yLabel);
end

if ~isempty(figTitle)
    tmp_hFig=gcf;
    tmp_hFig.Name=figTitle;
end

drawnow;

end