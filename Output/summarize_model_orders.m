function [meanOrder]=summarize_model_orders(ar)
%% [meanOrder]=summarize_model_orders(ar)
%
%  Given the ar output for all trials in a particular condition, print out the range
%  (min-max), and mean ± std for model orders
%
%   Inputs:
%    - ar: Struct containing the output from mvar, for all the trials in a particular
%       condition
%       mdl: Model struct, containing relevant fields
%           order: This is the only necessary field from the mdl output, should be an integer
%               denoting the model order that was used for the AR model
%
%   Outputs:
%    - meanOrder: Optional output. If an output variable is defined, the mean value of the
%       model orders is returned
%    Prints the output to the command window if there are no outputs defined
%
%  See also: mvar, dtf, estimate_ar_coefficients
%

if ~isfield(ar,'mdl')
    fields=fieldnames(ar);
    
    for i=1:length(fields)
        fprintf('%s: ',fields{i});
        summarize_model_orders(ar.(fields{i}));
    end
    
    return
end

numOrders=length(ar);
orders=zeros(numOrders,1);

for i=1:numOrders
    orders(i)=ar(i).mdl.order;
end

if nargout~=0
    meanOrder=mean(orders);
else
    fprintf('Range - [%2d, %2d], average model order - (%5.2f ± %-5.2f)\n',...
        min(orders),max(orders),mean(orders),std(orders));
end


end