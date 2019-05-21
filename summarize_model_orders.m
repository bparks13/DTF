function summarize_model_orders(cond)
%% summarize_model_orders(cond)
%
%  Given the mdl output for all trials in a particular condition, print out the range
%  (min-max), and mean ± std for model orders
%
%   Inputs:
%    - cond: Struct containing the output from mvar, for all the trials in a particular
%       condition
%       order: This is the only necessary field from the mdl output, should be an integer
%           denoting the model order that was used for the AR model
%
%   Outputs:
%    Prints the output to the command window
%
%  See also: mvar, dtf, estimate_ar_coefficients
%

numOrders=length(cond);
orders=zeros(numOrders,1);

for i=1:numOrders
    orders(i)=cond(i).mdl.order;
end

fprintf('Range - [%d, %d], average model order - (%.2f ± %.2f)\n',...
    min(orders),max(orders),mean(orders),std(orders));

end