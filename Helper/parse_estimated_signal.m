function x_hat=parse_estimated_signal(ar)
%% x_hat=parse_estimated_signal(ar)
%
%  Helper function that takes in the struct 'ar' containing the 'mdl' field, and returns a
%  simple matrix containing all the values from x_hat for a specific condition
%
%   Inputs:
%    - ar: Struct that is given by the function mvar, and contains the 'mdl' field with
%       all of the values from the model, but specifically contains the 'x_hat' field
%
%   Outputs:
%    - x_hat: Matrix containing the estimated signal values as created by mvar. Size is
%       [n x c x t], where n is the number of samples, c is the number of channels, and t
%       is the number of trials   
%
%  See also: mvar, plot_time_series, estimate_residuals
%

numTrials=length(ar);
numChannels=size(ar(1).mdl.x_hat,2);

max_size=size(ar(1).mdl.x_hat,1)+ar(1).mdl.order;

x_hat=nan(max_size,numChannels,numTrials);

for i=1:numTrials
    order=ar(i).mdl.order;
    
    x_hat(:,:,i)=[zeros(order,numChannels);ar(i).mdl.x_hat];
end

end