function print_whiteness(h,p,labels)
%% print_whiteness(h,p,labels)
%
%  Prints the individual trials p-values when the null hypothesis of uncorrelated noise is
%  rejected 
%   
%   Inputs:
%    - h: Hypothesis values, where 0 means the null hypothesis is not rejected, and 1 is
%       that the null hypothesis is rejected. Size is [t x s], where t is the number of
%       trials, and s is the number of series
%    - p: P-values associated with the corresponding hypothesis value. Size is [t x s],
%       where t is the number of trials, and s is the number of series
%    - labels: Cell array of the channel names
%
%   Outputs:
%    Prints in the command window the individual trials p-values if the null is rejected
%
%   See also: test_model, summarize_whiteness
%

numTrials=size(h,1);
numChannels=size(h,2);

for i=1:numTrials
    pass=~any(h(i,:));
    
    if ~pass
        fprintf('WARNING: Null hypothesis of uncorrelated errors rejected for trial %d\n',i);
        for k=1:numChannels
            fprintf('\t%s: h = %d with p = %.4f\n',labels{k},h(i,k),p(i,k));
        end
    end
end

end