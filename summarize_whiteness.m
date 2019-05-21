function summarize_whiteness(h,labels)
%% summarize_whiteness(h,labels)
%
%  Given the hypothesis values and the p-values, prints in the command window the number
%  of times each channel rejects the null
%
%   Inputs:
%    - h: Hypothesis values, where 0 means the null hypothesis is not rejected, and 1 is
%       that the null hypothesis is rejected. Size is [t x s], where t is the number of
%       trials, and s is the number of series
%    - labels: Cell array of the channel names
%
%   Outputs:
%    Prints in the command window a summary of the number of times each channel rejects
%       the null
%
%   See also: test_model, print_whiteness
%

numTrialsRejected=sum(h);
numTrials=size(h,1);
numChannels=length(labels);

fprintf('\n\t\t');

for i=1:numChannels
    fprintf('\t%s',labels{i});
end

fprintf('\nRejected: ')

for i=1:numChannels
    fprintf('\t%4d/%-6d',numTrialsRejected(i),numTrials);
end

fprintf('\n')
    
end