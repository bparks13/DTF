function plot_surrogate_histogram(distribution)
%% plot_surrogate_histogram(distribution)
%
%  For surrogate analysis that is done between conditions (i.e. method == 'combine'), the
%  distribution returned is a cell array. This function puts the distribution into a
%  format that histogram likes, and then plots it
%
%   Inputs:
%    - distribution: Either a cell array or a struct containing cell arrays for each
%       condition (each condition is assumed to be the same cell array, and only the first
%       condition will be plotted). Size of the cell array is [n x c], where n is the
%       number of iterations run, and c is the number of channels
%
%   Outputs:
%    Figure containing a histogram of the distribution of trials that were randomly chosen
%       for all the analyses
%
%  See also: surrogate_analysis
%

if isstruct(distribution)
    fields=fieldnames(distribution);
    distribution=distribution.(fields{1});
end

C=categorical(distribution(:));
figure; histogram(C);

% trialNames=unique(distribution);
% numTrials=length(trialNames);

% trialCount=nan(numTrials,1);

% for i=1:numTrials
%     trialCount(i)=sum(sum(strcmp(distribution,trialNames{i})));
% end

end