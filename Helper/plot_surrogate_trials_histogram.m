function plot_surrogate_trials_histogram(distribution,condition)
%% plot_surrogate_trials_histogram(distribution,condition)
%
%  For surrogate analysis that is done between conditions (i.e. method == 'combine'), the
%  distribution returned is a cell array. This function puts the distribution into a
%  format that histogram likes, and then plots it
%
%   Inputs:
%    - distribution: Either a cell array or a struct containing cell arrays for each
%       condition. Size of the cell array is [n x c], where n is the
%       number of iterations run, and c is the number of channels
%    - condition: String defining the condition of the current surrogate; used for the
%       title
%
%   Outputs:
%    Figure containing a histogram of the distribution of trials that were randomly chosen
%       for all the analyses
%
%  See also: surrogate_analysis
%

if isstruct(distribution)
    fields=fieldnames(distribution);
    
    for i=1:length(fields)
        plot_surrogate_trials_histogram(distribution.(fields{i}),fields{i});
    end
    
    return
end

C=categorical(distribution(:));
figure; histogram(C);
xlabel('Trial #'); ylabel('# of Times Used');
title(sprintf('Trials Used for Surrogate Analysis - %s',condition));

end