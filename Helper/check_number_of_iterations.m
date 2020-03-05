function config=check_number_of_iterations(config,x,minIterations)
%% config=check_number_of_iterations(config,x,minIterations)
%
%   Change the number of iterations for each condition so that each trial is represented a
%   minimum number of iterations
%
%   Inputs:
%    - config: Configuration struct for the surrogate_analysis.m function. If a number of
%       iterations is specified, that will be the minimum number of iterations for the
%       entire analysis, unless there are enough realizations to justify increasing the
%       number of iterations
%    - x: Struct containing the input data for the analysis. Only used to determine the
%       number of realizations in each condition. Size is [n x c x r], where n is the
%       number of samples, c is the number of channels, and r is the number of
%       realizations
%    - minIterations: Defines the minimum number of iterations that a single realization
%       should be represented in the distribution. This is different from the minimum
%       number of iterations in the config struct
%
%   Outputs:
%    - config: Same as the input, but now the iterations will be another struct that
%       contains the conditions, since the number of iterations will be
%       condition-dependent
%

totalIterations=1000;

if isfield(config,'iterations')
    if isstruct(config.iterations)
        warning('Iterations is already defined per condition. Assuming they are correct.');
        return
    end
    
    totalIterations=config.iterations;
end

fields=fieldnames(x);
numFields=length(fields);
numChannels=size(x.(fields{1}),2);

config.iterations=[];

for i=1:numFields
    currRealizations=size(x.(fields{i}),3);
    
    if ((currRealizations * minIterations) / numChannels) < totalIterations
        config.iterations.(fields{i})=totalIterations;
    else
        config.iterations.(fields{i})=round((currRealizations * minIterations) / numChannels);
    end
end

end