function trialsToDrop=drop_trials(file,condition,numRealizations)
%% trialsToDrop=drop_trials(file,condition,numRealizations)
%
%  Return the indices of trials to drop due to artifact or other problems. If there are no
%  trials to drop, trialsToDrop will be empty
%
%   Inputs:
%    - file: String containing the absolute file path. Not matching any particular
%       portions, but rather the whole name
%    - condition: Single value pertaining to the condition from which to drop the trials
%    - numRealizations: Int defining the number of realizations that each trial was split
%       into. If no number is given, the default value is 1
%
%   Outputs:
%    - trialsToDrop: Vector containing the indices of trials to drop. If there are none to
%       drop, this will be empty
%
%  See also: processing_pipeline, load_data
%

% warning('load_data?drop_trials. ADD A CHECK HERE IF THE TRIALS ARE SEPARATED INTO REALIZATIONS OR NOT');

if nargin==2
    numRealizations=1;
end

trialsToDrop=[];

if strcmp(file,'\\gunduz-lab.bme.ufl.edu\Study_Tourette\TS04 Double DBS Implantation\2017_03_01\preproc\run16')
    if condition==1
        trials=9;
        trialsToDrop=(trials-1)*numRealizations+1:trials*numRealizations;
    end
end

end