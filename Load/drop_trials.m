function trialsToDrop=drop_trials(file,condition)
%% trialsToDrop=drop_trials(file,condition)
%
%  Return the indices of trials to drop due to artifact or other problems. If there are no
%  trials to drop, trialsToDrop will be empty
%
%   Inputs:
%    - file: String containing the absolute file path. Not matching any particular
%       portions, but rather the whole name
%    - condition: Single value pertaining to the condition from which to drop the trials
%
%   Outputs:
%    - trialsToDrop: Vector containing the indices of trials to drop. If there are none to
%       drop, this will be empty
%

trialsToDrop=[];

if strcmp(file,'\\gunduz-lab.bme.ufl.edu\Study_Tourette\TS04 Double DBS Implantation\2017_03_01\preproc\run16')
    if condition==1
        trialsToDrop=9;
    end
end

end