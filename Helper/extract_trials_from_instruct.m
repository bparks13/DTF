function [ind_curr_start,ind_curr_end,numTrials]=extract_trials_from_instruct(instruct,condition)
%% [ind_curr_start,ind_curr_end,numTrials]=extract_trials_from_instruct(instruct)
%
%  Helper function to retrieve the start/end indices of the current condition
%
%   Inputs:
%    - instruct: Vector containing integer values defining which condition is currently
%       occurring during the trial. Type is single
%    - condition: Integer value defining the state matching a specific condition
%
%   Outputs:
%    - ind_curr_start: All indices where the current condition starts
%    - ind_curr_end: All indices where the current condition ends
%    - numTrials: Number of trials found
%
%  See also: force_instruct_to_zero, parse_postacq_type, extract_visual_stim_for_intraop,
%  extract_visual_stim_for_closed_loop, assign_task_values_to_instruct,
%  remove_brief_epochs
%

ind_change_all=find(diff(instruct) ~= 0)+1;
ind_change_curr=find(diff(instruct) ~= 0 & instruct(2:end) == condition)+1;

numTrials=length(ind_change_curr);

ind_curr_start=zeros(numTrials,1);
ind_curr_end=zeros(numTrials,1);

currTrial=1;

for i=1:length(ind_change_all)
    if instruct(ind_change_all(i)) == condition
        if i ~= length(ind_change_all)
            ind_curr_start(currTrial)=ind_change_all(i);
            ind_curr_end(currTrial)=ind_change_all(i+1);
        end
        currTrial=currTrial+1;
    end
end

if ind_curr_start(end) == 0
    ind_curr_start(end)=[];
    ind_curr_end(end)=[];
    numTrials=numTrials-1;
end


end