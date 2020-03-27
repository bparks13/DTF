function instruct=force_instruct_to_zero(instruct)
%% instruct=force_instruct_to_zero(instruct)
%
%  Helper function that forces the instruct variable to zero anytime there is a change in
%  the task. Assumed that this is called after assign_task_values_to_instruct.m. 
%
%   Inputs:
%    - instruct: Vector containing integer values defining which condition is currently
%       occurring during the trial. Type is single. Should only be zero during the first X
%       seconds of the run.
%   
%   Outputs:
%    - instruct: Vector containing integer values defining which condition is currently
%       occurring during the trial. Type is single. Value of zero occurs anytime there is
%       a change in task for easier separation of trials. 
%
%  See also: assign_task_values_to_instruct, parse_postacq_type, extract_visual_stim_for_intraop,
%  extract_visual_stim_for_closed_loop
%

posDiff=find(diff(instruct)>0);
negDiff=find(diff(instruct)<0);

for i=1:length(posDiff)
    instruct(posDiff(i)-1:posDiff(i))=0;
end

for i=1:length(negDiff)
    instruct(negDiff(i):negDiff(i)+1)=0;
end

end