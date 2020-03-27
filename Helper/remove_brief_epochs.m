function instruct=remove_brief_epochs(instruct,minLength)
%% instruct=remove_brief_epochs(instruct,minLength)
%
%  Helper function to find any trials of any task that are shorter than the minimum length
%  given (which is the realization length), and remove them so as not to confound further
%  processing steps.
%
%   Inputs:
%    - instruct: Vector containing integer values defining which condition is currently
%       occurring during the trial. Type is single
%    - minLength: Minimum length of a realization in samples
%
%   Outputs:
%    - instruct: Vector containing integer values defining which condition is currently
%       occurring during the trial. Type is single. Brief epochs removed
%
%  See also: force_instruct_to_zero, parse_postacq_type, extract_visual_stim_for_intraop,
%  extract_visual_stim_for_closed_loop, assign_task_values_to_instruct
%

values=unique(instruct);
values(values==0)=[];

for i=1:length(values)
    [ind_curr_start,ind_curr_end,~]=extract_trials_from_instruct(instruct,values(i));
    
    trial_lengths=ind_curr_end-ind_curr_start;
    
    too_brief=trial_lengths < minLength;
    
    if any(too_brief)
        for j=1:length(too_brief)
            if too_brief(j)
                instruct(ind_curr_start(j):ind_curr_end(j))=0;
            end
        end
    end
end

end