function instruct=assign_task_values_to_instruct(numSamples,fs,field,value,evt)
%% instruct=assign_task_values_to_instruct(numSamples,fs,field,value,evt)
%
%  Helper function to assign the value of each 'field' to instruct. 
%
%   Inputs:
%    - numSamples: Number of samples (length) to make the instruct variable
%    - fs: Sampling frequency. Only needed if Rest is a given field, but is a required
%       input
%    - field: Cell array defining the name of the field in datastorage_postacq.evt that
%       corresponds to postacq_type task in question
%    - value: Value of the field at the same index in the instruct variable
%    - evt: Event subfield from the data_postacq struct
%
%   Outputs:
%    - instruct: Vector containing integer values defining which condition is currently
%       occurring during the trial. Type is single
%
%  See also: force_instruct_to_zero, parse_postacq_type, extract_visual_stim_for_intraop,
%  extract_visual_stim_for_closed_loop
%

instruct=zeros(numSamples,1,'single');

for i=1:length(field)
    if contains(field{i},'Rest')
        indStartRest=9*fs;
        instruct(indStartRest:end-1)=value(contains(field,'Rest'));
    else
        instruct(evt.(field{i}))=value(i);
    end
end

end