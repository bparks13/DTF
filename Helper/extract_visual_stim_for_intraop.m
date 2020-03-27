function instruct=extract_visual_stim_for_intraop(data,cues_only,data_postacq,minLength)
%% instruct=extract_visual_stim_for_intraop(data,cues_only,data_postacq,minLength)
%
%  Helper function to grab the correct form of the visual_stim variable, depending on
%  whether or not it should be cues based or acceleration based.
%
%   Inputs: 
%    - data: Struct containing the datastorage field, which houses all relevant data
%    - cues_only: Boolean defining whether or not to use the trials based on acceleration
%       data (false) or to go based on the cues only (true, default)
%    - data_postacq: Struct containing two fields; one is the datastorage_postacq struct,
%       for aligning data according to acceleration, and the other is the postacq_type
%       string defining the types of tasks run, as well as the corresponding value of that
%       task for creating an instruct variable not based on cues 
%    - minLength: Minimum length of a realization in samples
%
%   Outputs
%    - instruct: Vector containing integer values defining which condition is currently
%       occurring during the trial. Type is uint16
%
%  See also: load_data, load_variables, extract_visual_stim_for_closed_loop
%

if cues_only
    instruct=data.datastorage.src.visual_stim.data;
    return
end

[field,value]=parse_postacq_type(data_postacq.postacq_type);

% Rest is not given in the postacq file; therefore i will define it as any period not
% overlapping with a movement/cue, and not within 9 seconds of the beginning of the run

evt_fields=fieldnames(data_postacq.datastorage_postacq.evt);
numSamples=length(data_postacq.datastorage_postacq.evt.(evt_fields{1}));
fs=extract_sampling_frequency(data);

instruct=assign_task_values_to_instruct(numSamples,fs,field,value,data_postacq.datastorage_postacq.evt);

instruct=force_instruct_to_zero(instruct);

instruct=remove_brief_epochs(instruct,minLength);

instruct=uint16(instruct);

end