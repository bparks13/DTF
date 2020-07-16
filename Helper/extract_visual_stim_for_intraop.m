function instruct=extract_visual_stim_for_intraop(data,data_postacq,cues_only,postacq_type,minLength)
%% instruct=extract_visual_stim_for_intraop(data,data_postacq,cues_only,postacq_type,minLength)
%
%  Helper function to grab the correct form of the visual_stim variable, depending on
%  whether or not it should be cues based or acceleration based.
%
%   Inputs: 
%    - data: Struct containing all relevant data
%    -- src: Struct containing all data from all sources
%    --- visual_stim: Struct from all subjects prior to ET_CL_06, containing the cues
%    --- LFP: Struct with all LFP data
%    ---- state: Struct from all subjects after and including ET_CL_06, containing the
%          cues
%    - data_postacq: Struct containing the datastorage_postacq struct,
%       for aligning data according to acceleration
%    - cues_only: Boolean defining whether or not to use the trials based on acceleration
%       data (false) or to go based on the cues only (true, default)
%    - postacq_type: String defining the types of cues given in this run
%    - minLength: Minimum length of a realization in samples
%
%   Outputs
%    - instruct: Vector containing integer values defining which condition is currently
%       occurring during the trial. Type is uint16
%
%  See also: load_data, load_variables, extract_visual_stim_for_closed_loop
%

if cues_only
    if isfield(data.src,'visual_stim')
        instruct=data.src.visual_stim.data;
    elseif isfield(data.src.LFP,'state')
        instruct=data.src.LFP.state;
    end
    
    return
end

[field,value]=parse_postacq_type(postacq_type);

% Rest is not given in the postacq file; therefore i will define it as any period not
% overlapping with a movement/cue, and not within 9 seconds of the beginning of the run

evt_fields=fieldnames(data_postacq.evt);
numSamples=length(data_postacq.evt.(evt_fields{1}));
fs=extract_sampling_frequency(data);

instruct=assign_task_values_to_instruct(numSamples,fs,field,value,data_postacq.evt);

instruct=force_instruct_to_zero(instruct);

instruct=remove_brief_epochs(instruct,minLength);

instruct=uint16(instruct);

end