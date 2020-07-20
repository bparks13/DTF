function meta=get_structure_names(meta)
%% meta=get_structure_names(meta)
%
%  Given a particular subject, return a cell array containing the names of the bipolar
%  combinations of channels.
%
%   Inputs:
%    - meta: Struct containing all data relevant to this run. Relevant fields highlighted
%       here
%    -- path: Struct containing details related to the path to the file
%    --- patID: String containing the patient ID
%
%   Outputs:
%    - names: Cell array of character vectors that correspond to the cortical/subcortical
%       structures instead of the bipolar channels
%
%  See also: load_variables.m
%

if strcmp(meta.path.patID,'ET_CL_002') || strcmp(meta.path.patID,'S01')
    meta.vars.contactNames={'L Vim','L PSA','L M1/S1','L S1/S2'};
elseif strcmp(meta.path.patID,'ET_CL_004') || strcmp(meta.path.patID,'S02')
    meta.vars.contactNames={'L Thal','L VOp','L M1/S1','L S2'};
elseif strcmp(meta.path.patID,'ET_OR_STIM_018') || strcmp(meta.path.patID,'S03')
    meta.vars.contactNames={'L VOp','L VOp/PSA','L Vim','L PSA','L M1','L M1/S1','L S1/S2'};
elseif strcmp(meta.path.patID,'TS04 Double DBS Implantation') || strcmp(meta.path.patID,'S04')
    meta.vars.contactNames={'R Vim','R PSA','L Vim','L PSA','R PM/M1','R M1/S1','L M1','L S1/S2'};
elseif strcmp(meta.path.patID,'ET_CL_001') || strcmp(meta.path.patID,'S05')
    meta.vars.contactNames={'R VOp','R PSA','R M1/S1','R S1/S2'};
else
    warning('No contact names are set for this subject. Using the contact labels instead');
    meta.vars.contactNames=meta.vars.labels;
end

end