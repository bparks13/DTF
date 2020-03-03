function names=get_structure_names(subjID)
%% names=get_structure_names(subjID)
%
%  Given a particular subject, return a cell array containing the names of the bipolar
%  combinations of channels.
%
%   Inputs:
%    - subjID: Subject ID, can be either the full title (i.e. ET_CL_XXX) or the internally
%       defined subject ID 
%
%   Outputs:
%    - names: Cell array of character vectors that correspond to the cortical/subcortical
%       structures instead of the bipolar channels
%
%  See also: load_variables.m
%

if strcmp(subjID,'ET_CL_002') || strcmp(subjID,'S01')
    names={'L Vim','L PSA','L M1/S1','L S1/S2'};
elseif strcmp(subjID,'ET_CL_004') || strcmp(subjID,'S02')
    names={'L Thal','L Vim','L M1/S1','L S2'};
elseif strcmp(subjID,'ET_OR_STIM_018') || strcmp(subjID,'S03')
    names={'L VOp','L VOp/PSA','L Vim','L PSA','L M1','L M1/S1','L S1/S2'};
elseif strcmp(subjID,'ET_CL_001') || strcmp(subjID,'S05')
    names={'R Vim','R PSA','R M1/S1','R S1/S2'};
else
    names={};
end

end