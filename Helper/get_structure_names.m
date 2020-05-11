function names=get_structure_names(subjID,preset)
%% names=get_structure_names(subjID,preset)
%
%  Given a particular subject, return a cell array containing the names of the bipolar
%  combinations of channels.
%
%   Inputs:
%    - subjID: Subject ID, can be either the full title (i.e. ET_CL_XXX) or the internally
%       defined subject ID 
%    - present: Optional input - integer value defining a preset array of contacts that is
%       not the full array
%
%   Outputs:
%    - names: Cell array of character vectors that correspond to the cortical/subcortical
%       structures instead of the bipolar channels
%
%  See also: load_variables.m
%

if nargin==1
    preset=-1;
end

if strcmp(subjID,'ET_CL_002') || strcmp(subjID,'S01')
    names={'L Vim','L PSA','L M1/S1','L S1/S2'};
elseif strcmp(subjID,'ET_CL_004') || strcmp(subjID,'S02')
    names={'L Thal','L VOp','L M1/S1','L S2'};
elseif strcmp(subjID,'ET_OR_STIM_018') || strcmp(subjID,'S03')
    names={'L VOp','L VOp/PSA','L Vim','L PSA','L M1','L M1/S1','L S1/S2'};
elseif strcmp(subjID,'TS04 Double DBS Implantation') || strcmp(subjID,'S04')
    if preset==-1
        names={'R Vim','R PSA','L Vim','L PSA','R PM/M1','R M1/S1','L M1','L S1/S2'};
    elseif preset==1
        names={'R Vim','L Vim','L PSA','R PM/M1','R M1/S1','L M1','L S1/S2'};
    end
elseif strcmp(subjID,'ET_CL_001') || strcmp(subjID,'S05')
    names={'R VOp','R PSA','R M1/S1','R S1/S2'};
else
    names={};
end

end