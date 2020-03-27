function plot_postacq_variables(data)
%% plot_postacq_variables(data)
%
%  Given the data struct containing the evt field, plot all events that are listed as
%  fields in 'evt'
%
%   Inputs:
%    - data: This should be the datastorage_postacq struct, containing the following
%       subfields
%    -- base_path: Irrelevant for me
%    -- run_name: Irrelevant for me
%    -- run_i: Possibly relevant to check that runID is the same for files?
%    -- evt: Struct containing an assortment of fields pertaining to particular events
%        (such as right_hand_go, etc.)
%
%   Outputs:
%    Figure containing all subplots of each individual event, with a value of 0 or 1
%     depending on when an event occurs. Check for overlaps, missing data, etc.
%
%  See also: extract_visual_stim_for_closed_loop, extract_visual_stim_for_intraop
%

if ~isstruct(data)
    if ischar(data)
        data=load(data);
        data=data.datastorage_postacq;
    else
        error('Invalid input');
    end    
end

fields=fieldnames(data);

if ~any(contains(fields,'evt'))
    error('No evt field detected');
end

fields=sort(fieldnames(data.evt));
numFields=length(fields);

figure;
hAx=nan(numFields,1);

for i=1:numFields
    hAx(i)=subplot_tight(numFields,1,i);
    
    plot(data.evt.(fields{i})); title(fields{i},'Interpreter','none');
end

linkaxes(hAx,'x');

end