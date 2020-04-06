function [field,value]=parse_postacq_type(postacq_type)
%% [field,value]=parse_postacq_type(postacq_type)
%
%  Figure out which postacquisition fields are needed to create the instruct variable, and
%  the value of that particular field, based on the string postacq_type
%
%   Inputs:
%    - postacq_type: String defining the types of tasks run, as well as the corresponding
%       value of that task for creating an instruct variable not based on cues 
%
%   Outputs:
%    - field: Cell array defining the name of the field in datastorage_postacq.evt that
%       corresponds to postacq_type task in question
%    - value: Value of the field at the same index in the instruct variable
%
%  See also: load_variables, extract_visual_stim_for_intraop, load_data,
%  extract_visual_stim_for_closed_loop
%

tasks=regexp(postacq_type,'[_]','split');

field=cell(length(tasks),1);
value=1:length(tasks);

for i=1:length(tasks)
    switch tasks{i}
        case 'R'
            field{i}='Rest';
            
        case 'CR'
            field{i}='right_hand_cue';
            
        case 'CL'
            field{i}='left_hand_cue';
            
        case 'MR'
            field{i}='right_hand_go';
            
        case 'ML'
            field{i}='left_hand_go';
            
        case 'FR'
            field{i}='right_feel_go';
%             warning('Check if there is a feel_gocue here')
            
        case 'FL'
            field{i}='left_feel_go';
%             warning('Check if there is a feel_gocue here')
            
        case 'CuR'
            field{i}='right_hand_go';   % Cup Right
            warning('Check if this makes sense from the video')
            
        case 'CuL'
            field{i}='left_hand_go';    % Cup Left
            warning('Check if this makes sense from the video')
            
        case 'IR'
            field{i}='right_hand_gocue';% Imagine Right
            warning('Check if this makes sense')
            
        case 'IL'
            field{i}='left_hand_gocue'; % Imagine Left
            warning('Check if this makes sense')
            
        otherwise
            error('Invalid postacq_type identifier');
    end
end

end