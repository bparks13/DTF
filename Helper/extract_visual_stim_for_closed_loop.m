function instruct=extract_visual_stim_for_closed_loop(file)
%% instruct=extract_visual_stim_for_closed_loop(file)
%
%  For closed-loop recordings that have visual_cues from the DAC, extract the preprocessed
%  stim cues from the 'process' folder on the server, with additional processing to ensure
%  accuracy
%
%   Inputs
%    - file: Filename specifying the full file path to the file
%
%   Outputs
%    - instruct: Vector containing integer values defining which condition is currently
%       occurring during the trial
%
%  See also: load_data, load_channels_labels_conditions
%

prepath='\\gunduz-lab.bme.ufl.edu\Study_ET_Closed_Loop\process';

[newPath,run_ID]=fileparts(file);
[newPath]=fileparts(newPath);
[newPath,date_ID]=fileparts(newPath);
[~,subj_ID]=fileparts(newPath);

if strcmp(subj_ID,'ET_CL_004')
    subj_ID='ET04';
elseif strcmp(subj_ID,'ET_CL_002')
    subj_ID='ET02';
else
    error('Undefined subject');
end

postacqFileName=fullfile(prepath,subj_ID,date_ID,sprintf('postacq_%s.mat',run_ID));

data=load(postacqFileName);

events=data.datastorage_postacq.evt;
instruct=zeros(length(events.right_hand_go),1);

% Here, right logic for when GO and CUE are separate vs. overlap, right vs. left, length
% of the condition, etc.

end