function instruct=extract_visual_stim_for_closed_loop(data_postacq,cues_only)
%% instruct=extract_visual_stim_for_closed_loop(data_postacq,cues_only)
%
%  For closed-loop recordings that have visual_cues from the DAC, extract the preprocessed
%  stim cues from the 'process' folder on the server, with additional processing to ensure
%  accuracy
%
%   Inputs
%    - data_postacq: Struct containing the datastorage_postacq field, which houses all
%       relevant data for time-aligning, and the postacq_type string, defining which
%       condition is which value
%    - cues_only: Boolean defining whether or not to use the trials based on acceleration
%       data (false) or to go based on the cues only (true, default)
%
%   Outputs
%    - instruct: Vector containing integer values defining which condition is currently
%       occurring during the trial
%
%  See also: load_data, load_variables, extract_visual_stim_for_intraop
%

error('need to rewrite the logic here to utilize the postacq_type variable, and the postacq struct');

%% Definitions

prepath='\\gunduz-lab.bme.ufl.edu\Study_ET_Closed_Loop\process';

REST=1;
R_CUE=2;
R_GO=4;
L_CUE=3;
L_GO=5;

if nargin == 1
    cues_only=true;
end

%% Load postacq file

[newPath,run_ID]=fileparts(data_postacq);
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

data_postacq=load(postacqFileName);

events=data_postacq.datastorage_postacq.evt;

%% Create the instruct variable and set the values

instruct=zeros(length(events.right_hand_go),1)+REST;

% Here, write logic for when GO and CUE are separate vs. overlap, right vs. left, length
% of the condition, etc.

if cues_only
    instruct(events.right_hand_cue)=R_CUE;
    instruct(events.right_hand_gocue)=R_GO;
    instruct=fill_instruct(instruct,events.right_hand_cue,events.right_hand_gocue);
    instruct(events.left_hand_cue)=L_CUE;
    instruct(events.left_hand_gocue)=L_GO;   
    instruct=fill_instruct(instruct,events.left_hand_cue,events.left_hand_gocue);
else
    instruct(events.right_hand_cue)=R_CUE;
    instruct(events.right_hand_go)=R_GO;
    instruct=clean_instruct(instruct);
%     instruct=fill_instruct(instruct,events.right_hand_cue,events.right_hand_go);
    instruct(events.left_hand_cue)=L_CUE;
    instruct(events.left_hand_go)=L_GO;
%     instruct=fill_instruct(instruct,events.left_hand_cue,events.left_hand_go);
end

%% Internal functions
    function instruct=fill_instruct(instruct,cue1,cue2)
    %% Internal function to interpolate from the end of one cue to the beginning of another cue
    
    cue1_end=find(diff(cue1)==-1)+1;
    cue2_begin=find(diff(cue2)==1)+1;
    
    for i=1:length(cue1_end)
        instruct(cue1_end(i):cue2_begin(i))=instruct(cue1_end(i)-1);
    end
        
    end

    function instruct=clean_instruct(instruct)
    %% Internal function to clean the instruct variable of miscellaneous movements
        
    error('have not created this function yet');
    end

end