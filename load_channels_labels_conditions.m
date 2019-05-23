function [channels,labels,conditions,cond_labels]=load_channels_labels_conditions(patient,record_date,run_id,config)
%% [channels,labels,conditions,cond_labels]=load_channels_and_labels(patient,record_date,run_id,config)
%
%  Given the relevant information, return the channels and labels corresponding to those
%  channels. Additional options can be given using 'config'.
%
%   Inputs:
%    - patient: String containing the patient identifier
%    - record_date: String containing the recording date in yyyy_mm_dd
%    - run_id: String containing the run ID, without the file extension appended
%    - config: Optional struct input that allows for additional parameters to be given
%       default: Bool denoting whether to use the predetermined set of channels and labels
%           or not. Default value is true; if true, overrides all following fields even if
%           they exist. Must be defined as false for any preset/custom inputs to be used
%       preset: Int defining a preset channels/labels combination. Will have a different
%           meaning depending on the specific patient/date/run combination used. If this
%           is defined it will override the following fields even if they exist. 
%           Example may be for ET_CL_004/2018_06_20/run5, where '1' == bipolar channels
%           Vim (1-0) and Cort (3-2).
%       custom: Matrix/vector defining the specific channels to use, and then returns the
%           corresponding labels. Vector should contain monopolar combinations requested
%           (not recommended), while a matrix [c x 2] contains the bipolar combinations
%           where the channel in the second column is subtracted from the channel in the
%           first column.
%       
%   Outputs:
%    - channels: Depending on the inputs given, is either a vector of monopolar channels,
%       or a matrix of bipolar channels, with size [c x 2], where c is the number of
%       bipolar channels, and the second channel is subtracted from the first channel
%    - labels: Cell array consisting of strings corresponding to the channels given back.
%       If monopolar channels are given, it lists the structure and number (i.e. Vim0), if
%       bipolar channels are given, it lists the structure and the order of subtraction
%       (i.e. Vim (3-2)).
%    - conditions: Vector containing all possible conditions available for this specific
%       set of inputs. Size is [1 x m], where m is the number of conditions returned
%    - cond_labels: Cell array consisting of strings denoting what each condition is
%

if nargin == 3
    bool_default=true;
    preset_value=[];
    bool_custom=false;
else
    if isstruct(config)
        if isfield(config,'default')
            bool_default=config.default;
        else
            bool_default=true;
        end
        
        if isfield(config,'preset')
            preset_value=config.preset;
        else
            preset_value=[];
        end
        
        if isfield(config,'custom')
            bool_custom=any(config.custom);
            channels=config.custom; %#ok<NASGU>
        else
            bool_custom=false;
        end
    end
end

if bool_default
    % set all default parameters here for patients/dates/runs
    if strcmp(patient,'ET_CL_004')
        if strcmp(record_date,'2018_06_20')
            if strcmp(run_id,'run5') || strcmp(run_id,'run5_fs600') || strcmp(run_id,'run5_fs300')
                channels=[8,7;6,5;4,3;2,1];
                labels={'Vim (3-2)','Vim (1-0)','Cort (3-2)','Cort (1-0)'};
                conditions=[1,2,3,4,5];
                cond_labels={'Rest','CueRight','CueLeft','MoveRight','MoveLeft'};
                return
            else
                disp('WARNING: Invalid run number given. Please set channels/labels for this combination');
            end
        else
            disp('WARNING: Invalid recording date given. Please set channels/labels for this combination');
        end
    else
        disp('WARNING: Invalid patient ID given. Please set channels/labels for this combination');
    end
elseif ~isempty(preset_value)
    % create some presets here as I need them; this is not expected to be a very large
    % section, but used for testing specific subsets/combinations of channels
    if strcmp(patient,'ET_CL_004')
        if strcmp(record_date,'2018_06_20')
            if strcmp(run_id,'run5')
                if preset_value == 1
                    channels=[8,7;6,5;4,3;2,1];
                    labels={'Vim (3-2)','Vim (1-0)','Cort (3-2)','Cort (1-0)'};
                    conditions=1;
                    cond_labels={'Rest'};
                    return
                elseif preset_value == 2
                    channels=[8;7;6;5;4;3;2;1];
                    labels={'Vim 3','Vim 2','Vim 1','Vim 0','Cort 3','Cort 2','Cort 1','Cort 0'};
                    conditions=[1,2,3,4,5];
                    cond_labels={'Rest','CueRight','CueLeft','MoveRight','MoveLeft'};
                    return
                elseif preset_value == 3
                    channels=[5;4];
                    labels={'Vim 0','Cort 3'};
                    conditions=1;
                    cond_labels={'Rest'};
                    return
                end
            else
                disp('WARNING: Invalid run number given. Please set channels/labels for this combination');
            end
        else
            disp('WARNING: Invalid recording date given. Please set channels/labels for this combination');
        end
    else
        disp('WARNING: Invalid patient ID given. Please set channels/labels for this combination');
    end
elseif bool_custom
    % define a way to intake channels themselves
else
    disp('WARNING: Invalid set of inputs given in config. No channels/labels returned.')
end

channels=[];
labels=[];
conditions=[];
cond_labels=[];

end