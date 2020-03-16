function [channels,labels,conditions,cond_labels,visit_type]=load_variables(patient,record_date,run_id,config)
%% [channels,labels,conditions,cond_labels,visit_type]=load_variables(patient,record_date,run_id,config)
%
%  Given the relevant information, return the channels, labels, and conditions
%  corresponding to those channels. Additional options can be given using 'config'.
%
%   Inputs:
%    - patient: String containing the patient identifier
%    - record_date: String containing the recording date in yyyy_mm_dd
%    - run_id: String containing the run ID, without the file extension appended
%    - config: Optional struct input that allows for additional parameters to be given
%    -- default: Bool denoting whether to use the predetermined set of channels and labels
%        or not. Default value is true; if true, overrides all following fields even if
%        they exist. Must be defined as false for any preset/custom inputs to be used
%    -- preset: Int defining a preset channels/labels combination. Will have a different
%        meaning depending on the specific patient/date/run combination used. If this
%        is defined it will override the following fields even if they exist. 
%        Example may be for ET_CL_004/2018_06_20/run5, where '1' == bipolar channels
%        Vim (1-0) and Cort (3-2).
%    -- custom: Matrix/vector defining the specific channels to use, and then returns the
%        corresponding labels. Vector should contain monopolar combinations requested
%        (not recommended), while a matrix [c x 2] contains the bipolar combinations
%        where the channel in the second column is subtracted from the channel in the
%        first column.
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
%    - visit_type: String defining whether the run is 'intraop' or 'closed-loop'
%
%  See also: processing_pipeline
%

if nargin == 3
    bool_default=true;
    preset_value=[];
    bool_custom=false;
else
    if isstruct(config)
        if isfield(config,'default') && ~isempty(config.default)
            bool_default=config.default;
        else
            bool_default=true;
        end
        
        if isfield(config,'preset') && ~isempty(config.preset)
            preset_value=config.preset;
            bool_default=false;
        else
            preset_value=[];
        end
        
        if isfield(config,'custom') && ~isempty(config.custom)
            bool_custom=any(config.custom);
            channels=config.custom; %#ok<NASGU>
        else
            bool_custom=false;
        end
    end
end

% Template
% elseif strcmp(patient,'')
%     if strcmp(record_date,'')
%         if strcmp(run_id,'')
%             channels=[];
%             labels={};
%             conditions=[];
%             cond_labels={};
%             visit_type='';
%             return
%         else
%             error('WARNING: Invalid run number given. Please set channels/labels for this combination');
%         end
%     else
%         error('WARNING: Invalid recording date given. Please set channels/labels for this combination');
%     end

if bool_default
    % set all default parameters here for patients/dates/runs
    if strcmp(patient,'ET_CL_002') % S01
        if strcmp(record_date,'2018_02_01')
            visit_type='intraop';
            channels=[8,7;6,5;4,3;2,1];
            labels={'Vim (3-2)','Vim (1-0)','Cort (3-2)','Cort (1-0)'};
            
            if strcmp(run_id,'run9')
                conditions=[1,2,3,4];
                cond_labels={'Rest','MoveRight','MoveLeft','FeelRight'};
                return
            elseif strcmp(run_id,'run10')
                conditions=[1,2,3];
                cond_labels={'Rest','CupReachRight','CupReachLeft'};
                return
            elseif strcmp(run_id,'run11')
                conditions=[1,2,3];
                cond_labels={'Rest','MoveRight','MoveLeft'};
                return
            elseif strcmp(run_id,'run12')
                conditions=[1,2,3,4,5];
                cond_labels={'Rest','MoveRight','MoveLeft','ImagineRight','ImagineLeft'};
                return
            else
                error('WARNING: Invalid run number given. Please set channels/labels for this combination');
            end
        else
            error('WARNING: Invalid recording date given. Please set channels/labels for this combination');
        end
    elseif strcmp(patient,'ET_CL_004') % S02
        if strcmp(record_date,'2018_06_20')
            visit_type='intraop';
            channels=[8,7;6,5;4,3;2,1];
            labels={'Vim (3-2)','Vim (1-0)','Cort (3-2)','Cort (1-0)'};
                
            if strcmp(run_id,'run5') || strcmp(run_id,'run5_fs600') || strcmp(run_id,'run5_fs300')
                conditions=[1,2,3,4,5];
                cond_labels={'Rest','CueRight','CueLeft','MoveRight','MoveLeft'};
                return
            else
                error('WARNING: Invalid run number given. Please set channels/labels for this combination');
            end
        elseif strcmp(record_date,'2018_07_26')
            if strcmp(run_id,'run3') || strcmp(run_id,'run4')
                channels=[1;2];
                labels={'E1-E3','E9-E11'};
                conditions=[];
                cond_labels={};
                visit_type='device-implantation';
                warning('Device implantation surgery. Skip this data');
                return             
            end
        elseif strcmp(record_date,'2018_08_23')
            visit_type='closed-loop';
            channels=[1;2];
            labels={'E1-E3','E9-E11'};
            
            if strcmp(run_id,'run1')
                conditions=[1,2,3,4,5];
                cond_labels={'Rest','CueRight','MoveRight','CueLeft','MoveLeft'};
                return
            end
        else
            error('WARNING: Invalid recording date given. Please set channels/labels for this combination');
        end
    elseif strcmp(patient,'ET_OR_STIM_018') % S03
        if strcmp(record_date,'2018_11_28')
            if strcmp(run_id,'run12')
                channels=[14,13;12,11;10,9;8,7;6,5;4,3;2,1];
                labels={'VO (3-2)','VO (1-0)','Vim (3-2)','Vim (1-0)','Cort (5-4)','Cort (3-2)','Cort (1-0)'};
                conditions=[1,2,3,4];
                cond_labels={'Rest','MoveRight','MoveLeft','Caress'};
                visit_type='intraop';
                return
            else
                error('WARNING: Invalid run number given. Please set channels/labels for this combination');
            end
        else
            error('WARNING: Invalid recording date given. Please set channels/labels for this combination');
        end
    elseif strcmp(patient,'TS04 Double DBS Implantation') % S04
        if strcmp(record_date,'2017_03_01')
            if strcmp(run_id,'run16')
                channels=[16,15;14,13;12,11;10,9;8,7;6,5;4,3;2,1];
                labels={'R Thal (3-2)','R Thal (1-0)','L Thal (3-2)','L Thal (1-0)',...
                    'R Cort (3-2)','R Cort (1-0)','L Cort (3-2)','L Cort (1-0)'};
                conditions=[1,2,3];
                cond_labels={'Rest','MoveRight','MoveLeft'};
                visit_type='intraop';
                return
            else
                error('WARNING: Invalid run number given. Please set channels/labels for this combination');
            end
        else
            error('WARNING: Invalid recording date given. Please set channels/labels for this combination');
        end
    elseif strcmp(patient,'ET_CL_001') % S05
        if strcmp(record_date,'2017_05_17')
            visit_type='intraop';
            channels=[8,7;6,5;4,3;2,1];
            labels={'Vim (3-2)','Vim (1-0)','Cort (3-2)','Cort (1-0)'};
            
            if strcmp(run_id,'run12')
                conditions=[1,2,3,4,5];
                cond_labels={'Rest','CueRight','CueLeft','MoveRight','MoveLeft'};
                return
            elseif strcmp(run_id,'run13')
                conditions=[1,2];
                cond_labels={'CupReachRight','CupReachLeft'};
                return
            elseif strcmp(run_id,'run15')
                conditions=[1,2];
                cond_labels={'Finger_CupToNoseRight','Finger_CupToNoseLeft'};
                return
            elseif strcmp(run_id,'run17')
                conditions=[1,2,3,4,5];
                cond_labels={'Rest','MoveRight','MoveLeft','ImagineRight','ImagineLeft'};
                return
            else
                error('WARNING: Invalid run number given. Please set channels/labels for this combination');
            end
        else
            error('WARNING: Invalid recording date given. Please set channels/labels for this combination');
        end
    else
        error('WARNING: Invalid patient ID given. Please set channels/labels for this combination');
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
                error('WARNING: Invalid run number given. Please set channels/labels for this combination');
            end
        else
            error('WARNING: Invalid recording date given. Please set channels/labels for this combination');
        end
        
    elseif strcmp(patient,'TS04 Double DBS Implantation') % S04
        if strcmp(record_date,'2017_03_01')
            if strcmp(run_id,'run16')
                if preset_value == 1
                    % This is the same as the default, just without the Right PSA (R Thal (1-0))
                    channels=[16,15;12,11;10,9;8,7;6,5;4,3;2,1];
                    labels={'R Thal (3-2)','L Thal (3-2)','L Thal (1-0)',...
                        'R Cort (3-2)','R Cort (1-0)','L Cort (3-2)','L Cort (1-0)'};
                    conditions=[1,2,3];
                    cond_labels={'Rest','MoveRight','MoveLeft'};
                    visit_type='intraop';
                    return
                end
            else
                error('WARNING: Invalid run number given. Please set channels/labels for this combination');
            end
        else
            error('WARNING: Invalid recording date given. Please set channels/labels for this combination');
        end
    else
        error('WARNING: Invalid patient ID given. Please set channels/labels for this combination');
    end
elseif bool_custom
    % define a way to intake channels themselves
else
    error('WARNING: Invalid set of inputs given in config. No channels/labels returned.')
end

channels=[];
labels=[];
conditions=[];
cond_labels=[];

end