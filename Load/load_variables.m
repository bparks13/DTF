function meta=load_variables(meta,config)
%% meta=load_variables(meta,config)
%
%  Given the relevant information, return the channels, labels, and conditions
%  corresponding to those channels. Additional options can be given using 'config'.
%
%   Inputs:
%    - meta: Struct containing relevant information for this run
%    -- path: Struct containing information of the filepath, which includes the
%        information requested in this function
%    --- patient: String containing the patient identifier
%    --- record_date: String containing the recording date in yyyy_mm_dd
%    --- run_num: String containing the run ID, without the file extension appended
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
%    - meta: Same as the input, but with the addition of a vars field
%    -- vars: Struct that contains the variables used in this run
%    --- channels: Depending on the inputs given, is either a vector of monopolar channels,
%         or a matrix of bipolar channels, with size [c x 2], where c is the number of
%         bipolar channels, and the second channel is subtracted from the first channel
%    --- labels: Cell array consisting of strings corresponding to the channels given back.
%         If monopolar channels are given, it lists the structure and number (i.e. Vim0), if
%         bipolar channels are given, it lists the structure and the order of subtraction
%        (i.e. Vim (3-2)).
%    --- conditions: Vector containing all possible conditions available for this specific
%         set of inputs. Size is [1 x m], where m is the number of conditions returned
%    --- cond_labels: Cell array consisting of strings denoting what each condition is
%    --- visit_type: String defining whether the run is 'intraop' or 'closed-loop'
%    --- postacq_type: String used as a regular expression of the order/numbering of conditions,
%         where the string R_CR_CL_MR_ML means there is Rest, Cue Right/Left, Move
%         Right/Left, and the corresponding condition values are 1,2,3,4,5, respectively. 
%
%  See also: processing_pipeline, extract_visual_stim_for_intraop,
%   extract_visual_stim_for_closed_loop, simplify_filename, get_structure_names
%

if nargin == 1
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
%{
elseif strcmp(patient,'')
    if strcmp(record_date,'')
        visit_type='';
        channels=[];
        labels={};

        if strcmp(run_num,'')
            conditions=[];
            cond_labels={};
            postacq_type='';
            return
        else
            error('WARNING: Invalid run number given. Please set channels/labels for this combination');
        end
    else
        error('WARNING: Invalid recording date given. Please set channels/labels for this combination');
    end
%}

if bool_default
    % set all default parameters here for patients/dates/runs
    if strcmp(meta.path.patID,'ET_CL_002') % S01
        % For S01, Right is Contralateral, Left is ipsilateral
        if strcmp(meta.path.record_date,'2018_02_01') % D01
            meta.vars.visit_type='intraop';
            meta.vars.channels=[8,7;6,5;4,3;2,1];
            meta.vars.labels={'Vim (3-2)','Vim (1-0)','Cort (3-2)','Cort (1-0)'};
            
            if strcmp(meta.path.run_num,'run9') % R01
                meta.vars.conditions=[1,2,3,4];
                meta.vars.cond_labels={'Rest','MoveContra','MoveIpsi','FeelContra'};
                meta.vars.postacq_type='R_MR_ML_FR';
                return
            elseif strcmp(meta.path.run_num,'run10') % R02
                meta.vars.conditions=[1,2,3];
                meta.vars.cond_labels={'Rest','CupReachContra','CupReachIpsi'};
                meta.vars.postacq_type='R_CuR_CuL';
                return
            elseif strcmp(meta.path.run_num,'run11') % R03
                meta.vars.conditions=[1,2,3];
                meta.vars.cond_labels={'Rest','MoveContra','MoveIpsi'};
                meta.vars.postacq_type='R_MR_ML';
                return
            elseif strcmp(meta.path.run_num,'run12') % R04
                meta.vars.conditions=[1,2,3,4,5];
                meta.vars.cond_labels={'Rest','MoveContra','MoveIpsi','ImagineContra','ImagineIpsi'};
                meta.vars.postacq_type='R_MR_ML_IR_IL';
                return
            else
                error('WARNING: Invalid run number given. Please set channels/labels for this combination');
            end
        else
            error('WARNING: Invalid recording date given. Please set channels/labels for this combination');
        end
    elseif strcmp(meta.path.patID,'ET_CL_004') % S02
        % For S02, Right is Contra, Left is Ipsi
        if strcmp(meta.path.record_date,'2018_06_20') % D01
            meta.vars.visit_type='intraop';
            meta.vars.channels=[8,7;6,5;4,3;2,1];
            meta.vars.labels={'Vim (3-2)','Vim (1-0)','Cort (3-2)','Cort (1-0)'};
                
            if strcmp(meta.path.run_num,'run5') % R01
                meta.vars.conditions=[1,2,3,4,5];
                meta.vars.cond_labels={'Rest','CueContra','CueIpsi','MoveContra','MoveIpsi'};
                meta.vars.postacq_type='R_CR_CL_MR_ML';
                return
            else
                error('WARNING: Invalid run number given. Please set channels/labels for this combination');
            end
        elseif strcmp(meta.path.record_date,'2018_07_26') % NOT LISTED
            if strcmp(meta.path.run_num,'run3') || strcmp(meta.path.run_num,'run4')
                meta.vars.channels=[1;2];
                meta.vars.labels={'E1-E3','E9-E11'};
                meta.vars.conditions=[];
                meta.vars.cond_labels={};
                meta.vars.visit_type='device-implantation';
                meta.vars.postacq_type='';
                warning('Device implantation surgery. Skip this data');
                return            
            else
                error('WARNING: Invalid run number given. Please set channels/labels for this combination'); 
            end
        elseif strcmp(meta.path.record_date,'2018_08_23') % D02
            meta.vars.visit_type='closed-loop'; 
            meta.vars.channels=[1;2];
            meta.vars.labels={'E1-E3','E9-E11'};
            
            if strcmp(meta.path.run_num,'run1') % R01
                meta.vars.conditions=[1,2,3,4,5];
                meta.vars.cond_labels={'Rest','CueContra','CueIpsi','MoveContra','MoveIpsi'};
                meta.vars.postacq_type='R_CR_CL_MR_ML';
                return
            else
                error('WARNING: Invalid run number given. Please set channels/labels for this combination');
            end
        else
            error('WARNING: Invalid recording date given. Please set channels/labels for this combination');
        end
    elseif strcmp(meta.path.patID,'ET_OR_STIM_018') % S03
        % For S03, Right is Contra, Left is Ipsi
        if strcmp(meta.path.record_date,'2018_11_28') % D01
            meta.vars.channels=[14,13;12,11;10,9;8,7;6,5;4,3;2,1];
            meta.vars.labels={'VO (3-2)','VO (1-0)','Vim (3-2)','Vim (1-0)','Cort (5-4)','Cort (3-2)','Cort (1-0)'};
            meta.vars.visit_type='intraop';
            
            if strcmp(meta.path.run_num,'run12') % R01
                meta.vars.conditions=[1,2,3,4];
                meta.vars.cond_labels={'Rest','MoveContra','MoveIpsi','FeelContra'};
                meta.vars.postacq_type='R_MR_ML_FR';
                return
            else
                error('WARNING: Invalid run number given. Please set channels/labels for this combination');
            end
        else
            error('WARNING: Invalid recording date given. Please set channels/labels for this combination');
        end
    elseif strcmp(meta.path.patID,'TS04 Double DBS Implantation') % S04
        if strcmp(meta.path.record_date,'2017_03_01') % D01
            meta.vars.channels=[16,15;14,13;12,11;10,9;8,7;6,5;4,3;2,1];
            meta.vars.labels={'R Thal (3-2)','R Thal (1-0)','L Thal (3-2)','L Thal (1-0)',...
                'R Cort (3-2)','R Cort (1-0)','L Cort (3-2)','L Cort (1-0)'};
            meta.vars.visit_type='intraop';
            
            if strcmp(meta.path.run_num,'run16') % R01
                meta.vars.conditions=[1,2,3];
                meta.vars.cond_labels={'Rest','MoveRight','MoveLeft'};
                meta.vars.postacq_type='R_MR_ML';
                return
            else
                error('WARNING: Invalid run number given. Please set channels/labels for this combination');
            end
        else
            error('WARNING: Invalid recording date given. Please set channels/labels for this combination');
        end
    elseif strcmp(meta.path.patID,'ET_CL_001') % S05
        % For S05, Right is Ipsi, Left is Contra
        if strcmp(meta.path.record_date,'2017_05_17') % D01
            meta.vars.visit_type='intraop';
            meta.vars.channels=[8,7;6,5;4,3;2,1];
            meta.vars.labels={'Vim (3-2)','Vim (1-0)','Cort (3-2)','Cort (1-0)'};
            
            if strcmp(meta.path.run_num,'run12') % R01
                meta.vars.conditions=[1,2,3,4,5];
                meta.vars.cond_labels={'Rest','CueIpsi','CueContra','MoveIpsi','MoveContra'};
                meta.vars.postacq_type='R_CR_CL_MR_ML';
                return
            elseif strcmp(meta.path.run_num,'run13') % R02
                meta.vars.conditions=[1,2];
                meta.vars.cond_labels={'CupReachIpsi','CupReachContra'};
                meta.vars.postacq_type='CuR_CuL';
                return
            elseif strcmp(meta.path.run_num,'run15') % R03
                meta.vars.conditions=[1,2];
                meta.vars.cond_labels={'Finger_CupToNoseIpsi','Finger_CupToNoseContra'};
                meta.vars.postacq_type='CuR_CuL';
                return
            elseif strcmp(meta.path.run_num,'run17') % R04
                meta.vars.conditions=[1,2,3,4,5];
                meta.vars.cond_labels={'Rest','MoveIpsi','MoveContra','ImagineIpsi','ImagineRight'};
                meta.vars.postacq_type='R_MR_ML_IR_IL';
                return
            else
                error('WARNING: Invalid run number given. Please set channels/labels for this combination');
            end
        else
            error('WARNING: Invalid recording date given. Please set channels/labels for this combination');
        end
    elseif strcmp(meta.path.patID,'ET_CL_006') % S06
        % For S06, Right is Contra and Left is Ipsi
        if strcmp(meta.path.record_date,'2020_05_21') % D01
            meta.vars.visit_type='intraop';
            meta.vars.channels=[8,7;6,5;4,3;2,1];
            meta.vars.labels={'Vim (3-2)','Vim (1-0)','Cort (3-2)','Cort (1-0)'};
            
            if strcmp(meta.path.run_num,'run3') % R01
                meta.vars.conditions=[1,2,3,4,5];
                meta.vars.cond_labels={'Rest','CueContra','CueIpsi','MoveContra','MoveIpsi'};
                meta.vars.postacq_type='R_CR_CL_MR_ML';
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
%     if strcmp(patient,'ET_CL_004')
%         if strcmp(record_date,'2018_06_20')
%             if strcmp(run_num,'run5')
%                 if preset_value == 1
%                     channels=[8,7;6,5;4,3;2,1];
%                     labels={'Vim (3-2)','Vim (1-0)','Cort (3-2)','Cort (1-0)'};
%                     conditions=1;
%                     cond_labels={'Rest'};
%                     return
%                 elseif preset_value == 2
%                     channels=[8;7;6;5;4;3;2;1];
%                     labels={'Vim 3','Vim 2','Vim 1','Vim 0','Cort 3','Cort 2','Cort 1','Cort 0'};
%                     conditions=[1,2,3,4,5];
%                     cond_labels={'Rest','CueRight','CueLeft','MoveRight','MoveLeft'};
%                     return
%                 elseif preset_value == 3
%                     channels=[5;4];
%                     labels={'Vim 0','Cort 3'};
%                     conditions=1;
%                     cond_labels={'Rest'};
%                     return
%                 end
%             else
%                 error('WARNING: Invalid run number given. Please set channels/labels for this combination');
%             end
%         else
%             error('WARNING: Invalid recording date given. Please set channels/labels for this combination');
%         end
%         
%     elseif strcmp(patient,'TS04 Double DBS Implantation') % S04
%         if strcmp(record_date,'2017_03_01')
%             if strcmp(run_num,'run16')
%                 if preset_value == 1
%                     % This is the same as the default, just without the Right PSA (R Thal (1-0))
%                     channels=[16,15;12,11;10,9;8,7;6,5;4,3;2,1];
%                     labels={'R Thal (3-2)','L Thal (3-2)','L Thal (1-0)',...
%                         'R Cort (3-2)','R Cort (1-0)','L Cort (3-2)','L Cort (1-0)'};
%                     conditions=[1,2,3];
%                     cond_labels={'Rest','MoveRight','MoveLeft'};
%                     visit_type='intraop';
%                     return
%                 end
%             else
%                 error('WARNING: Invalid run number given. Please set channels/labels for this combination');
%             end
%         else
%             error('WARNING: Invalid recording date given. Please set channels/labels for this combination');
%         end
%     else
%         error('WARNING: Invalid patient ID given. Please set channels/labels for this combination');
%     end
elseif bool_custom
    % define a way to intake channels themselves
else
    error('WARNING: Invalid set of inputs given in config. No channels/labels returned.')
end

end