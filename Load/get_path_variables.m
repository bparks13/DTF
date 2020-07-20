function meta=get_path_variables(subjID,dateID,runID,addon,notes)
%% meta=get_path_variables(subjID,dateID,runID,addon,notes)
%
%  Given the simplified IDs, returns all the necessary and patient-specific variables to
%  load the current run
%
%   Inputs:
%    - subjID: Integer defining the value of the subject; an input of 1 will
%       correspond to S01 (ET_CL_02)
%    - dateID: Integer defining the value of the recording date for this subject; an input
%       of 1 for S01 corresponds to 2018_02_01
%    - runID: Integer defining the specific run of the given date; an input of 1 for S01
%       on 2018_02_01 corresponds to run9
%
%   Outputs:
%    - meta: Struct containing the following strings defining meta information about the
%       specific run
%    -- path: Struct containing the path variables needed to open the data file
%    --- prepath: String defining the prepath to the raw data. Typically on the server, either
%         the closed-loop folder or the Tourettes folder
%    --- patID: String defining the patient ID as given by the folder on the server
%    --- patID_postacq: String defining the patient ID as given by the folder in the 'process'
%         folder (different than the original patient ID for some reason)
%    --- record_date: String defining the date given in the format YYYY_MM_DD, corresponds to
%         the folder on the server for the runs from that day
%    --- run_num: String defining the run number given in the format 'runX' or 'runXX', where
%         X is a number defining the particular run corresponding to a single .mat file for
%         the specific run
%    --- addon: String containing any additional info to append to the end of the filename
%           *** MIGHT NOT NEED THIS IF I'M SAVING THE FILE IN THE PREPROC FOLDER
%    --- notes: String containing any relevant notes for this run
%

meta=struct;

switch subjID
    case 1
        meta.path.pre='\\gunduz-lab.bme.ufl.edu\\Study_ET_Closed_Loop';
        meta.path.patID='ET_CL_002';
        meta.path.patID_postacq='ET02';
        
        switch dateID
            case 1
                meta.path.record_date='2018_02_01';
                
                switch runID
                    case 1
                        meta.path.run_num='run9';
                        
                    case 2 
                        meta.path.run_num='run10';
                        
                    case 3 
                        meta.path.run_num='run11';
                        
                    case 4
                        meta.path.run_num='run12';
                        
                    otherwise
                        error('Invalid run for this subject/date combination');
                end
                
            otherwise
                error('Invalid date for this subject');
        end
        
    case 2
        meta.path.pre='\\gunduz-lab.bme.ufl.edu\\Study_ET_Closed_Loop';
        meta.path.patID='ET_CL_004';
        meta.path.patID_postacq='ET04';
        
        switch dateID
            case 1
                meta.path.record_date='2018_06_20';
                
                switch runID
                    case 1
                        meta.path.run_num='run5';
                        
                    otherwise
                        error('Invalid run for this subject/date combination');
                end
                
            case 2
                meta.path.record_date='2018_08_23';
                
                switch runID
                    case 1 
                        meta.path.run_num='run1';
                        
                    otherwise
                        error('Invalid run for this subject/date combination');
                end
                
            otherwise
                error('Invalid date for this subject');
        end
        
    case 3
        meta.path.pre='\\gunduz-lab.bme.ufl.edu\\Study_ET\\OR\\with_DBS';
        meta.path.patID='ET_OR_STIM_018';
        meta.path.patID_postacq='';
        
        switch dateID
            case 1
                meta.path.record_date='2018_11_28';
                
                switch runID
                    case 1
                        meta.path.run_num='run12';
                        
                    otherwise
                        error('Invalid run for this subject/date combination');
                end
                
            otherwise
                error('Invalid date for this subject');
        end
        
    case 4
        meta.path.pre='\\gunduz-lab.bme.ufl.edu\\Study_Tourette';
        meta.path.patID='TS04 Double DBS Implantation';
        meta.path.patID_postacq='';
        
        switch dateID
            case 1
                meta.path.record_date='2017_03_01';
                
                switch runID
                    case 1
                        meta.path.run_num='run16';
                        
                    otherwise
                        error('Invalid run for this subject/date combination');
                end
                
            otherwise
                error('Invalid date for this subject');
        end
        
    case 5
        meta.path.pre='\\gunduz-lab.bme.ufl.edu\\Study_ET_Closed_Loop';
        meta.path.patID='ET_CL_001';
        meta.path.patID_postacq='ET01';
        
        switch dateID
            case 1
                meta.path.record_date='2017_05_17';
                
                switch runID
                    case 1
                        meta.path.run_num='run12';
                        
                    case 2
                        meta.path.run_num='run13';
                        
                    case 3
                        meta.path.run_num='run15';
                        
                    case 4
                        meta.path.run_num='run17';
                        
                    otherwise
                        error('Invalid run for this subject/date combination');
                end
                
            otherwise
                error('Invalid date for this subject');
        end
        
    case 6
        meta.path.pre='\\gunduz-lab.bme.ufl.edu\\Study_ET_Closed_Loop';
        meta.path.patID='ET_CL_006';
        meta.path.patID_postacq='ET06';
        
        switch dateID
            case 1
                meta.path.record_date='2020_05_21';
                
                switch runID
                    case 1
                        meta.path.run_num='run3';
                        
                    case 2
                        meta.path.run_num='run4';
                        
                    case 3
                        meta.path.run_num='run7';
                        
                    case 4
                        meta.path.run_num='run8';
                        
                    case 5
                        meta.path.run_num='run9';
                        
                    case 6
                        meta.path.run_num='run10';
                        
                    otherwise
                        error('Invalid run for this subject/date combination');
                end
                
            otherwise
                error('Invalid date for this subject');
        end
        
    otherwise
        error('Invalid subject')
end

if ~isfield(meta.path,'patID_postacq') || isempty(meta.path.patID_postacq)
    warning('Postacq patient ID not defined. Ignore if there is no postacq file for this subject');
end

meta.path.midpath='preproc';
meta.path.midpath_postacq='process';
meta.path.addon=addon;
meta.notes=notes;

end