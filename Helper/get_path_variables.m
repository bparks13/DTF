function [pre,patID,patID_post,recordDate,runNum]=get_path_variables(subjID,dateID,runID)
%% [pre,patID,patID_post,recordDate,runNum]=get_path_variables(subjID,dateID,runID)
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
%    - pre: String defining the prepath to the raw data. Typically on the server, either
%       the closed-loop folder or the Tourettes folder
%    - patID: String defining the patient ID as given by the folder on the server
%    - patID_post: String defining the patient ID as given by the folder in the 'process'
%       folder (different than the original patient ID for some reason)
%    - recordDate: String defining the date given in the format YYYY_MM_DD, corresponds to
%       the folder on the server for the runs from that day
%    - runNum: String defining the run number given in the format 'runX' or 'runXX', where
%       X is a number defining the particular run corresponding to a single .mat file for
%       the specific run
%

switch subjID
    case 1
        pre='\\gunduz-lab.bme.ufl.edu\\Study_ET_Closed_Loop';
        patID='ET_CL_002';
        patID_post='ET02';
        
        switch dateID
            case 1
                recordDate='2018_02_01';
                
                switch runID
                    case 1
                        runNum='run9';
                        
                    case 2 
                        runNum='run10';
                        
                    case 3 
                        runNum='run11';
                        
                    case 4
                        runNum='run12';
                        
                    otherwise
                        runNum='';
                end
                
            otherwise
                recordDate='';
        end
        
    case 2
        pre='\\gunduz-lab.bme.ufl.edu\\Study_ET_Closed_Loop';
        patID='ET_CL_004';
        patID_post='ET04';
        
        switch dateID
            case 1
                recordDate='2018_06_20';
                
                switch runID
                    case 1
                        runNum='run5';
                        
                    otherwise
                        runNum='';
                end
                
            case 2
                recordDate='2018_08_23';
                
                switch runID
                    case 1 
                        runNum='run1';
                        
                    otherwise
                        runNum='';
                end
                
            otherwise
                recordDate='';
        end
        
    case 3
        pre='\\gunduz-lab.bme.ufl.edu\\Study_ET\\OR\\with_DBS';
        patID='ET_OR_STIM_018';
        patID_post='';
        
        switch dateID
            case 1
                recordDate='2018_11_28';
                
                switch runID
                    case 1
                        runNum='run12';
                        
                    otherwise
                        runNum='';
                end
                
            otherwise
                recordDate='';
        end
        
    case 4
        pre='\\gunduz-lab.bme.ufl.edu\\Study_Tourette';
        patID='TS04 Double DBS Implantation';
        patID_post='';
        
        switch dateID
            case 1
                recordDate='2017_03_01';
                
                switch runID
                    case 1
                        runNum='run16';
                        
                    otherwise
                        runNum='';
                end
                
            otherwise
                recordDate='';
        end
        
    case 5
        pre='\\gunduz-lab.bme.ufl.edu\\Study_ET_Closed_Loop';
        patID='ET_CL_001';
        patID_post='ET01';
        
        switch dateID
            case 1
                recordDate='2017_05_17';
                
                switch runID
                    case 1
                        runNum='run12';
                        
                    case 2
                        runNum='run13';
                        
                    case 3
                        runNum='run15';
                        
                    case 4
                        runNum='run17';
                        
                    otherwise
                        runNum='';
                end
                
            otherwise
                recordDate='';
        end
        
    otherwise
        pre='';
        patID='';
        patID_post='';
end

if isempty(pre) || isempty(patID) || isempty(recordDate) || isempty(runNum)
    error('Some variable was unassigned. Check that the current subject/date/run combination exists here');
end

if isempty(patID_post)
    warning('Postacq patient ID not defined. Ignore if there is no postacq file for this subject');
end

end