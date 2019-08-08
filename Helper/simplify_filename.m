function [newFile]=simplify_filename(subjID,dateID,runID,addOn)
%% [newFile]=simplify_filename(subjID,dateID,runID,addOn)
%
%  Instead of using the entire subject ID/date of recording/run ID, simplify the new
%  filename by assigning a number to each input; i.e. subjID == 'ET_CL_04' becomes S02,
%  dateID == '2018_06_20' becomes D01, etc.
%
%   Inputs:
%    - subjID: Subject ID as a string
%    - dateID: Date of recording as a string
%    - runID: Run number for the particular date given as a string
%    - addOn: Optional input, contains any additional information needed in the filename
%
%   Outputs:
%    - newFile: String containing the absolute path to the folder that contains all of the
%       processed files; guaranteed to be unique
%
%  See also: processing_pipeline
%

if nargin==3
    addOn='';
end

subjID_simp='';
dateID_simp='';
runID_simp='';

if strcmp(subjID,'ET_CL_001')
    subjID_simp='S01';
    
    if strcmp(dateID,'2017_05_17')
        dateID_simp='D01';
        
        if strcmp(runID,'run12')
            runID_simp='R01';
        end
    end
elseif strcmp(subjID,'ET_CL_004')
    subjID_simp='S02';
    
    if strcmp(dateID,'2018_06_20')
        dateID_simp='D01';
        
        if strcmp(runID,'run5')
            runID_simp='R01';
        end
    end
end

if isempty(subjID_simp) || isempty(dateID_simp) || isempty(runID_simp)
    error('Invalid inputs given. Check that all inputs match expected values. Check that these inputs exist in the current context');
end

counter=1;

newFile=fullfile(get_root_path,'Files',sprintf('%s_%s_%s%s.mat',subjID_simp,dateID_simp,runID_simp,addOn));

while exist(newFile,'file') == 2
    newFile=fullfile(get_root_path,'Files',sprintf('%s_%s_%s%s_(%d).mat',subjID_simp,dateID_simp,runID_simp,addOn,counter));
    counter=counter+1;
end

end