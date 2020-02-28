function [newFile,subjID_simple,dateID_simple,runID_simple]=simplify_filename(subjID,dateID,runID,addOn)
%% [newFile,subjID_simple,dateID_simple,runID_simple]=simplify_filename(subjID,dateID,runID,addOn)
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
%    - subjID_simple: Simplified subject ID, in the format SXX, where XX is the ID number
%    - dateID_simple: Simplified date, in the format DXX, where XX is the date value
%    - runID_simple: Simplified run ID, in the format RXX, where XX is the run number
%
%  See also: processing_pipeline
%

if nargin==3
    addOn='';
end

subjID_simple='';
dateID_simple='';
runID_simple='';

%  TEMPLATE 
%
% elseif strcmp(subjID,'')
%     subjID_simp='Sxx';
%     
%     if strcmp(dateID,'')
%         dateID_simp='Dxx';
%         
%         if strcmp(runID,'')
%             runID_simp='Rxx';
%         end
%     end

if strcmp(subjID,'ET_CL_002')
    subjID_simple='S01';
    
    if strcmp(dateID,'2018_02_01')
        dateID_simple='D01';
        
        if strcmp(runID,'run9')
            runID_simple='R01';
        end
    end
elseif strcmp(subjID,'ET_CL_004')
    subjID_simple='S02';
    
    if strcmp(dateID,'2018_06_20')
        dateID_simple='D01';
        
        if strcmp(runID,'run5')
            runID_simple='R01';
        end
    elseif strcmp(dateID,'2018_08_23')
        dateID_simple='D02';
        
        if strcmp(runID,'run1')
            runID_simple='R01';
        end
    end
elseif strcmp(subjID,'ET_OR_STIM_018')
    subjID_simple='S03';
    
    if strcmp(dateID,'2018_11_28')
        dateID_simple='D01';
        
        if strcmp(runID,'run12')
            runID_simple='R01';
        end
    end
elseif strcmp(subjID,'TS04 Double DBS Implantation')
    subjID_simple='S04';
    
    if strcmp(dateID,'2017_03_01')
        dateID_simple='D01';
        
        if strcmp(runID,'run16')
            runID_simple='R01';
        end
    end
end

if isempty(subjID_simple) || isempty(dateID_simple) || isempty(runID_simple)
    error('Invalid inputs given. Check that all inputs match expected values. Check that these inputs exist in the current context');
end

counter=1;

newFile=fullfile(get_root_path,'Files',sprintf('%s_%s_%s%s.mat',subjID_simple,dateID_simple,runID_simple,addOn));

while exist(newFile,'file') == 2
    newFile=fullfile(get_root_path,'Files',sprintf('%s_%s_%s%s_(%d).mat',subjID_simple,dateID_simple,runID_simple,addOn,counter));
    counter=counter+1;
end

end