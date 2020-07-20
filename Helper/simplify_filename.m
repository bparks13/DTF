function [newFile,meta]=simplify_filename(meta,original_file)
%% [newFile,meta]=simplify_filename(meta,original_file)
%
%  Instead of using the entire subject ID/date of recording/run ID, simplify the new
%  filename by assigning a number to each input; i.e. meta.patID == 'ET_CL_04' becomes S02,
%  meta.record_date == '2018_06_20' becomes D01, etc.
%
%   Inputs:
%    - meta: Struct containing all data relevant to this run
%    -- path: Struct containing data pertaining to the location of the original run
%    --- patID: String containing the patient ID
%    --- record_date: String containing the recording date
%    --- run_num: String containing the run number
%    --- addon: String containing anything to add at the end of the file name
%    - original_file: String containing the absolute path to the original file; used to
%       place the new file in the same folder with a different name
%
%   Outputs:
%    - newFile: String containing the absolute path to the new file; guaranteed to be unique
%    - meta: Same as the input, but includes a field defining the IDs
%    -- subject: Struct containing the simplified IDs, separate from the original IDs
%    --- subjID: String defining the subject ID
%    --- dateID: String defining the date ID
%    --- runID: String defining the run ID
%
%  See also: processing_pipeline
%

if strcmp(meta.path.patID,'ET_CL_002')
    meta.subject.subjID='S01';
    
    if strcmp(meta.path.record_date,'2018_02_01')
        meta.subject.dateID='D01';
        
        if strcmp(meta.path.run_num,'run9')
            meta.subject.runID='R01';
        elseif strcmp(meta.path.run_num,'run10')
            meta.subject.runID='R02';
        elseif strcmp(meta.path.run_num,'run11')
            meta.subject.runID='R03';
        elseif strcmp(meta.path.run_num,'run12')
            meta.subject.runID='R04';
        else
            error('Invalid run ID given');
        end
    else
        error('Invalid date ID given');
    end
elseif strcmp(meta.path.patID,'ET_CL_004')
    meta.subject.subjID='S02';
    
    if strcmp(meta.path.record_date,'2018_06_20')
        meta.subject.dateID='D01';
        
        if strcmp(meta.path.run_num,'run5')
            meta.subject.runID='R01';
        else
            error('Invalid run ID given');
        end
    elseif strcmp(meta.path.record_date,'2018_08_23')
        meta.subject.dateID='D02';
        
        if strcmp(meta.path.run_num,'run1')
            meta.subject.runID='R01';
        else
            error('Invalid run ID given');
        end
    else
        error('Invalid date ID given');
    end
elseif strcmp(meta.path.patID,'ET_OR_STIM_018')
    meta.subject.subjID='S03';
    
    if strcmp(meta.path.record_date,'2018_11_28')
        meta.subject.dateID='D01';
        
        if strcmp(meta.path.run_num,'run12')
            meta.subject.runID='R01';
        else
            error('Invalid run ID given');
        end
    else
        error('Invalid date ID given');
    end
elseif strcmp(meta.path.patID,'TS04 Double DBS Implantation')
    meta.subject.subjID='S04';
    
    if strcmp(meta.path.record_date,'2017_03_01')
        meta.subject.dateID='D01';
        
        if strcmp(meta.path.run_num,'run16')
            meta.subject.runID='R01';
        else
            error('Invalid run ID given');
        end
    else
        error('Invalid date ID given');
    end
elseif strcmp(meta.path.patID,'ET_CL_001')
    meta.subject.subjID='S05';
    
    if strcmp(meta.path.record_date,'2017_05_17')
        meta.subject.dateID='D01';
        
        if strcmp(meta.path.run_num,'run12')
            meta.subject.runID='R01';
        elseif strcmp(meta.path.run_num,'run13')
            meta.subject.runID='R02';
        elseif strcmp(meta.path.run_num,'run15')
            meta.subject.runID='R03';
        elseif strcmp(meta.path.run_num,'run17')
            meta.subject.runID='R04';
        else
            error('Invalid run ID given');
        end
    else
        error('Invalid date ID given');
    end
elseif strcmp(meta.path.patID,'ET_CL_006')
    meta.subject.subjID='S06';
    
    if strcmp(meta.path.record_date,'2020_05_21')
        meta.subject.dateID='D01';
        
        if strcmp(meta.path.run_num,'run3')
            meta.subject.runID='R01';
        else
            error('Invalid run ID given');
        end
    else
        error('Invalid date ID given');
    end
else
    error('Invalid patient ID given');
end

counter=1;

[folder,name,~]=fileparts(original_file);

newFile=fullfile(folder,sprintf('dtf_%s.mat',name));

while exist(newFile,'file') == 2
    fullfile(folder,sprintf('dtf_%s_(%s).mat',name,counter))
    newFile=fullfile(get_root_path,'Files',sprintf('%s_%s_%s%s_(%d).mat',subjID,dateID,runID,addOn,counter));
    counter=counter+1;
end

% newFile=fullfile(get_root_path,'Files',sprintf('%s_%s_%s%s.mat',subjID,dateID,runID,meta.path.addon));
% 
% while exist(newFile,'file') == 2
%     newFile=fullfile(get_root_path,'Files',sprintf('%s_%s_%s%s_(%d).mat',subjID,dateID,runID,meta.path.addOn,counter));
%     counter=counter+1;
% end

end