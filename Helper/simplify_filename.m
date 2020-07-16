function newFile=simplify_filename(meta)
%% newFile=simplify_filename(meta)
%
%  Instead of using the entire subject ID/date of recording/run ID, simplify the new
%  filename by assigning a number to each input; i.e. meta.patID == 'ET_CL_04' becomes S02,
%  meta.record_date == '2018_06_20' becomes D01, etc.
%
%   Inputs:
%    - meta: Struct containing all data relevant to this run
%    -- patID: String containing the patient ID
%    -- record_date: String containing the recording date
%    -- run_num: String containing the run number
%    -- addon: String containing anything to add at the end of the file name
%
%   Outputs:
%    - newFile: String containing the absolute path to the folder that contains all of the
%       processed files; guaranteed to be unique
%
%  See also: processing_pipeline
%

if strcmp(meta.path.patID,'ET_CL_002')
    subjID='S01';
    
    if strcmp(meta.path.record_date,'2018_02_01')
        dateID='D01';
        
        if strcmp(meta.path.run_num,'run9')
            runID='R01';
        elseif strcmp(meta.path.run_num,'run10')
            runID='R02';
        elseif strcmp(meta.path.run_num,'run11')
            runID='R03';
        elseif strcmp(meta.path.run_num,'run12')
            runID='R04';
        else
            error('Invalid run ID given');
        end
    else
        error('Invalid date ID given');
    end
elseif strcmp(meta.path.patID,'ET_CL_004')
    subjID='S02';
    
    if strcmp(meta.path.record_date,'2018_06_20')
        dateID='D01';
        
        if strcmp(meta.path.run_num,'run5')
            runID='R01';
        else
            error('Invalid run ID given');
        end
    elseif strcmp(meta.path.record_date,'2018_08_23')
        dateID='D02';
        
        if strcmp(meta.path.run_num,'run1')
            runID='R01';
        else
            error('Invalid run ID given');
        end
    else
        error('Invalid date ID given');
    end
elseif strcmp(meta.path.patID,'ET_OR_STIM_018')
    subjID='S03';
    
    if strcmp(meta.path.record_date,'2018_11_28')
        dateID='D01';
        
        if strcmp(meta.path.run_num,'run12')
            runID='R01';
        else
            error('Invalid run ID given');
        end
    else
        error('Invalid date ID given');
    end
elseif strcmp(meta.path.patID,'TS04 Double DBS Implantation')
    subjID='S04';
    
    if strcmp(meta.path.record_date,'2017_03_01')
        dateID='D01';
        
        if strcmp(meta.path.run_num,'run16')
            runID='R01';
        else
            error('Invalid run ID given');
        end
    else
        error('Invalid date ID given');
    end
elseif strcmp(meta.path.patID,'ET_CL_001')
    subjID='S05';
    
    if strcmp(meta.path.record_date,'2017_05_17')
        dateID='D01';
        
        if strcmp(meta.path.run_num,'run12')
            runID='R01';
        elseif strcmp(meta.path.run_num,'run13')
            runID='R02';
        elseif strcmp(meta.path.run_num,'run15')
            runID='R03';
        elseif strcmp(meta.path.run_num,'run17')
            runID='R04';
        else
            error('Invalid run ID given');
        end
    else
        error('Invalid date ID given');
    end
elseif strcmp(meta.path.patID,'ET_CL_006')
    subjID='S06';
    
    if strcmp(meta.path.record_date,'2020_05_21')
        dateID='D01';
        
        if strcmp(meta.path.run_num,'run3')
            runID='R01';
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

newFile=fullfile(get_root_path,'Files',sprintf('%s_%s_%s%s.mat',subjID,dateID,runID,meta.path.addon));

while exist(newFile,'file') == 2
    newFile=fullfile(get_root_path,'Files',sprintf('%s_%s_%s%s_(%d).mat',subjID,dateID,runID,addOn,counter));
    counter=counter+1;
end

end