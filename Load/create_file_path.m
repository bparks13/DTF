function [file,file_postacq]=create_file_path(meta)
%% [file,file_postacq]=create_file_path(meta)
% 
%  Given the meta variable, create the filenames for the original data and the postacq data
%
%   Inputs:
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
%   Outputs:
%    - file: String defining the absolute filepath to the original data file
%    - file_postacq: String defining the absolute filepath to the postacquisition data
%       file
%

file=fullfile(meta.path.pre,meta.path.patID,meta.path.record_date,meta.path.midpath,meta.path.run_num);

file_postacq=fullfile(meta.path.pre,meta.path.midpath_postacq,meta.path.patID_postacq,meta.path.record_date,['postacq_',meta.path.run_num]);

end