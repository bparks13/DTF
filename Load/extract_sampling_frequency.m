function fs=extract_sampling_frequency(file)
%% fs=extract_sampling_frequency(file)
%
%  Simple function designed to open the file given, and extract the sampling frequency to
%  be used for creating filters. The datastorage struct can also be given as the input to
%  minimize the number of times the file needs to be opened
%
%   Inputs:
%    - file: Complete file path to the file. OR The datastorage struct. If the datastorage
%       struct, it is assumed that it is the datastorage struct itself
%
%   Outputs:
%    - fs: Sampling frequency in Hz
%
%  See also: processing_pipeline
%

if ischar(file)
    tmp=load(file,'datastorage');
    fs=tmp.datastorage.src.LFP.Fs;
elseif isstruct(file)
    fs=file.src.LFP.Fs;
else
    error('Invalid input');
end

end