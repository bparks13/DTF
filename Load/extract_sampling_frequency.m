function fs=extract_sampling_frequency(file)
%% fs=extract_sampling_frequency(file)
%
%  Simple function designed to open the file given, and extract the sampling frequency to
%  be used for creating filters. The datastorage struct can also be given as the input to
%  minimize the number of times the file needs to be opened
%
%   Inputs:
%    - file: Complete file path to the file. OR The datastorage struct. If the datastorage
%       struct, it is assumed that the first field level contains only 'datastorage' as a
%       field.
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
    if isfield(file,'datastorage')  
        fs=file.datastorage.src.LFP.Fs;
    else
        error('Wrong struct given.');
    end
else
    error('Invalid input');
end

end