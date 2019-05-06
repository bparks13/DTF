function fs=extract_sampling_frequency(file)
%% fs=extract_sampling_frequency(file)
%
%  Simple function designed to open the file given, and extract the sampling frequency to
%  be used for creating filters
%
%   Inputs:
%    - file: Complete file path to the file
%
%   Outputs:
%    - fs: Sampling frequency in Hz
%

tmp=load(file,'datastorage');
fs=tmp.datastorage.src.LFP.Fs;

end