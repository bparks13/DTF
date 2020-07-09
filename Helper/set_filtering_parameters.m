function filtering=set_filtering_parameters(fs_init,fs_new)
%% filtering=set_filtering_parameters(fs_init,fs_new)
%
%  Helper function to set the filtering parameters and return them in a struct for ease of
%  use
%
%   Inputs:
%    - fs_init: Initial sampling frequency in Hz
%    - fs_new: Sampling frequency to downsample the data to, in Hz
%
%   Outputs:
%    - filtering: Struct containing all filtering parameters
%    -- hpf: Struct containing the high-pass filtering parameters
%    --- num: Numerator coefficients for the high pass filter
%    --- den: Denominator coefficients for the high pass filter
%    --- notes: String defining the filter in words, with dynamically created filtering
%         parameters embedded in the string
%    -- downsample: Sampling frequency to downsample the data to, same as fs_new
%    -- normalize: String defining the method of normalization. Often is 'z-score'
%    -- realizations: Struct containing any relevant information about the realizations
%    --- length: Integer defining the length of each realization in samples at the initial
%         sampling frequency (fs_init)
%    -- notch: Struct containing the notch filtering parameters
%    --- num: Numerator coefficients for the notch filter
%    --- den: Denominator coefficients for the notch filter
%    --- notes: String defining the filter in words, with dynamically created filtering
%         parameters embedded in the string
%    -- lpf: Struct containing the low-pass filtering parameters
%    --- num: Numerator coefficients for the low pass filter
%    --- den: Denominator coefficients for the low pass filter
%    --- notes: String defining the filter in words, with dynamically created filtering
%         parameters embedded in the string

filtering=struct;

order_hp=3;
cutoff_hp=4;
[filtering.hpf.num,filtering.hpf.den]=CreateHPF_butter(fs_init,order_hp,cutoff_hp);
filtering.hpf.note=sprintf('High-Pass Butterworth filter. Order = %d, cutoff = %d Hz',order_hp,cutoff_hp);

filtering.downsample=fs_new;
filtering.normalize='z-score';
realizationLengthInSeconds=1;
filtering.realizations.length=floor(realizationLengthInSeconds*fs_init);

order_notch=3;
cutoff_notch=[58,62];
[filtering.notch.num,filtering.notch.den]=CreateBSF_butter(fs_init,order_notch,cutoff_notch);
filtering.notch.note=sprintf('Notch Butterworth filter from %d-%d Hz, order = %d',cutoff_notch(1),cutoff_notch(2),order_notch);

order_lp=8;
cutoff_lp=floor(filtering.downsample/2);
[filtering.lpf.num,filtering.lpf.den]=CreateLPF_butter(fs_init,order_lp,cutoff_lp);
filtering.lpf.note=sprintf('Low-Pass Butterworth filter. Order = %d, cutoff = %d Hz',order_lp,cutoff_lp);


end