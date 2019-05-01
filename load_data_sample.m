function [x,fs]=load_data_sample()
%% [x,fs]=load_data_sample()
%
%  Sample load_data function that only loads one file from the server, and filters it for
%  demonstrative purposes.
%
%   Inputs: none
%
%   Outputs:
%    - x: Matrix of two signals from the vim and the cortex. Size is [12151 x 2].
%    - fs: Sampling frequency in Hz. fs = 2400
%
%  See also: varm, estimate, mvar, sample_script
%

FILE='\\gunduz-lab.bme.ufl.edu\Study_ET_Closed_Loop\ET_CL_004\2018_06_20\preproc\run5.mat';
% In visual_stim, stimulus 1 is Rest, 2 is Right hand, 3 is Left hand
data=load(FILE);

fs=data.datastorage.src.LFP.Fs;

vim=data.datastorage.src.LFP.data(:,7)-data.datastorage.src.LFP.data(:,6);
cort=data.datastorage.src.LFP.data(:,4)-data.datastorage.src.LFP.data(:,3);

% Filtering, high pass filter and notch filter

order_hp=4;
cutoff_hp=1;
[num_hp,den_hp]=CreateHPF_butter(fs,order_hp,cutoff_hp);

vim=filtfilt(num_hp,den_hp,vim);
cort=filtfilt(num_hp,den_hp,cort);

f0=60;
w0=f0/(fs/2);
qFactor=35;
bw=w0/qFactor;

[num_comb,den_comb]=iircomb(fs/f0,bw,'notch');

vim=filtfilt(num_comb,den_comb,vim);
cort=filtfilt(num_comb,den_comb,cort);

% Extract one trial to look at

instruct=data.datastorage.src.visual_stim.data;
ind_rest=find(instruct==1);
ind_right=find(instruct==2);

ind_trial=ind_rest(1):ind_right(1);

oneTrial_vim=vim(ind_trial);
oneTrial_cort=cort(ind_trial);

x=[oneTrial_vim,oneTrial_cort];

end