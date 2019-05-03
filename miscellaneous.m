% FILE='\\gunduz-lab.bme.ufl.edu\Study_ET_Closed_Loop\ET_CL_004\2018_06_20\preproc\run5.mat';
% data=load(FILE);
% 
% fs=data.datastorage.src.LFP.Fs;
% 
% tmp_x=data.datastorage.src.LFP.data(1:end-1,:);
% 
% order_hp=4;
% cutoff_hp=2;
% [num_hp,den_hp]=CreateHPF_butter(fs,order_hp,cutoff_hp);
% 
% tmp_x=filtfilt(num_hp,den_hp,tmp_x);
% 
% f0=60;
% w0=f0/(fs/2);
% qFactor=35;
% bw=w0/qFactor;
% [num_comb,den_comb]=iircomb(fs/f0,bw,'notch');
% 
% tmp_x=filtfilt(num_comb,den_comb,tmp_x);
% 
