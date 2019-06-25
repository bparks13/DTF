%% Miscellaneous code pieces for testing/debugging purposes

FILE='\\gunduz-lab.bme.ufl.edu\Study_ET_Closed_Loop\ET_CL_004\2018_06_20\preproc\run5.mat';
data=load(FILE);

fs=data.datastorage.src.LFP.Fs;

tmp_x=data.datastorage.src.LFP.data(1:end-1,:);

order_hp=4;
cutoff_hp=2;
[num_hp,den_hp]=CreateHPF_butter(fs,order_hp,cutoff_hp);

tmp_x=filtfilt(num_hp,den_hp,tmp_x);

f0=60;
w0=f0/(fs/2);
qFactor=35;
bw=w0/qFactor;
[num_comb,den_comb]=iircomb(fs/f0,bw,'notch');

tmp_x=filtfilt(num_comb,den_comb,tmp_x); %#ok<NASGU>

%% Testing a moving average filter on the data

CCC;

FILE='\\gunduz-lab.bme.ufl.edu\Study_ET_Closed_Loop\ET_CL_004\2018_06_20\preproc\run5.mat';
data=load(FILE);

fs=data.datastorage.src.LFP.Fs;
ma=10;

tmp_x=data.datastorage.src.LFP.data(1:end-1,6)-data.datastorage.src.LFP.data(1:end-1,5);

order_hp=4;
cutoff_hp=5;
[num_hp,den_hp]=CreateHPF_butter(fs,order_hp,cutoff_hp);

tmp_x=filtfilt(num_hp,den_hp,tmp_x);

f0=60;
w0=f0/(fs/2);
qFactor=35;
bw=w0/qFactor;
[num_comb,den_comb]=iircomb(fs/f0,bw,'notch');

tmp_x=filtfilt(num_comb,den_comb,tmp_x);

t=(0:length(tmp_x)-1)/fs;

hFig1=figure;
plot(t,tmp_x); hold on;

hFig2=figure;
pxx=pwelch(tmp_x,fs,fs/2,1:(fs/2),fs);
plot(10*log10(pxx)); hold on;

a=1;
b=ones(ma,1)/ma;

tmp_x=filter(b,a,tmp_x);

figure(hFig1);
plot(t,tmp_x)

pxx=pwelch(tmp_x,fs,fs/2,1:(fs/2),fs);
figure(hFig2);
plot(10*log10(pxx)); hold on;

%% Code for plotting the periodogram of the signal

pxx=pwelch(X,fs,fs/2,1:(fs/2),fs); 
figure; plot(10*log10(pxx));

%% Comparing the PSD of the original signal and the estimated signal

figure;

for i=1:numChannels
    subplot(numChannels,1,i);
    plot(freqRange,pxx_sig(freqRange,i)); hold on;
    plot(freqRange,pxx_ar(freqRange,i));
end

%% Testing using no filtering except for high pass filtering

CCC;

FILE='\\gunduz-lab.bme.ufl.edu\Study_ET_Closed_Loop\ET_CL_004\2018_06_20\preproc\run5.mat';
modelOrder=15;

filtering=struct('NO_FILTERING',true);

X=load_data(FILE,[6,5],1,filtering);
X=X(:,1,1);

N=length(X);

[AR,C]=estimate_ar_coefficients(X,modelOrder);
[E,x_hat]=estimate_residuals(X,AR);
logL=calculate_loglikelihood(E,C);
bic=calculate_bic(logL,modelOrder,N-modelOrder);

[pass,h,pVal]=test_model(E,N);

%% Code for visualizing multivariate data

t=(0:length(x.Rest)-1)/2400;
modelOrder=ar.Rest.mdl.order;

figure;

for i=1:numChannels
    subplot(numChannels,1,i);
    plot(t,tmp_x(:,i),'b'); hold on;
    plot(t(modelOrder+1:end),ar.Rest.mdl.x_hat(:,i),'g');
    plot(t(modelOrder+1:end),res.Rest.E(:,i),'r');
    legend('X','x-hat','E');
end


%% Testing getting the estimated values without actual data

file=fullfile(get_root_path,'Files','ET_CL_004__2018_06_20__run5__200Hz__Z_SCORE__BIC_(1).mat');

data=load(file);
fs=data.fs;
x=data.x.Rest((2.5*fs):(4.5*fs),4,4); % Grabbing WSS data
t=(0:length(x)-1)/fs;
conf=data.config_crit;

[mdl,E,~]=mvar(x,conf);

%% cont'd

numSamples=length(E);
AR=mdl.AR;
m=mdl.order;
noise=randn(numSamples,1);

tmp_xhat=zeros(1,numSamples+m);

for i=1+m:numSamples+m
    tmp_xhat(i)=tmp_xhat(i-1:-1:i-m)*AR+randn(1)*mdl.C;
end

tmp_xhat=tmp_xhat(1+m:end);
    
figure; plot(t,x,'b',t(1+m:end),mdl.x_hat,'g',t(1+m:end),tmp_xhat,'r');
legend('Original','mdl.x\_hat','tmp\_xhat');







