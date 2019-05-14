%% Miscellaneous code pieces for testing/debugging purposes

%#ok<*UNRCH>

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

return;

%% Code for plotting the periodogram of the signal

pxx=pwelch(X,fs,fs/2,1:(fs/2),fs); 
figure; plot(10*log10(pxx))

return

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

t=(0:length(X)-1)/2400;

figure;
subplot(211);
plot(t,X(:,1),'b'); hold on;
plot(t(modelOrder+1:end),x_hat(:,1),'g');
plot(t(modelOrder+1:end),E(:,1),'r');
legend('X','x-hat','E');

subplot(212);
plot(t,X(:,2),'b'); hold on;
plot(t(modelOrder+1:end),x_hat(:,2),'g');
plot(t(modelOrder+1:end),E(:,2),'r');
legend('X','x-hat','E');



