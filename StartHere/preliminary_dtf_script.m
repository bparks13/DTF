CCC;

%% Load files, and separate variables (ET_CL_04, WORKS)

% FILE='\\gunduz-lab.bme.ufl.edu\Study_ET_Closed_Loop\ET_CL_004\2019_01_24\preproc\run3.mat';
% 
% data=load(FILE);
% 
% fs=data.datastorage.src.LFP.Fs;
% 
% vim=data.datastorage.src.LFP.data(:,1) * 1e3;
% cort=data.datastorage.src.LFP.data(:,2) * 1e3;
% 
% acc_rh=sqrt(data.datastorage.src.ACC.data(:,1).^2 + data.datastorage.src.ACC.data(:,2).^2 + ...
%     data.datastorage.src.ACC.data(:,3).^2);
% acc_rh=acc_rh-mean(acc_rh);
% 
% acc_lh=sqrt(data.datastorage.src.ACC.data(:,4).^2 + data.datastorage.src.ACC.data(:,5).^2 + ...
%     data.datastorage.src.ACC.data(:,6).^2);
% acc_lh=acc_lh-mean(acc_lh);
% 
% dac=data.datastorage.src.DAC.data;
% 
% % Cue Right 1
% cue_rh=(dac > 1.2 & dac < 1.4);
% start=find(diff(cue_rh)==1);
% stop=find(diff(cue_rh)==-1);
% 
% for i=1:length(start)
%     if stop(i)-start(i) < 100
%         cue_rh(start(i):stop(i))=0;
%     end
% end
% 
% % Move Right 1
% move_rh=(dac > 1.4 & dac < 1.6);
% start=find(diff(move_rh)==1);
% stop=find(diff(move_rh)==-1);
% 
% for i=1:length(start)
%     if stop(i)-start(i) < 100
%         move_rh(start(i):stop(i))=0;
%     end
% end
% 
% % Cue Left 1
% cue_lh=(dac > 0.8 & dac < 1.0);
% start=find(diff(cue_lh)==1);
% stop=find(diff(cue_lh)==-1);
% 
% for i=1:length(start)
%     if stop(i)-start(i) < 100
%         cue_lh(start(i):stop(i))=0;
%     end
% end
% 
% % Move Left 1
% move_lh=(dac > 1.0 & dac < 1.2);
% start=find(diff(move_lh)==1);
% stop=find(diff(move_lh)==-1);
% 
% for i=1:length(start)
%     if stop(i)-start(i) < 100
%         move_lh(start(i):stop(i))=0;
%     end
% end
% 
% ind_trial_start=find(move_rh==1,1);
% ind_trial_end=find(move_rh(ind_trial_start:end)==0,1)+ind_trial_start;
% 
% oneTrial_vim=vim(ind_trial_start:ind_trial_end);
% oneTrial_cort=cort(ind_trial_start:ind_trial_end);

%% Load files, and separate variables (TOO MUCH NOISE)

% % % % % % FILE='\\gunduz-lab.bme.ufl.edu\Study_ET\OR\withDBS\ET_OR_STIM_018\2018_11_28\preproc\run14.mat'; 
% % % % % % % in visual_stim, stimulus 1 is Rest, 2 is Right hand, 3 is Left hand, 4 is Caress
% % % % % % % yellow contact is 1, blue is motor cortex
% % % % % % 
% % % % % % data=load(FILE);
% % % % % % 
% % % % % % fs=data.datastorage.src.LFP.Fs;
% % % % % % 
% % % % % % vim_3=data.datastorage.src.LFP.data(1:end-40000,10); 
% % % % % % vim_2=data.datastorage.src.LFP.data(1:end-40000,9); 
% % % % % % vim_1=data.datastorage.src.LFP.data(1:end-40000,8); 
% % % % % % vim_12=vim_2-vim_1;
% % % % % % cort=data.datastorage.src.LFP.data(1:end-40000,4);
% % % % % % 
% % % % % % % Filtering
% % % % % % 
% % % % % % % order_hp=4;
% % % % % % % cutoff_hp=0.5;
% % % % % % % [num_hp,den_hp]=CreateHPF_butter(fs,order_hp,cutoff_hp);
% % % % % % % 
% % % % % % % vim_3=filtfilt(num_hp,den_hp,vim_3);
% % % % % % % vim_2=filtfilt(num_hp,den_hp,vim_2);
% % % % % % % vim_1=filtfilt(num_hp,den_hp,vim_1);
% % % % % % % vim_12=filtfilt(num_hp,den_hp,vim_12);
% % % % % % % cort=filtfilt(num_hp,den_hp,cort);
% % % % % % 
% % % % % % % [num_bs,den_bs]=CreateBSF_butter(fs,2,[55 65]);
% % % % % % % 
% % % % % % % vim_3=filtfilt(num_bs,den_bs,vim_3);
% % % % % % % vim_2=filtfilt(num_bs,den_bs,vim_2);
% % % % % % % vim_1=filtfilt(num_bs,den_bs,vim_1);
% % % % % % % vim_12=filtfilt(num_bs,den_bs,vim_12);
% % % % % % % cort=filtfilt(num_bs,den_bs,cort);
% % % % % % % 
% % % % % % % order_lp=6;
% % % % % % % cutoff_lp=30;
% % % % % % % [num_lp,den_lp]=CreateLPF_butter(fs,order_lp,cutoff_lp);
% % % % % % % 
% % % % % % % vim_3=filtfilt(num_lp,den_lp,vim_3);
% % % % % % % vim_2=filtfilt(num_lp,den_lp,vim_2);
% % % % % % % vim_1=filtfilt(num_lp,den_lp,vim_1);
% % % % % % % vim_12=filtfilt(num_lp,den_lp,vim_12);
% % % % % % % cort=filtfilt(num_lp,den_lp,cort);
% % % % % % 
% % % % % % instruct=data.datastorage.src.visual_stim.data;
% % % % % % 
% % % % % % ind_rest=find(instruct==1);
% % % % % % ind_right=find(instruct==2);
% % % % % % ind_left=find(instruct==3);
% % % % % % ind_caress=find(instruct==4);
% % % % % % 
% % % % % % ind_trial=ind_rest(1):ind_right(1);
% % % % % % 
% % % % % % oneTrial_vim=vim_2(ind_trial);
% % % % % % oneTrial_cort=cort(ind_trial);

%% Load files, and separate variables (ET_CL_04,INTRAOP)

FILE='\\gunduz-lab.bme.ufl.edu\Study_ET_Closed_Loop\ET_CL_004\2018_06_20\preproc\run5.mat';
% In visual_stim, stimulus 1 is Rest, 2 is Right hand, 3 is Left hand
data=load(FILE);

fs=data.datastorage.src.LFP.Fs;

vim=data.datastorage.src.LFP.data(:,7)-data.datastorage.src.LFP.data(:,6);
cort=data.datastorage.src.LFP.data(:,4)-data.datastorage.src.LFP.data(:,3);

% Filtering

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

% [num_bs,den_bs]=CreateBSF_butter(fs,2,[58 62]);
% 
% vim=filtfilt(num_bs,den_bs,vim);
% cort=filtfilt(num_bs,den_bs,cort);
% 
% [num_bs,den_bs]=CreateBSF_butter(fs,2,[118 122]);
% 
% vim=filtfilt(num_bs,den_bs,vim);
% cort=filtfilt(num_bs,den_bs,cort);

% [num_bs,den_bs]=CreateBSF_butter(fs,2,[130 135]);
% 
% vim=filtfilt(num_bs,den_bs,vim);
% cort=filtfilt(num_bs,den_bs,cort);

% [num_bs,den_bs]=CreateBSF_butter(fs,2,[178 182]);
% 
% vim=filtfilt(num_bs,den_bs,vim);
% cort=filtfilt(num_bs,den_bs,cort);

% Extract one trial

instruct=data.datastorage.src.visual_stim.data;
ind_rest=find(instruct==1);
ind_right=find(instruct==2);

ind_trial=ind_rest(1):ind_right(1);

oneTrial_vim=vim(ind_trial);
oneTrial_cort=cort(ind_trial);

x=[oneTrial_vim,oneTrial_cort];

%% Calculate the MVAR model

numSeries=2;
orderRange=1:30;

BIC=zeros(max(orderRange),1);
minBIC=inf;

fprintf('Beginning order estimation:\n');

for i=orderRange
    mdl=varm(numSeries,i);

    [estMdl,~,~,~]=estimate(mdl,[oneTrial_vim,oneTrial_cort]);

    results=summarize(estMdl);
    
    BIC(i)=results.BIC;
    
    if results.BIC < minBIC
        minBIC=results.BIC;
    end
    
    fprintf('%d,',i);
end

modelOrder=find(BIC==minBIC);
% modelOrder=15;

fprintf('\nDone: Minimum BIC at Model Order = %d\n',modelOrder);

mdl=varm(numSeries,modelOrder);

[estMdl,estSE,logL,E]=estimate(mdl,[oneTrial_vim,oneTrial_cort]);

%% Test for whiteness

lags=round(log(length(oneTrial_cort)));
[h_1,pValue_1,stat_1,cValue_1]=lbqtest(E(:,1),'lags',lags); 
[h_2,pValue_2,stat_2,cValue_2]=lbqtest(E(:,2),'lags',lags); 
% Lags changed based on Matlab documentation at https://www.mathworks.com/help/econ/ljung-box-q-test.html
% Consider using a higher number of lags, and also changing the degrees of freedom as
% described in the help function for the test based on the number of parameters of the
% model

% for i=1:20
%     [h_1,pValue_1,stat_1,cValue_1]=lbqtest(E(:,1),'lags',i); 
%     [h_2,pValue_2,stat_2,cValue_2]=lbqtest(E(:,2),'lags',i); 
%     fprintf('%d - 1: h = %d, c = %.2f. 2: h = %d, c = %.2f\n',i,h_1,cValue_1,h_2,cValue_2);
% end

%% Calculate DTF

freqRange=2:100;
nFreqs=length(freqRange);
I=eye(numSeries);

AR=zeros([size(estMdl.AR{1}),length(estMdl.AR)]);  % Extract the coefficients from the model, [numseries x numseries x modelorder]
for i=1:length(estMdl.AR)
    AR(:,:,i)=estMdl.AR{i};
end

tmp_pdc=zeros(numSeries,numSeries,nFreqs);

% Calculate PDC from individual calculations
for i=1:nFreqs
%     for k=1:modelOrder
%         tmp_pdc(:,:,i)=tmp_pdc(:,:,i)+AR(:,:,k).*exp(-(2*pi*1i/fs)*k*freqRange(i));
%     end
    tmp_pdc(:,:,i)=I-reshape(sum(bsxfun(@times,reshape(AR,numSeries^2,modelOrder),...
        exp(-(2*pi*1i/fs)*(1:modelOrder)*freqRange(i))),2),numSeries,numSeries);
% This is the matrix calculation from SIFT
end

% Calculate H from the inverse of the PDC above
H=zeros(size(tmp_pdc));

for i=1:nFreqs
    H(:,:,i)=inv(tmp_pdc(:,:,i));
end

% Calculate theta^2

theta=zeros(size(H));

for i=1:size(H,1)
    for j=1:size(H,2)
        for k=1:size(H,3)
            theta(i,j,k)=abs(H(i,j,k))^2;
        end
    end
end

%% Calculate gamma^2

gamma=zeros(size(H));

for k=1:size(H,3)
    for i=1:size(H,1)
        den=sum(abs(H(i,:,k)).^2);
        for j=1:size(H,2)
            gamma(i,j,k)=abs(H(i,j,k))^2/den;
        end
    end
end

%% Plotting

figure;

ax=zeros(4,1);

ax(1)=subplot(221); plot(squeeze(gamma(1,1,:))); 
ax(2)=subplot(222); plot(squeeze(gamma(1,2,:))); title('Vim ? cortex');
ax(3)=subplot(223); plot(squeeze(gamma(2,1,:))); title('cortex ? Vim');
ax(4)=subplot(224); plot(squeeze(gamma(2,2,:)));

linkaxes(ax)














