CCC;

%% 

file='\\gunduz-lab.bme.ufl.edu\Study_ET_Closed_Loop\ET_CL_004\2018_06_20\preproc\run5';

fs=extract_sampling_frequency(file);

PATIENT_ID='ET_CL_004';
RECORDING_DATE='2018_06_20';
RUN_ID='run5';

filtering=struct;
filtering.NO_FILTERING=true;
[filtering.hpf.num,filtering.hpf.den]=CreateHPF_butter(fs,3,4);
% filtering.normalize='z-score';
% order_notch=4;
% cutoff_notch=[54,66;114,126;176,184;236,244];
% 
% for i=1:length(cutoff_notch)
%     [filtering.notch(i).num,filtering.notch(i).den]=CreateBSF_butter(fs,order_notch,cutoff_notch(i,:));
% end
% 
% [filtering.lpf.num,filtering.lpf.den]=CreateLPF_butter(fs,8,300);

freqRange=1:1200;
[channels,labels,conditions,cond_labels]=load_channels_labels_conditions(PATIENT_ID,RECORDING_DATE,RUN_ID);

[x,~,x_all]=load_data(file,channels,conditions(1),filtering);

numChannels=length(channels);

orderRange=1:30;

window=round(fs);
overlap=round(window/2);

currTrial=6;

pxx_sig=pwelch(x(:,:,currTrial),window,overlap,freqRange,fs);

%%

% figure;
% subplot(221); plot(freqRange,10*log10(pxx_sig(:,1)),'k'); hold on; title(labels{1}); 
% subplot(222); plot(freqRange,10*log10(pxx_sig(:,2)),'k'); hold on; title(labels{2}); 
% subplot(223); plot(freqRange,10*log10(pxx_sig(:,3)),'k'); hold on; title(labels{3}); 
% subplot(224); plot(freqRange,10*log10(pxx_sig(:,4)),'k'); hold on; title(labels{4}); 

diff_pxx=nan(length(freqRange),numChannels,length(orderRange));

for i=orderRange
    ar_coeff=estimate_ar_coefficients(x(:,:,currTrial),orderRange(i));
    [E,C,x_hat]=estimate_residuals(x(:,:,currTrial),ar_coeff);
    
    pxx_ar=pwelch(x_hat,window,overlap,freqRange,fs);
    
    diff_pxx(:,:,i)=abs(10*log(pxx_sig)-10*log10(pxx_ar)); disp(i);
    
%     subplot(221); plot(freqRange,10*log10(pxx_ar(:,1)),'r'); 
%     subplot(222); plot(freqRange,10*log10(pxx_ar(:,2)),'r');  
%     subplot(223); plot(freqRange,10*log10(pxx_ar(:,3)),'r'); 
%     subplot(224); plot(freqRange,10*log10(pxx_ar(:,4)),'r'); 
%     
%     drawnow
%     waitforbuttonpress;
end

return;

%%
%#ok<*UNRCH>
% [E,C,x_hat]=estimate_residuals(x(:,:,1),ar_coeff);

t=(0:length(x(:,1,1))-1)/fs;

figure;

subplot(411); plot(t,x(:,1,currTrial)); hold on; 
plot(t(orderRange(i)+1:end),x_hat(:,1)); plot(t(orderRange(i)+1:end),E(:,1));
subplot(412); plot(t,x(:,2,currTrial)); hold on; 
plot(t(orderRange(i)+1:end),x_hat(:,2)); plot(t(orderRange(i)+1:end),E(:,2));
subplot(413); plot(t,x(:,3,currTrial)); hold on; 
plot(t(orderRange(i)+1:end),x_hat(:,3)); plot(t(orderRange(i)+1:end),E(:,3));
subplot(414); plot(t,x(:,4,currTrial)); hold on; 
plot(t(orderRange(i)+1:end),x_hat(:,4)); plot(t(orderRange(i)+1:end),E(:,4));


%% Plot overall differences in PSD

figure

plot(squeeze(mean(mean(diff_pxx,1))))
find(mean(mean(diff_pxx))==min(mean(mean(diff_pxx))))

%% Plot individual differences in PSD

figure;

subplot(411); plot(squeeze(mean(diff_pxx(:,1,:),1)));
subplot(412); plot(squeeze(mean(diff_pxx(:,2,:),1)));
subplot(413); plot(squeeze(mean(diff_pxx(:,3,:),1)));
subplot(414); plot(squeeze(mean(diff_pxx(:,4,:),1)));





