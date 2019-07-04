CCC;

%% Load file

load(fullfile(get_root_path,'Files',uigetfile(fullfile(get_root_path,'Files\*.mat'))));

%% Split a single trial into multiple realizations

lengthInSecs=1;
lengthInSamples=lengthInSecs*fs;
trial=x.Rest(:,:,2);
numSamples=length(trial);
numChannels=size(trial,2);
numRealizations=floor(numSamples/(lengthInSamples));

realizations=zeros(lengthInSamples,numChannels,numRealizations);

for i=1:numRealizations
    realizations(:,:,i)=trial(lengthInSamples*(i-1)+1:i*lengthInSamples,:);
end

%% Model it

tmp_ar=struct('mdl',[]);
tmp_res=struct('E',[]);
tmp_crit=struct('bic',[]);
c=config_crit;

for i=1:numRealizations
    [tmp_ar(i).mdl,tmp_res(i).E,tmp_crit(i).(config_crit.crit)]=mvar(realizations(:,:,i),c);
end

avgCrit=zeros(length(c.orderRange),1);

for i=1:numRealizations
    avgCrit=avgCrit+tmp_crit(i).(config_crit.crit);
end

avgCrit=avgCrit/numRealizations;

%% Plot the realizations

t=(0:lengthInSamples-1)/fs;

hFig=figure;

for i=1:numRealizations
    figure(hFig);
    modelOrder=tmp_ar(i).mdl.order;
    
    for j=1:numChannels
        subplot(numChannels,1,j); cla; 
        plot(t,realizations(:,j,i),'b'); hold on;
        plot(t(modelOrder+1:end),tmp_ar(i).mdl.x_hat(:,j),'g');
        plot(t(modelOrder+1:end),tmp_res(i).E(:,j),'r');
%         plot(t(1:end-modelOrder),tmp_ar(i).mdl.x_hat(:,j),'g');
%         plot(t(1:end-modelOrder),tmp_res(i).E(:,j),'r');
        legend('X','x\_hat','E');
    end
    
    waitforbuttonpress
end
