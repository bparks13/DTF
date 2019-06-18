CCC;

%% Create a multivariate model

colors=linspecer(6);

numSamples=10000;
numSeries=2;
N=[numSamples,numSeries];
modelOrder=[12,30];
stdZ=2;
% a=0.4;
% rho=0.2;
% [X,noise,a,rho]=create_data_correlated_noise(N,modelOrder,a,rho,stdZ); 
[X,noise,a,rho]=create_data_correlated_noise(N,modelOrder,[],[],stdZ); 
% [X,noise,a,rho]=create_data_correlated_noise(N,modelOrder,[],0,stdZ); 

fs=1;

% FILE='\\gunduz-lab.bme.ufl.edu\Study_ET_Closed_Loop\ET_CL_004\2018_06_20\preproc\run5.mat';
% 
% % fs=extract_sampling_frequency(FILE);
% % 
% % order_notch=4;
% % cutoff_notch=[54,66;114,126;176,184;236,244];
% 
% filtering=struct('NO_FILTERING',true);
% 
% % for i=1:length(cutoff_notch)
% %     [filtering.notch(i).num,filtering.notch(i).den]=CreateBSF_butter(fs,order_notch,cutoff_notch(i,:));
% % end
% % 
% % order_hp=4;
% % cutoff_hp=6;
% % [num_hp,den_hp]=CreateHPF_butter(fs,order_hp,cutoff_hp);
% % 
% % filtering.hpf.num=num_hp;
% % filtering.hpf.den=den_hp;
% % 
% % filtering.ma=20;
% 
% X=load_data(FILE,[6,5;4,3],1,filtering);
% X=X(:,:,1);
% 
% numSamples=length(X);
% numSeries=size(X,2);

t=(0:(length(X)-1))/fs;

%%
% 
% modelOrder=6;
% a=nan(numSeries,numSeries,modelOrder);

%% Calculate model coefficients using Yule Walker

tic_yule=tic;
[AR]=estimate_ar_coefficients(X,modelOrder(1));
[E,C,x_hat]=estimate_residuals(X,AR);
logL=calculate_loglikelihood(E,C);
bic=calculate_bic(logL,modelOrder(1),numSamples-modelOrder(1));

% config=struct('changeDOF',true,'numParameters',modelOrder);
% [pass,h,pVal]=test_model(E,length(X),config);
[pass,h,pVal]=test_model(E,length(X));

%% Calculate model coefficients for the correlated noise using Yule Walker

AR_noise=estimate_ar_coefficients(E,modelOrder(2));
[E_noise,C_noise,x_hat_noise]=estimate_residuals(E,AR_noise);
logL_noise=calculate_loglikelihood(E_noise,C_noise);
bic_noise=calculate_bic(logL_noise,modelOrder(2),numSamples-modelOrder(2));

[pass_noise,h_noise,pVal_noise]=test_model(E_noise,length(E));

%% Capture timing

fprintf('Yule-Walker: '); toc(tic_yule);

%% Calculate model coefficients using varm

tic_varm=tic;
mdl=varm(numSeries,modelOrder(1));
[estMdl,~,logL_varm,E_varm]=estimate(mdl,X);
results=summarize(estMdl);
bic_varm=results.BIC;

% config=struct('changeDOF',true,'numParameters',modelOrder);
% [pass_varm,h_varm,pVal_varm]=test_model(E_varm,length(X),config);
[pass_varm,h_varm,pVal_varm]=test_model(E_varm,length(X));

%% Calculate model coefficients for the correlated noise using varm

mdl_noise=varm(numSeries,modelOrder(2));
[estMdl_noise,~,logL_varm_noise,E_varm_noise]=estimate(mdl_noise,E_varm);
results_noise=summarize(estMdl_noise);
bic_varm_noise=results_noise.BIC;

% config=struct('changeDOF',true,'numParameters',modelOrder);
% [pass_varm,h_varm,pVal_varm]=test_model(E_varm,length(X),config);
[pass_varm_noise,h_varm_noise,pVal_varm_noise]=test_model(E_varm_noise,length(E_varm));

%% Capture timing

fprintf('varm: '); toc(tic_varm);

%% Plot original signal, estimated signal, original noise, and estimated noise for Yule Walker

figure('Name','Yule-Walker'); 

for i=1:numSeries
    subplot(numSeries,1,i);
    plot(t,X(:,i),'Color',colors(1,:)); hold on;
    plot(t(modelOrder(1)+1:end),x_hat(:,i),'Color',colors(2,:));
    plot(t,noise(:,i),'Color',colors(3,:));
    plot(t(modelOrder(1)+1:end),E(:,i),'Color',colors(4,:));
    plot(t(modelOrder(1)+modelOrder(2)+1:end),x_hat_noise(:,i),'Color',colors(5,:));
    plot(t(modelOrder(1)+modelOrder(2)+1:end),E_noise(:,i),'Color',colors(6,:));

    legend('Given X','X-hat','Given Noise','E-hat','X-hat-noise','E-noise'); drawnow;
end

%% Plot original signal, estimated signal, original noise, and estimated noise for varm

figure('Name','Varm'); 

for i=1:numSeries
    subplot(numSeries,1,i);
    plot(t,X(:,i),'Color',colors(1,:)); hold on;
    % plot(t(modelOrder+1:end),x_hat,'Color',colors(2,:));
    plot(t,noise(:,i),'Color',colors(3,:));
    plot(t(modelOrder(1)+1:end),E_varm(:,i),'Color',colors(4,:));

    legend('X','E','E-hat'); drawnow;
end

%% Print out the differences between the models

AR_varm=zeros(size(AR));
AR_varm_noise=zeros(size(AR_noise));

for i=1:modelOrder(1)
    AR_varm(:,:,i)=estMdl.AR{i};
end

for i=1:modelOrder(2)
    AR_varm_noise(:,:,i)=estMdl_noise.AR{i};
end

fprintf('\nCoefficients comparison for m = %d, d = %d. Model Coefficients:\n\ta\t\t\tYule\t\t(a-Y)\t\tvarm\t\t(a-v)\t\t(Y-v)\n',modelOrder(1),numSeries);
for i=1:(size(AR,1) * size(AR,2) * size(AR,3))
    fprintf('\t%-5.3f\t\t%-5.3f\t\t%-5.3f\t\t%-5.3f\t\t%-5.3f\t\t%-.6f\n',...
        a(i),AR(i),a(i)-AR(i),AR_varm(i),a(i)-AR_varm(i),AR(i)-AR_varm(i));
end

fprintf('\nCoefficients comparison for m = %d, d = %d. Correlated Noise Coefficients:\n\trho\t\t\tYule\t\t(rho-Y)\t\tvarm\t\t(rho-v)\t\t(Y-v)\n',modelOrder(2),numSeries);
for i=1:(size(AR_noise,1) * size(AR_noise,2) * size(AR_noise,3))
    fprintf('\t%-5.3f\t\t%-5.3f\t\t%-5.3f\t\t%-5.3f\t\t%-5.3f\t\t%-.6f\n',...
        rho(i),AR_noise(i),rho(i)-AR_noise(i),AR_varm_noise(i),rho(i)-AR_varm_noise(i),...
        AR_noise(i)-AR_varm_noise(i));
end

fprintf('\nResiduals comparison:\n\tRoot-Mean-Square of model residual differences = %.5f\n',...
    rms(E-E_varm));

fprintf('\nResiduals comparison:\n\tRoot-Mean-Square of correlated noise residual differences = %.5f\n',...
    rms(E_noise-E_varm_noise));

fprintf('\nCovariance matrix comparison for model:\n\tYule:\t\t\tVarm:\t\t\tDiff:\n');

for i=1:numSeries^2
    fprintf('\t%-6.2f\t\t\t%-6.2f\t\t\t%.6f\n',C(i),estMdl.Covariance(i),C(i)-estMdl.Covariance(i));
end

fprintf('\nCovariance matrix comparison for correlated noise:\n\tYule:\t\t\tVarm:\t\t\tDiff:\n');

for i=1:numSeries^2
    fprintf('\t%-6.2f\t\t\t%-6.2f\t\t\t%.6f\n',C_noise(i),estMdl_noise.Covariance(i),...
        C_noise(i)-estMdl_noise.Covariance(i));
end

fprintf('\nLog-Likelihood of model:\n\tlogL_y\t\t\tlogL_v\t\t\t(logL_y-logL_v)\n\t%-10.2f\t\t%-10.2f\t\t%.5f\n',...
    logL,logL_varm,logL-logL_varm);

fprintf('\nLog-Likelihood of correlated noise:\n\tlogL_y\t\t\tlogL_v\t\t\t(logL_y-logL_v)\n\t%-10.2f\t\t%-10.2f\t\t%.5f\n',...
    logL_noise,logL_varm_noise,logL_noise-logL_varm_noise);

fprintf('\nBIC comparison of model:\n\tBIC_y\t\t\tBIC_v\t\t\t(BIC_y-BIC_v)\n\t%-10.2f\t\t%-10.2f\t\t%.5f\n',...
    bic,bic_varm,bic-bic_varm);

fprintf('\nBIC comparison of correlated noise:\n\tBIC_y\t\t\tBIC_v\t\t\t(BIC_y-BIC_v)\n\t%-10.2f\t\t%-10.2f\t\t%.5f\n',...
    bic_noise,bic_varm_noise,bic_noise-bic_varm_noise);

fprintf('\nWhiteness comparison of model:\n\t\t\t\tYule\tvarm\n');
fprintf('\tPass:\t\t%d\t\t%d\n',pass,pass_varm);

for i=1:numSeries
    fprintf('\th:\t\t\t%d\t\t%d\n\tP-value:\t%.3f\t%.3f\n',h(i),h_varm(i),pVal(i),...
        pVal_varm(i));
end

fprintf('\nWhiteness comparison of correlated noise:\n\t\t\t\tYule\tvarm\n');
fprintf('\tPass:\t\t%d\t\t%d\n',pass_noise,pass_varm_noise);

for i=1:numSeries
    fprintf('\th:\t\t\t%d\t\t%d\n\tP-value:\t%.3f\t%.3f\n',h_noise(i),h_varm_noise(i),...
        pVal_noise(i),pVal_varm_noise(i));
end








