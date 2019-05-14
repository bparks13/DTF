CCC;

%% Create a multivariate model

numSamples=10000;
numSeries=4;
N=[numSamples,numSeries];
modelOrder=2;
stdZ=1;
[X,a]=create_data(N,modelOrder,stdZ); 

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
% X=load_data(FILE,[2,1],1,filtering);
% X=X(:,1,1);
% 
% modelOrder=30;
% N=length(X);
% a=nan(modelOrder,1);

%% Calculate model coefficients and estimated variances using Yule Walker

tic_yule=tic;
[AR,C]=estimate_ar_coefficients(X,modelOrder);
[E,x_hat]=estimate_residuals(X,AR);
logL=calculate_loglikelihood(E,C,modelOrder);
bic=calculate_bic(logL,modelOrder,numSamples-modelOrder);

% config=struct('changeDOF',true,'numParameters',modelOrder);
% [pass,h,pVal]=test_model(E,length(X),config);
[pass,h,pVal]=test_model(E,length(X));
fprintf('Yule-Walker: '); toc(tic_yule);

%% Calculate model coefficients using varm

tic_varm=tic;
mdl=varm(numSeries,modelOrder);
[estMdl,~,logL_varm,E_varm]=estimate(mdl,X);
results=summarize(estMdl);
bic_mdl=results.BIC;

% config=struct('changeDOF',true,'numParameters',modelOrder);
% [pass_varm,h_varm,pVal_varm]=test_model(E_varm,length(X),config);
[pass_varm,h_varm,pVal_varm]=test_model(E_varm,length(X));
fprintf('varm: '); toc(tic_varm);

%% Print out the differences between the models

AR_varm=zeros(size(a));

for i=1:modelOrder
    AR_varm(:,:,i)=estMdl.AR{i};
end

fprintf('\nCoefficients comparison for m = %d, d = %d:\n\ta\t\t\tYule\t\t(a-Y)\t\tvarm\t\t(a-v)\t\t(Y-v)\n',modelOrder,numSeries);
for i=1:(size(a,1) * size(a,2) * size(a,3))
    fprintf('\t%.3f\t\t%.3f\t\t%.3f\t\t%.3f\t\t%.3f\t\t%.6f\n',...
        a(i),AR(i),a(i)-AR(i),AR_varm(i),a(i)-AR_varm(i),AR(i)-AR_varm(i));
end

fprintf('\nResiduals comparison:\n\tAverage difference in residuals = %.5f\n',mean(E-E_varm));

fprintf('\nLog-Likelihood comparison:\n\tlogL_y\t\t\tlogL_v\t\t(logL_y-logL_v)\n\t%.2f\t\t%.2f\t\t%.5f\n',logL,logL_varm,logL-logL_varm);

fprintf('\nBIC comparison:\n\tBIC_y\t\t\tBIC_v\t\t(BIC_y-BIC_v)\n\t%.2f\t\t%.2f\t\t%.5f\n',bic,bic_mdl,bic-bic_mdl);

fprintf('\nWhiteness comparison:\n\t\t\t\tYule\tvarm\n');
fprintf('\tPass:\t\t%d\t\t%d\n',pass,pass_varm);

for i=1:modelOrder
    fprintf('\th:\t\t\t%d\t\t%d\n\tP-value:\t%.3f\t%.3f\n',h(i),h_varm(i),pVal(i),pVal_varm(i));
end

