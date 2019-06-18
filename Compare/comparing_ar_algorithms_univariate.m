CCC;

%% Create a univariate model

% N=10000;
% modelOrder=30;
% stdZ=2;
% [X,a]=create_data(N,modelOrder,stdZ); % parameters are automatically: 100000, 2, 2, [0.5 0.2]'

FILE='\\gunduz-lab.bme.ufl.edu\Study_ET_Closed_Loop\ET_CL_004\2018_06_20\preproc\run5.mat';

% fs=extract_sampling_frequency(FILE);
% 
% order_notch=4;
% cutoff_notch=[54,66;114,126;176,184;236,244];

filtering=struct('NO_FILTERING',true);

% for i=1:length(cutoff_notch)
%     [filtering.notch(i).num,filtering.notch(i).den]=CreateBSF_butter(fs,order_notch,cutoff_notch(i,:));
% end
% 
% order_hp=4;
% cutoff_hp=6;
% [num_hp,den_hp]=CreateHPF_butter(fs,order_hp,cutoff_hp);
% 
% filtering.hpf.num=num_hp;
% filtering.hpf.den=den_hp;
% 
% filtering.ma=20;

X=load_data(FILE,[2,1],1,filtering);
X=X(:,1,1);

%%

modelOrder=20;
N=length(X);
a=nan(modelOrder,1);

%% Calculate model coefficients and estimated variances using Yule Walker

tic_yule=tic;
[AR,C]=estimate_ar_coefficients(X,modelOrder);
[E,C_2]=estimate_residuals(X,AR);
logL=calculate_loglikelihood(E,C);
bic=calculate_bic(logL,modelOrder,N-modelOrder);

% config=struct('changeDOF',true,'numParameters',modelOrder);
% [pass,h,pVal]=test_model(E,length(X),config);
[pass,h,pVal]=test_model(E,length(X));
fprintf('Yule-Walker: '); toc(tic_yule);

%% Calculate model coefficients using varm

tic_varm=tic;
mdl=varm(1,modelOrder);
[estMdl,~,logL_varm,E_varm]=estimate(mdl,X);
results=summarize(estMdl);
bic_mdl=results.BIC;

% config=struct('changeDOF',true,'numParameters',modelOrder);
% [pass_varm,h_varm,pVal_varm]=test_model(E_varm,length(X),config);
[pass_varm,h_varm,pVal_varm]=test_model(E_varm,length(X));
fprintf('varm: '); toc(tic_varm);

%% Print out the differences between the models

fprintf('\nCoefficients comparison for m = %d:\n\ta\t\t\tYule\t\t(a-Y)\t\tvarm\t\t(a-v)\t\t(Y-v)\n',modelOrder);
for i=1:length(a)
    fprintf('\t%-5.3f\t\t%-5.3f\t\t%-5.3f\t\t%-5.3f\t\t%-5.3f\t\t%.6f\n',...
        a(i),AR(i),a(i)-AR(i),estMdl.AR{i},a(i)-estMdl.AR{i},AR(i)-estMdl.AR{i});
end

fprintf('\nResiduals comparison:\n\tAverage difference in residuals = %.5f\n',mean(E-E_varm));

fprintf('\nCovariance matrix comparison:\n\tYule 1:\t\t\tYule 2:\t\t\tVarm:\t\t\tDiff 1:\t\t\tDiff 2\n');

fprintf('\t%6.3f\t\t\t%6.3f\t\t\t%6.3f\t\t\t%.6f\t\t%.6f\n',C,C_2,estMdl.Covariance,C-estMdl.Covariance,C_2-estMdl.Covariance);

fprintf('\nLog-Likelihood comparison:\n\tlogL_y\t\t\tlogL_v\t\t(logL_y-logL_v)\n\t%.2f\t\t%.2f\t\t%.5f\n',logL,logL_varm,logL-logL_varm);

fprintf('\nBIC comparison:\n\tBIC_y\t\t\tBIC_v\t\t(BIC_y-BIC_v)\n\t%.2f\t\t%.2f\t\t%.5f\n',bic,bic_mdl,bic-bic_mdl);

fprintf('\nWhiteness comparison:\n\t\t\t\tYule\tvarm\n');
fprintf('\tPass:\t\t%d\t\t%d\n\th:\t\t\t%d\t\t%d\n\tP-value:\t%.3f\t%.3f\n',pass,pass_varm,h,h_varm,pVal,pVal_varm);

