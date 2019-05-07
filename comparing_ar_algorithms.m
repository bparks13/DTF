CCC;

%% Create a model

X=create_data(); % parameters are automatically: 100000, 2, 2, [0.5 0.2]'

%% Calculate model coefficients and estimated variances using Yule Walker

[AR,C]=estimate_ar_coefficients(X,2);
[E,x_hat]=estimate_residuals(X,AR);
logL=calculate_loglikelihood(E,C);
bic=calculate_bic(logL,2,100000-2);

%% Calculate model coefficients using varm

mdl=varm(1,2);
[estMdl,~,logL_varm,E_varm]=estimate(mdl,X);
results=summarize(estMdl);
bic_mdl=results.BIC;





