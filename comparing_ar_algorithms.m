CCC;

%% Create a model

X=create_data(); % parameters are automatically: 100000, 2, 2, [0.5 0.2]'

%% Calculate model coefficients and estimated variances using Yule Walker

[Theta_hat,SigmaZ_hat]=estimate_ar_coefficients(X,2);

%% Calculate model coefficients using varm

mdl=varm(1,2);

[estMdl,~,~,~]=estimate(mdl,X);




