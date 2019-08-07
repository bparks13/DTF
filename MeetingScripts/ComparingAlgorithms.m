%% 
% Comparing between developed algorithms, Bijurika/Ding's published data, and the Signal
% Processing for Neuroscientists book
%

CCC;

%% Run 'ExternalFunctions/pr14_2.m' and save figures

% Make sure relevant variables from that script are still loaded in memory

%% Run the same signal through my code

fs=1000;
freqRange=1:0.2:500;

config=struct(...
'orderRange',1:20,...
'crit','bic',...
'fs',fs,...
'freqRange',freqRange,...
'method','arfit',...
'orderSelection','min',...
'logLikelihoodMethod',2 ...
);

mdl_orig=struct('AR',AR,'order',13,'numSeries',3);
[conn_orig,H_orig]=dtf(mdl_orig,freqRange,fs,config);

config_plot=struct('seriesType',6);
plot_connectivity(conn_orig,H_orig,freqRange,{'Sig1','Sig2','Sig3'},config_plot);

hFig=gcf;
hFig.Position=[hFig.Position(1:2) 850 650];

[mdl,E,crit]=mvar(X',config);

[connectivity,H]=dtf(mdl,freqRange,fs);

config_plot=struct('seriesType',6);
plot_connectivity(connectivity,H,freqRange,{'Sig1','Sig2','Sig3'},config_plot);
hFig=gcf;
hFig.Position=[hFig.Position(1:2) 850 650];


%% Bijurika 2 variable AR(2) model - Figure 1

fs=200;
freqRange=4:100;
spectral_range=freqRange(1):0.1:freqRange(end);
orderRange=1:20;

N=1000;
m=2;
stdZ=sqrt([1,0.7]); C_orig=diag(stdZ.^2);
a=[0.9,0;0.16,0.8]; a=cat(3,a,[-0.5,0;-0.2,-0.5]);

X1=create_data(N,m,stdZ,a);

config=struct(...
    'orderRange',orderRange,...
    'crit','aic',...
    'fs',fs,...
    'method','arfit',...
    'orderSelection','min',...
    'freqRange',freqRange,...
    'logLikelihoodMethod',2 ...
    );

[mdl_aic,E_aic,crit_aic]=mvar(X1,config); %#ok<*ASGLU>

config.crit='bic';

[mdl_bic,E_bic,crit_bic]=mvar(X1,config);

config.crit='spectra';
config.spectral_range=spectral_range;
config.normalizeSpectra=true;
config.orderSelection='diff2';

[mdl_mse1,E_mse1,crit_mse1]=mvar(X1,config);

figure;
subplot(221); plot(orderRange,crit_aic); title('AIC'); xlabel('Model Order'); 
subplot(222); plot(orderRange,crit_bic); title('BIC'); xlabel('Model Order'); 
subplot(223); plot(orderRange,crit_mse1); title('MSE (1,000 Points)'); xlabel('Model Order'); 

N=10000;
X2=create_data(N,m,stdZ,a);
[mdl_mse2,E_mse2,crit_mse2]=mvar(X2,config);

subplot(224); plot(orderRange,crit_mse2); title('MSE (10,000 Points)'); xlabel('Model Order'); 

hFig=gcf;
hFig.Position=[hFig.Position(1:2) 720 500];

%% Bijurika 2 variable AR(2) model - Figure 2

mdl_orig=struct('AR',a,'order',2,'numSeries',2);
[conn_orig,H_orig]=dtf(mdl_orig,freqRange,fs);

S_orig=calculate_spectra(H_orig,C_orig);

config_plot=struct('seriesType',7);
plot_connectivity(conn_orig,S_orig,freqRange,{'Sig1','Sig2'},config_plot);

%% Bijurika 2 variable AR(2) model - Figure 3

% Using MSE 1 for the modeling

[conn_est,H_est]=dtf(mdl_mse1,freqRange,fs);

S_est=calculate_spectra(H_est,mdl_mse1.C);

config_plot=struct('seriesType',7);
plot_connectivity(conn_est,S_est,freqRange,{'Sig1','Sig2'},config_plot);

%% Bijurika 5 variable AR(3) model - Figure 1

orderRange=1:20;

a=zeros(5,5,3);
a(1,1,1)=0.95*sqrt(2);
a(4,4,1)=0.25*sqrt(2);
a(4,5,1)=0.25*sqrt(2);
a(5,5,1)=0.25*sqrt(2);
a(5,4,1)=-0.25*sqrt(2);
a(1,1,2)=-0.9025;
a(2,1,2)=0.5;
a(4,1,2)=-0.5;
a(3,1,3)=-0.4;

stdZ=sqrt([0.6 0.5 0.3 0.3 0.6]);
C_orig=diag(stdZ.^2);

fs=200;
freqRange=4:100;
spectral_range=freqRange(1):0.1:freqRange(end);

N=1000;
m=3;

X=create_data(N,m,stdZ,a);

config=struct(...
    'orderRange',orderRange,...
    'crit','aic',...
    'method','arfit',...
    'fs',fs,...
    'orderSelection','diff1',...
    'freqRange',freqRange,...
    'logLikelihoodMethod',2 ...
    );

[mdl_aic,E_aic,crit_aic]=mvar(X,config);

config.crit='bic';

[mdl_bic,E_bic,crit_bic]=mvar(X,config);

config.crit='spectra';
config.spectral_range=spectral_range;
config.normalizeSpectra=true;
% config.orderSelection='diff1';

[mdl_mse1,E_mse1,crit_mse1]=mvar(X,config);

figure;
subplot(221); plot(orderRange,crit_aic); title('AIC'); xlabel('Model Order'); 
subplot(222); plot(orderRange,crit_bic); title('BIC'); xlabel('Model Order'); 
subplot(223); plot(orderRange,crit_mse1); title('MSE (1,000 Points)'); xlabel('Model Order'); 

N=10000;
X2=create_data(N,m,stdZ,a);
[mdl_mse2,E_mse2,crit_mse2]=mvar(X2,config);

subplot(224); plot(orderRange,crit_mse2); title('MSE (10,000 Points)'); xlabel('Model Order'); 

hFig=gcf;
hFig.Position=[hFig.Position(1:2) 720 500];

%% Bijurika 5 variable AR(3) model - Figure 2

mdl_orig=struct('AR',a,'order',3,'numSeries',5);
[conn_orig,H_orig]=dtf(mdl_orig,freqRange,fs);

S_orig=calculate_spectra(H_orig,C_orig);

config_plot=struct('seriesType',7);
plot_connectivity(conn_orig,S_orig,freqRange,{'Sig1','Sig2','Sig3','Sig4','Sig5'},config_plot);
hFig=gcf;
hFig.Position=[hFig.Position(1:2) 1338 1026];

%% Bijurika 5 variable AR(3) model - Figure 3

% Using bic for modeling

[conn_bic,H_bic]=dtf(mdl_bic,freqRange,fs);

S_bic=calculate_spectra(H_bic,mdl_bic.C);

config_plot=struct('seriesType',7);
plot_connectivity(conn_bic,S_bic,freqRange,{'Sig1','Sig2','Sig3','Sig4','Sig5'},config_plot);
hFig=gcf;
hFig.Position=[hFig.Position(1:2) 1338 1026];

%% Directed Transfer Function: Unified Asymptotic Theory and Some of Its Implications, 2 Var AR(2)

fs=1;
freqRange=0:0.01:0.5;
spectral_range=freqRange(1):0.002:freqRange(end);
orderRange=1:20;

N=100;
m=2;
stdZ=sqrt([1,1]); C_orig=diag(stdZ.^2);
a=[0.95*sqrt(2),0;-0.5,0.5]; a=cat(3,a,[-0.9025,0;0,0]);

X1=create_data(N,m,stdZ,a);

config=struct(...
    'orderRange',orderRange,...
    'crit','aic',...
    'fs',fs,...
    'method','arfit',...
    'orderSelection','min',...
    'freqRange',freqRange,...
    'logLikelihoodMethod',2 ...
    );

[mdl_aic,E_aic,crit_aic]=mvar(X1,config); %#ok<*ASGLU>

config.crit='bic';

[mdl_bic,E_bic,crit_bic]=mvar(X1,config);

config.crit='spectra';
config.spectral_range=spectral_range;
config.normalizeSpectra=true;
config.orderSelection='diff';

[mdl_mse1,E_mse1,crit_mse1]=mvar(X1,config);

figure;
subplot(221); plot(orderRange,crit_aic); title('AIC');
subplot(222); plot(orderRange,crit_bic); title('BIC');
subplot(223); plot(orderRange,crit_mse1); title('MSE (1,000 Points)'); 

N=500;
X2=create_data(N,m,stdZ,a);
[mdl_mse2,E_mse2,crit_mse2]=mvar(X2,config);

subplot(224); plot(orderRange,crit_mse2); title('MSE (10,000 Points)'); 

%% Directed Transfer Function: Unified Asymptotic Theory and Some of Its Implications, 2 variable AR(2) model - Figure 2

mdl_orig=struct('AR',a,'order',2,'numSeries',2);
[conn_orig,H_orig]=dtf(mdl_orig,spectral_range,fs);

S_orig=calculate_spectra(H_orig,C_orig);
S_orig=resize_spectra(S_orig);
S_orig(:,1)=S_orig(:,1)/max(S_orig(:,1));
S_orig(:,2)=S_orig(:,2)/max(S_orig(:,2));

config_plot=struct('seriesType',7);
plot_connectivity(conn_orig,S_orig,spectral_range,{'Sig1','Sig2'},config_plot);
subplot(221); ylim([0 1]);
subplot(223); ylim([0 1]);

%% Directed Transfer Function: Unified Asymptotic Theory and Some of Its Implications, 2 variable AR(2) model - Figure 3

[conn_est,H_est]=dtf(mdl_bic,spectral_range,fs);

S_est=calculate_spectra(H_est,mdl_bic.C);
S_est=resize_spectra(S_est);
S_est(:,1)=S_est(:,1)/max(S_est(:,1));
S_est(:,2)=S_est(:,2)/max(S_est(:,2));

config_plot=struct('seriesType',7);
plot_connectivity(conn_est,S_est,spectral_range,{'Sig1','Sig2'},config_plot);
subplot(221); ylim([0 1]);
subplot(223); ylim([0 1]);


