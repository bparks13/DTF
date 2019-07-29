CCC;

% load(fullfile(get_root_path,'Files','S02_D01_R01__MATLABLL.mat'));

%% 2 variable AR model

freqRange=1:100;

N=1000;
m=2;
stdZ=[1,0.7];
a=[0.9,0;0.16,0.8]; a=cat(3,a,[-0.5,0;-0.2,-0.5]);

X=create_data(N,m,stdZ,a);

config=struct(...
    'orderRange',1:30,...
    'crit','aic',...
    'orderSelection','min',...
    'freqRange',freqRange,...
    'logLikelihoodMethod',2 ...
    );

[mdl_aic,E_aic,crit_aic]=mvar(X,config);

config.crit='bic';

[mdl_bic,E_bic,crit_bic]=mvar(X,config);

config.crit='spectra';
config.fs=200;

[mdl_psd,E_psd,crit_psd]=mvar(X,config);

