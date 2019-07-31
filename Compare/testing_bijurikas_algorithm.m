%% 

CCC;

%% 2 variable AR(2) model
freqRange=4:100;

N=1000;
m=2;
stdZ=sqrt([1,0.7]); C_orig=[stdZ(1),0;0,stdZ(2)].^2;
a=[0.9,0;0.16,0.8]; a=cat(3,a,[-0.5,0;-0.2,-0.5]);

X=create_data(N,m,stdZ,a);

config=struct(...
    'orderRange',1:20,...
    'crit','aic',...
    'orderSelection','min',...
    'freqRange',freqRange,...
    'logLikelihoodMethod',2 ...
    );

[mdl_aic,E_aic,crit_aic]=mvar(X,config);

config.crit='bic';

[mdl_bic,E_bic,crit_bic]=mvar(X,config);

config.crit='spectra';
config.simulated=struct('a',a,'C',C_orig);

[mdl_psd,E_psd,crit_psd]=mvar(X,config);

%% Testing the spectra/FFT calculations and normalizations

% freqRange=1:round(N/2);
% 
% Y=fft(X);
% P2=abs(Y/N);
% P1=P2(2:floor(N/2)+1,:);
% P1(1:end-1,:)=2*P1(1:end-1,:);
% f=freqRange(end)*(1:(N/2))/N;
% 
% AR=estimate_ar_coefficients(X,2);
% [tmp_E,C,x_hat]=estimate_residuals(X,AR);
% S=calculate_ar_spectra(AR,freqRange,2*freqRange(end),C);
% S_orig=calculate_ar_spectra(a,freqRange,2*freqRange(end),C_orig);

%% 

CCC;

%% 5 variable AR(3) model

a=zeros(5,5,3);
a(1,1,1)=0.95*sqrt(2);
% a(3,3,1)=0.25;
a(4,4,1)=0.25*sqrt(2);
a(4,5,1)=0.25*sqrt(2);
a(5,5,1)=0.25*sqrt(2);
a(5,4,1)=-0.25*sqrt(2);
a(1,1,2)=-0.9025;
a(2,1,2)=0.5;
a(4,1,2)=-0.5;
a(3,1,3)=-0.4;

stdZ=sqrt([0.6 0.5 0.3 0.3 0.6]);
C_orig=zeros(5);
C_orig=C_orig+diag(stdZ).^2;

freqRange=4:100;

N=1000;
m=3;

X=create_data(N,m,stdZ,a);

config=struct(...
    'orderRange',1:20,...
    'crit','aic',...
    'orderSelection','min',...
    'freqRange',freqRange,...
    'logLikelihoodMethod',2 ...
    );

[mdl_aic,E_aic,crit_aic]=mvar(X,config);

config.crit='bic';

[mdl_bic,E_bic,crit_bic]=mvar(X,config);

config.crit='spectra';
config.simulated=struct('a',a,'C',C_orig);

[mdl_psd,E_psd,crit_psd]=mvar(X,config);

%%

AR=estimate_ar_coefficients(X,2);
[tmp_E,C,x_hat]=estimate_residuals(X,AR);
S=calculate_ar_spectra(AR,freqRange,2*freqRange(end),C);
S_orig=calculate_ar_spectra(a,freqRange,2*freqRange(end),C_orig);

%% Plot the PSD values
figure
plot(abs(squeeze(S(1,1,:))))
hold on
plot(abs(squeeze(S_orig(1,1,:))))
figure
plot(abs(squeeze(S(2,2,:))))
hold on
plot(abs(squeeze(S_orig(2,2,:))))
figure
plot(abs(squeeze(S(3,3,:))))
hold on
plot(abs(squeeze(S_orig(3,3,:))))
figure
plot(abs(squeeze(S(4,4,:))))
hold on
plot(abs(squeeze(S_orig(4,4,:))))
figure
plot(abs(squeeze(S(5,5,:))))
hold on
plot(abs(squeeze(S_orig(5,5,:))))

%% Comparing Matlab to my code

mdl=varm(5,3);
[estMdl,~,logL_varm,E_varm]=estimate(mdl,X);
results=summarize(estMdl);

[AR,C1]=estimate_ar_coefficients(X,3);
[tmp_E,C2,x_hat]=estimate_residuals(X,AR);


