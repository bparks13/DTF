%% 

CCC;

%% 2 variable AR(2) model

fs=200;
freqRange=4:100;
spectral_range=freqRange(1):0.1:freqRange(end);
normalizeSpectra=true;

N=1000;
m=2;
stdZ=sqrt([1,0.7]); C_orig=[stdZ(1),0;0,stdZ(2)].^2;
a=[0.9,0;0.16,0.8]; a=cat(3,a,[-0.5,0;-0.2,-0.5]);

X=create_data(N,m,stdZ,a);

config=struct(...
    'orderRange',1:20,...
    'crit','aic',...
    'fs',fs,...
    'method','arfit',...
    'orderSelection','min',...
    'freqRange',freqRange,...
    'logLikelihoodMethod',2 ...
    );

%%
[mdl_aic,E_aic,crit_aic]=mvar(X,config);

config.crit='bic';

[mdl_bic,E_bic,crit_bic]=mvar(X,config);

config.crit='spectra';
config.spectral_range=spectral_range;
config.normalizeSpectra=true;
% config.simulated=struct('a',a,'C',C_orig);

[mdl_psd,E_psd,crit_psd]=mvar(X,config);

%% 

AR=estimate_ar_coefficients(X,4,'arfit');
[tmp_E,C,x_hat]=estimate_residuals(X,AR);
% [S_orig,tmp_sr]=calculate_fft(X,fs,normalizeSpectra);
tmp_sr=1:0.5:(fs/2);
S_orig=pwelch(X,fs,fs/2,tmp_sr,fs);
S=calculate_ar_spectra(AR,tmp_sr,fs,C,normalizeSpectra);
% S_orig=calculate_ar_spectra(a,spectral_range,fs,C_orig);

plot_spectra(S,tmp_sr,S_orig);

%% Playing around with using pwelch instead of fft

tmp_sr=1:0.5:(fs/2);
tmp_crit=nan(15,1);

for i=1:15
    AR=estimate_ar_coefficients(X,i,'arfit');
    [tmp_E,C,tmp_xhat]=estimate_residuals(X,AR);
    [tmp_X]=create_data(length(tmp_E),i,diag(C),AR);
    S_orig=pwelch(X,fs,fs/2,tmp_sr,fs);
    S=pwelch(tmp_X,fs,fs/2,tmp_sr,fs);
%     S=pwelch(tmp_xhat,fs,fs/2,tmp_sr,fs);
%     S=calculate_ar_spectra(AR,tmp_sr,fs,C,normalizeSpectra);

    tmp_crit(i)=mean(mean((S-sqrt(S_orig)).^2));
    plot_spectra(S,tmp_sr,S_orig);
    waitforbuttonpress
end

%% 

CCC;

%% 5 variable AR(3) model

a=zeros(5,5,3);
% a=zeros(3,3,3);
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
C_orig=zeros(5);
% stdZ=sqrt([0.6 0.5 0.3]);
% C_orig=zeros(3);
C_orig=C_orig+diag(stdZ).^2;

fs=200;
freqRange=4:100;
spectral_range=freqRange(1):0.1:freqRange(end);

N=1000;
m=3;

X=create_data(N,m,stdZ,a);

config=struct(...
    'orderRange',1:20,...
    'crit','aic',...
    'method','arfit',...
    'fs',fs,...
    'orderSelection','min',...
    'freqRange',freqRange,...
    'logLikelihoodMethod',2 ...
    );

%%

[mdl_aic,E_aic,crit_aic]=mvar(X,config);

config.crit='bic';

[mdl_bic,E_bic,crit_bic]=mvar(X,config);

config.crit='spectra';
config.spectral_range=spectral_range;
% config.simulated=struct('a',a,'C',C_orig);

[mdl_psd,E_psd,crit_psd]=mvar(X,config);

%%

AR=estimate_ar_coefficients(X,4,'arfit');
[tmp_E,C,x_hat]=estimate_residuals(X,AR);
[P1,tmp_sr]=calculate_fft(X,fs);
S=calculate_ar_spectra(AR,tmp_sr,fs,C,true);
% S=calculate_ar_spectra(AR,spectral_range,fs,C);
% S_orig=calculate_ar_spectra(a,spectral_range,fs,C_orig);

plot_spectra(S,tmp_sr,P1);

