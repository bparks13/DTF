function [mdl,E,criterion]=mvar(x,config)
%% [mdl,E,criterion]=mvar(x,config)
%  
%  Given a matrix of signals, and optionally a configuration struct, construct a
%  multi-variate autoregressive model of the signals.
%
%   Inputs:
%    - x: Signals in a matrix format [n x s], where n is the number of samples and s is
%       the number of series to model
%    - config: Struct containing optional parameters
%       orderRange: Vector containing the model orders to consider. Default [1:30]
%       crit: String defining which information criterion to use, 'aic', 'bic' [default],
%           'psd', or 'spectra'
%       output: Int defining level of verbosity for output. 0 (none), 1 (model number)
%           [default], 2 (model order and criterion tested)
%       method: String defining which method to use; Matlab's varm ('varm') or Yule-Walker
%           equations ('yule') [default]
%       orderSelection: String defining how to algorithmically choose the model order.
%           'min' [default] uses the minimum information criterion found, while 'diff'
%           uses the first model order that the abs(difference) between successive ICs is
%           smaller than some epsilon, which is either 0.01 or user specified
%       epsilon: Can be specified if orderSelection is set to 'diff', where epsilon is the
%           threshold for defining when differences have decreased to a small enough
%           degree to select the model order, defined as a percentage (i.e. 0.01 == 1%).
%           Default value is 0.1% (= 0.001)
%       fs: If crit is defined as 'psd', the sampling frequency is required
%       freqRange: If crit is defined as 'psd', the PSD will be calculated up to half the
%           sampling frequency. freqRange can be defined to restrict the analysis of the
%           PSD to a smaller range
%       logLikelihoodMethod: Int defining which equation to use for calculating the log
%           likelihood of the model. 1 is for Matlab [default], 2 for Ding, 3 for
%           Awareness paper (see calculate_loglikelihood)
%
%   Outputs:
%    - mdl: Struct containing the AR model fit for the data given
%       AR: Autoregressive coefficients as found by estimate_ar_coefficients
%       C: Covariance matrx as calculated by the Yule-Walker equations
%       logL: Log-likelihood of the model fit, used for calculating information criterion
%       order: Model order that is found to have the lowest information criterion
%       numSeries: Number of series in the model
%       x_hat: Estimated data from the AR model
%       pxx: Estimated Power Spectral Density of the x_hat, if crit is defined as 'psd'
%    - E: Residuals of the model fit, used for testing the whiteness of the model. Size is
%       [(n - o) x s], where n is the number of samples, o is the model order, and s is
%       the number of series
%    - criterion: Optional output, contains the information criterion chosen in config
%
% See also: estimate_ar_coefficients, estimate_residuals, calculate_loglikelihood,
%   calculate_bic
%

% Defaults
orderRange=1:30;
crit='bic';
output=1;
method='yule';
orderSelection='min';
epsilon=0.001;
pxx_sig=[];
pxx_ar=[];
ll_method=1; % Log-Likelihood method; 1 == Matlab, 2 == Ding

if nargin > 1 && isstruct(config)
    if isfield(config,'orderRange')
        orderRange=config.orderRange;
    end
    
    if isfield(config,'crit')
        crit=config.crit;
        
        if strcmp(crit,'psd')
            fs=config.fs;
            window=round(fs);
            overlap=round(fs/2);
            freqRange=1:round(fs/2);
            pxx_sig=pwelch(x,window,overlap,freqRange,fs);
            
            if isfield(config,'freqRange')
                freqForAnalysis=config.freqRange;
            else
                freqForAnalysis=freqRange;
            end
        elseif strcmp(crit,'spectra')
            if isfield(config,'simulated')
                N=length(x);
                freqRange=1:round(N/2);
                S_orig=calculate_ar_spectra(config.simulated.a,freqRange,2*freqRange(end),config.simulated.C);
                
                P1=nan(size(S_orig,3),size(S_orig,1)); 
            
                for j=1:size(S_orig,1)
                    P1(:,j)=squeeze(abs(S_orig(j,j,:)));
                end
            end
% %             fs=config.fs;
% %             freqRange=1:round(fs/2);
%             N=length(x);
%             freqRange=1:round(N/2);
%             
% %             if isfield(config,'freqRange')
% %                 freqForAnalysis=config.freqRange;
% %             else
% %                 freqForAnalysis=freqRange;
% %             end
%             
%             Y=fft(x);
%             P2=abs(Y/N);
%             P1=P2(2:floor(N/2)+1,:); % Ignore the DC component
%             P1(1:end-1,:)=2*P1(1:end-1,:);
% %             P1=P1(freqForAnalysis,:);
        end 
    end
    
    if isfield(config,'output')
        output=config.output;
    end
    
    if isfield(config,'method')
        method=config.method;
    end
    
    if isfield(config,'orderSelection')
        orderSelection=config.orderSelection;
        
        if isfield(config,'epsilon')
            epsilon=config.epsilon;
        end
    end
    
    if isfield(config,'logLikelihoodMethod')
        ll_method=config.logLikelihoodMethod;
    end
end

numSeries=size(x,2);
numOrders=length(orderRange);
numSamples=length(x);

criterion=zeros(numOrders,1);

minCrit=inf;
bool_minDiffFound=false;

mdl=struct('AR',[],'C',[],'logL',[],'order',[],'numSeries',numSeries);

if output ~= 0
    fprintf('Beginning order estimation:\n');
end

for i=1:numOrders
    if strcmp(method,'varm')
        mdl_varm=varm(numSeries,i);
        [estMdl,~,logL,tmp_E]=estimate(mdl_varm,x);
        results=summarize(estMdl);
        
        if strcmp(crit,'bic')
            criterion(i)=results.BIC;
            
            if strcmp(orderSelection,'min')
                result = results.BIC < minCrit;
            elseif strcmp(orderSelection,'diff')
                if i>1 && ~bool_minDiffFound
                    result = abs((criterion(i) - criterion(i-1)) / criterion(i)) < epsilon || ...
                        criterion(i) - criterion(i-1) > 0;
                    bool_minDiffFound=true;
                else
                    result=false;
                end
            end
            
            if result
                minCrit=results.BIC;
                mdl.AR=zeros(numSeries,numSeries,orderRange(i));
                for j=1:orderRange(i)
                    mdl.AR(:,:,j)=estMdl.AR{j};
                end
                mdl.C=estMdl.Covariance;
                mdl.logL=logL;
                mdl.order=orderRange(i);
                E=tmp_E;
                
                if strcmp(orderSelection,'diff')
                    bool_minDiffFound=true;
                end
            end
        elseif strcmp(crit,'aic')
            criterion(i)=results.AIC;
            
            if results.AIC < minCrit
                minCrit=results.AIC;
            end
        end
    elseif strcmp(method,'yule')
        [AR]=estimate_ar_coefficients(x,orderRange(i));
        [tmp_E,C,x_hat]=estimate_residuals(x,AR);
        logL=calculate_loglikelihood(tmp_E,C,ll_method);

        if strcmp(crit,'bic')
            criterion(i)=calculate_bic(logL,orderRange(i),numSeries,numSamples,ll_method);
            
            if strcmp(orderSelection,'min')
                result = criterion(i) < minCrit;
            elseif strcmp(orderSelection,'diff')
                if i>1 && ~bool_minDiffFound
                    result = abs((criterion(i) - criterion(i-1)) / criterion(i)) < epsilon || ...
                        criterion(i) - criterion(i-1) > 0;
                else
                    result=false;
                end
            end
        elseif strcmp(crit,'aic')
            criterion(i)=calculate_aic(logL,orderRange(i),numSeries,numSamples,ll_method);
            
            if strcmp(orderSelection,'min')
                result = criterion(i) < minCrit;
            elseif strcmp(orderSelection,'diff')
                if i>1 && ~bool_minDiffFound
                    result = abs((criterion(i) - criterion(i-1)) / criterion(i)) < epsilon || ...
                        criterion(i) - criterion(i-1) > 0;
                else
                    result=false;
                end
            end
        elseif strcmp(crit,'psd')
            pxx_ar=pwelch(x_hat,window,overlap,freqRange,fs);
            
            if numSeries == 1
                criterion(i)=mean(mean(abs(pxx_sig(freqForAnalysis) - pxx_ar(freqForAnalysis))));
            else
                criterion(i)=mean(mean(abs(pxx_sig(freqForAnalysis,:) - pxx_ar(freqForAnalysis,:))));
            end
            
            if strcmp(orderSelection,'min')
                result = criterion(i) < minCrit;
            elseif strcmp(orderSelection,'diff')
                if i>1 && ~bool_minDiffFound
                    result = abs((criterion(i) - criterion(i-1)) / criterion(i)) < epsilon || ...
                        criterion(i) - criterion(i-1) > 0;
                else
                    result=false;
                end
            end
        elseif strcmp(crit,'spectra')
%             S_ar=calculate_ar_spectra(AR,freqForAnalysis,fs,C);
            S_ar=calculate_ar_spectra(AR,freqRange,2*freqRange(end),C);
            
            pxx_ar=nan(size(S_ar,3),size(S_ar,1)); % Not actually pxx, but maintaining consistency of the naming scheme 
            
            for j=1:size(S_ar,1)
                pxx_ar(:,j)=squeeze(abs(S_ar(j,j,:)));
            end
            
            criterion(i)=mean(mean((P1-pxx_ar).^2));
            
            if strcmp(orderSelection,'min')
                result = criterion(i) < minCrit;
            elseif strcmp(orderSelection,'diff')
                if i>1 && ~bool_minDiffFound
                    result = abs((criterion(i) - criterion(i-1)) / criterion(i)) < epsilon || ...
                        criterion(i) - criterion(i-1) > 0;
                else
                    result=false;
                end
            end
            
        end
        
        if result
            minCrit=criterion(i);
            mdl.AR=AR;
            mdl.C=C;
            mdl.logL=logL;
            mdl.order=orderRange(i);
            mdl.x_hat=x_hat;
            mdl.pxx=pxx_ar;
            E=tmp_E;

            if strcmp(orderSelection,'diff')
                bool_minDiffFound=true;
            end
        end
    end
    
    if output == 1
        fprintf('%d,',i);
    elseif output == 2
        fprintf('%d: Current %s = %.2f. Minimum %s = %.2f\n',i,crit,criterion(i),crit,minCrit);
    end
end

if isempty(mdl.order)
    fprintf('\nWARNING: NO MODEL FOUND WITH ''%s'' AND ''%s''\nFINDING MINIMUM VALUE...\n',orderSelection,crit);
    
    minCritInd=find(criterion==min(criterion));
    
    [mdl.AR]=estimate_ar_coefficients(x,minCritInd);
    [E,mdl.C]=estimate_residuals(x,mdl.AR);
    mdl.logL=calculate_loglikelihood(E,mdl.C,ll_method);
    mdl.order=minCritInd;
end

if output ~= 0
    fprintf('\nDone: Minimum %s found at model order %d\n',crit,mdl.order);
end

end

