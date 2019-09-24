function [x_filt,filt_values]=filter_serial_correlation(x,res,h,config_crit,config)
%% [x_filt,filt_values]=filter_serial_correlation(x,res,h,config_crit,config)
%
%  If the AR model has correlated errors, filter the original signal to get rid of the
%  correlated errors. Iterates through error model orders sequentially until the signal
%  passes the Portmanteau Ljung Box Test for serial correlation. Returns the filtered
%  signal only; need to rerun mvar to reacqure new AR coefficients
%
%   Inputs:
%    - x: Struct containing the original time series data. The first fields are the
%       conditions from the runs, and each field is a matrix of size is [n x c], where n
%       is the number of samples  and c is the number of channels
%    - res: Struct containing the residuals of the model fit, used for testing the
%       whiteness of the model. Fields are the same as x. Size is [(n - o) x c], where n
%       is the number of samples, o is the model order, and c is the number of channels
%    - h: Struct containing the hypothesis rejection values of each series. Fields are the
%       same as x. '0' means that the null hypothesis is not rejected, while '1' means
%       that the null hypothesis is rejected . Size is [t x c], where t is the number of
%       trials in that condition and c is the number of channels
%    - config_crit: Struct containing additional parameters, as defined in mvar.m
%    - config: Optional struct containing additional parameters
%       maxIterations: Maximum number of times the signal will be filtered at a particular
%         error model order before the next model order is attempted. Default is 10
%       modelOrders: Vector containing the model orders to be attempted for the error AR
%         model. Default is 1:15
%
%   Outputs:
%    - x_filt: Filtered signal as a struct, in the same format as the struct 'x'.
%       Contains the same subfields and number of samples/trials/channels in each field
%    - filt_values: Struct that matches the format of x_filt, and contains the error
%       model order that was filtered, as well as the iteration it was filtered on
%
%  See also: filter_signal, test_model, mvar
%

maxIterations=10;
errorModelOrders=1:15;

if nargin == 5 && isstruct(config)
    if isfield(config,'maxIterations')
        maxIterations=config.maxIterations;
    end
    
    if isfield(config,'modelOrders')
        errorModelOrders=config.modelOrders;
    end
end

config_e=config_crit;
config_e.output=0;

config_filt=config_crit;
config_filt.orderRange=1:15;
config_filt.output=0;

fields=fieldnames(x);
c=cell(length(fields),1);

x_filt=cell2struct(c,fields,1);
filt_values=cell2struct(c,fields,1);

numChannels=size(x.(fields{1}),2);

for i=1:length(fields)
    numTrials=size(x.(fields{i}),3);
    filt_values.(fields{i})=struct('decorrelated',nan,'order',nan,'iteration',nan);
    
    for j=1:numTrials
        sig_orig=x.(fields{i})(:,:,j);
        E_orig=res.(fields{i})(j).E;
        
        tmp_h=h.(fields{i})(j,:);
        tmp_pass=all(tmp_h==0);
        
        bool_decorrelated=false;
        
        if tmp_pass
            x_filt.(fields{i})(:,:,j)=sig_orig;
            filt_values.(fields{i})(j).decorrelated=true;
            fprintf('\n%s: Trial %d already uncorrelated',fields{i},j);
        else
            for k=1:length(errorModelOrders)
                currModelOrder=errorModelOrders(k);
                config_e.orderRange=currModelOrder;
                sig=sig_orig;
                tmp_E=E_orig;
                
                fprintf('\n%s: Trial %d; Error model order %d - ',fields{i},j,k);
                
                for l=1:maxIterations
                    fprintf('%d, ',l);
                    
                    [e_mdl,~,~]=mvar(tmp_E,config_e);
                    order=e_mdl.order;
                    phi_hat=zeros(numChannels,numChannels,order);
                    
                    for m=1:order
                        phi_hat(:,:,m)=diag(diag(e_mdl.AR(:,:,m)));
                    end
                    
                    sig_filt=filter_signal(sig,phi_hat);
                    
                    [~,tmp_E,~]=mvar(sig_filt,config_filt);
                    [tmp_pass,~,~]=test_model(tmp_E,length(tmp_E));
                    
                    sig=sig_filt;
                    
                    if tmp_pass
                        x_filt.(fields{i})(:,:,j)=sig;
                        filt_values.(fields{i})(j).decorrelated=true;
                        filt_values.(fields{i})(j).order=currModelOrder;
                        filt_values.(fields{i})(j).iteration=l;
                        fprintf('\n%s: Trial %d decorrelated',fields{i},j);
                        bool_decorrelated=true;
                        break
                    end
                end
                
                if bool_decorrelated
                    break
                end
            end
                
            if ~bool_decorrelated
                fprintf('\nWARNING: %s Trial %d was not decorrelated\n',fields{i},j);
                x_filt.(fields{i})(:,:,j)=nan(size(sig));
                filt_values.(fields{i})(j).decorrelated=false;
            end
        end
    end
end

% Print a new line to flush the command line output
fprintf('\n');

end