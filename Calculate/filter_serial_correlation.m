function [x_filt,filt_values]=filter_serial_correlation(x,res,h,config_crit,config)
%% [x_filt,filt_values]=filter_serial_correlation(x,res,h,config_crit,config)
%
%  If the AR model has correlated errors, filter the original signal to get rid of the
%  correlated errors. Iterates through error model orders sequentially until the signal
%  passes the Portmanteau Ljung Box Test for serial correlation. Returns the filtered
%  signal only; need to rerun mvar to acquire new AR coefficients
%
%   Inputs:
%    - x: Struct containing the original time series data. The first fields are the
%       conditions from the runs, and each field is a matrix of size is [n x c x r], where
%       n is the number of samples, c is the number of channels, and r is the number of
%       realizations
%    - res: Struct containing the residuals of the model fit, used for testing the
%       whiteness of the model. Fields are the same as x. Size is [(n - o) x c x r], where
%       n is the number of samples, o is the model order, c is the number of channels, and
%       r is the number of realizations
%    - h: Struct containing the hypothesis rejection values of each series. Fields are the
%       same as x. '0' means that the null hypothesis is not rejected, while '1' means
%       that the null hypothesis is rejected . Size is [t x c], where t is the number of
%       trials in that condition and c is the number of channels
%    - config_crit: Struct containing additional parameters, as defined in mvar.m
%    - config: Optional struct containing additional parameters
%    -- maxIterations: Maximum number of times the signal will be filtered at a particular
%        error model order before the next model order is attempted. Default is 10
%    -- modelOrders: Vector containing the model orders to be attempted for the error AR
%        model. Default is 1:15
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

config_filt=config_crit;
config_filt.orderRange=1:15;
config_filt.output=0;

% if isfield(config_filt,'modelOrder')
%     config_filt=rmfield(config_filt,'modelOrder');
% end

fields=fieldnames(x);

numFields=length(fields);

x_filt_cell=cell(numFields,1);
filt_values_cell=cell(numFields,1);
res_cell=struct2cell(res);
h_cell=struct2cell(h);

numChannels=size(x.(fields{1}),2);

x_cell=struct2cell(x);

parfor i=1:numFields
% for i=1:numFields
    numTrials=size(x_cell{i},3);
    filt_values_cell{i}=struct('decorrelated',nan,'order',nan,'iteration',nan);
    x_filt_cell{i}=nan(size(x_cell{i}));
    
    config_e=config_crit;
    config_e.output=0;
    
    for j=1:numTrials
        sig_orig=x_cell{i}(:,:,j);
        E_orig=res_cell{i}(j).E;
        
        tmp_h=h_cell{i}(j,:);
        tmp_pass=all(tmp_h==0);
        
        bool_decorrelated=false;
        
        if tmp_pass
            x_filt_cell{i}(:,:,j)=sig_orig;
            filt_values_cell{i}(j).decorrelated=true;
            fprintf('%s: Trial %d uncorrelated\n',fields{i},j);
        else
            for k=1:length(errorModelOrders)
                currModelOrder=errorModelOrders(k);
                config_e.modelOrder=currModelOrder;
                sig=sig_orig;
                tmp_E=E_orig;
                
                for l=1:maxIterations
                    [e_mdl,~,~]=mvar(tmp_E,config_e);
                    order=e_mdl.order;
                    phi_hat=zeros(numChannels,numChannels,order);
                    
                    for m=1:order
                        phi_hat(:,:,m)=diag(diag(e_mdl.AR(:,:,m)));
                    end
                    
                    sig_filt=filter_signal(sig,phi_hat);
                    
                    [~,tmp_E,~]=mvar(sig_filt,config_filt);
                    [tmp_pass,~,~]=test_model(tmp_E,size(tmp_E,1));
                    
                    sig=sig_filt;
                    
                    if tmp_pass
                        x_filt_cell{i}(:,:,j)=sig;
                        filt_values_cell{i}(j).decorrelated=true;
                        filt_values_cell{i}(j).order=currModelOrder;
                        filt_values_cell{i}(j).iteration=l;
                        fprintf('%s: Trial %d decorrelated with order = %d on iteration = %d\n',fields{i},j,errorModelOrders(k),l);
                        bool_decorrelated=true;
                        break
                    end
                end
                
                if bool_decorrelated
                    break
                end
            end
                
            if ~bool_decorrelated
                fprintf('WARNING: %s Trial %d was not decorrelated\n',fields{i},j);
                x_filt_cell{i}(:,:,j)=nan(size(sig));
                filt_values_cell{i}(j).decorrelated=false;
            end
        end
    end
end

filt_values=cell2struct(filt_values_cell,fields,1);
x_filt=cell2struct(x_filt_cell,fields,1);

end