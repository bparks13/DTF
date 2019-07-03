function [x_filt,filt_values]=filter_serial_correlation(file,config)
%% [x_filt,filt_values]=filter_serial_correlation(file,config)
%
%  If the AR model has correlated errors, filter the original signal to get rid of the
%  correlated errors. Iterates through error model orders sequentially until the signal
%  passes the Portmanteau Ljung Box Test for serial correlation. Returns the filtered
%  signal only; need to rerun mvar to reacqure new AR coefficients
%
%   Inputs:
%    - file: String containing the name of the file only, not the path to the file. Can be
%       empty to open up the uigetfile interface
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

if nargin == 0 || isempty(file)
    file=uigetfile(fullfile(get_root_path,'Files'));
    
    if file == 0
        x_filt=[];
        return
    end
end

maxIterations=10;
errorModelOrders=1:15;

if nargin == 2 && isstruct(config)
    if isfield(config,'maxIterations')
        maxIterations=config.maxIterations;
    end
    
    if isfield(config,'modelOrders')
        errorModelOrders=config.modelOrders;
    end
end

data=load(fullfile(get_root_path,'Files',file));
x=data.x;
res=data.res;
h=data.h;

config_e=data.config_crit;
config_e.output=0;

config_filt=data.config_crit;
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

%% Either return the structs, or save the file

if nargout == 0
    save(fullfile(get_root_path,'Files',file),'-append',x_filt,filt_values);
    clear x_filt filt_values
end

end