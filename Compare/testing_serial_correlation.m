%% testing_serial_correlation
%
%  This is a script to test the estimation of serially correlated errors
%
%   Based on The Econometric Analysis of Time Series, by Andrew Harvey, as well as the
%   article found at https://www.osapublishing.org/boe/abstract.cfm?uri=boe-4-8-1366
%

CCC;

%% Load

FILE='ET_CL_004__2018_06_20__run5__200Hz__Z_SCORE__BIC_(1).mat';
load(fullfile(get_root_path,'Files',FILE));
currTrial=1;
numChannels=size(x.Rest,2);

%% Testing different model orders for estimating AR coefficients for serially correlated errors

maxOrder=15;
maxIterations=10;
epsilon=0.1;

config_e=config_crit;
config_e.output=0;
% config_e.orderRange=1:maxOrder;

e_crit=zeros(maxOrder,numChannels);
x_sig_orig=x.Rest(:,:,currTrial);
tmp_mdl_orig=ar.Rest(currTrial).mdl;
E_orig=res.Rest(currTrial).E;

t=(0:size(x_sig_orig,1)-1)/fs;

%%

% figure;
% ax=subplot(221);
% plot_acf(E_orig(:,1),[],[],ax); title('Ch. 1')
% ax=subplot(222);
% plot_acf(E_orig(:,2),[],[],ax); title('Ch. 2')
% ax=subplot(223);
% plot_acf(E_orig(:,3),[],[],ax); title('Ch. 3')
% ax=subplot(224);
% plot_acf(E_orig(:,4),[],[],ax); title('Ch. 4')

hFig1=figure('Name','Ch. 1');
hFig2=figure('Name','Ch. 2');
hFig3=figure('Name','Ch. 3');
hFig4=figure('Name','Ch. 4');

%% Iteratively search through all model orders, and choose the lowest BIC for the calculations

x_sig=x_sig_orig;
prevAR=tmp_mdl_orig.AR;
tmp_E=E_orig;

for i=1:maxIterations
    fprintf('Starting error estimation...\n');
    
    bic=nan(maxOrder,1);
    
    for j=1:maxOrder
        phi=calculate_serial_coefficients(tmp_E,j);
        E_filtered=filter_signal(tmp_E,phi);
        logL=calculate_loglikelihood(E_filtered);
        bic(j)=calculate_bic(logL,j,length(tmp_E));
        
%         config_e.orderRange=j;
%         [e_mdl,e_E,e_crit]=mvar(tmp_E,config_e);
%         [e_pass,e_h,e_pval]=test_model(e_E,length(e_E));
        
%         if e_pass
%             break
%         end
    end
    
    phi_min=calculate_serial_coefficients(tmp_E,find(bic==min(bic)));
    
%     x_sig=filter_signal(x_sig,e_mdl.AR);
    x_sig=filter_signal(x_sig,phi_min);
    
    fprintf('Starting signal estimation...\n');
    [tmp_mdl,tmp_E,tmp_crit]=mvar(x_sig,config_crit);
    
    if size(prevAR,3) == size(tmp_mdl.AR,3)
        if all(all(all(prevAR-tmp_mdl.AR < epsilon)))
            fprintf('Convergence.\n')
            break;
        else
            prevAR=tmp_mdl.AR;
        end
    else
        prevAR=tmp_mdl.AR;        
    end
end

plotSignal();

%%
% 
% %% Find an order where the error residuals are uncorrelated; iteratively filter the signal until convergence
% 
% for i=1:maxOrder
%     %%
%     config_crit.orderRange=i;
%     config_crit.output=0;
%     
%     e_mdl=struct;
%     e_E=zeros(length(E_orig)-i,numChannels);
%     
%     for j=1:numChannels
%     	[e_mdl(j).mdl,e_E(:,j),e_crit(i,j)]=mvar(E_orig(:,j),config_crit);
%     end
%     
%     [tmp_pass1,~,~]=test_model(e_E,size(E_orig,1));
%     
%     %%
%     if tmp_pass1
%         x=struct('sig',[]);
%         
%         for j=1:numChannels
%             %%
%             tmp_x_sig=x_sig_orig(:,j);
%             tmp_mdl=e_mdl(j).mdl;
%             
%             for k=1:maxIterations
%                 %%
%                 config_crit.orderRange=1:50;
%                 tmp_x_sig=filter_signal(tmp_x_sig,tmp_mdl.AR);
%                 [~,tmp_E1,~]=mvar(tmp_x_sig,config_crit);
%                 [tmp_pass2,~,~]=test_model(tmp_E1,length(tmp_E1));
%                 
%                 if tmp_pass2
%                     break
%                 end
%                 
%                 for l=1:maxOrder
%                     %%
%                     config_crit.orderRange=l;
%                     [tmp_mdl,tmp_E2,~]=mvar(tmp_E1,config_crit);
% 
%                     [tmp_pass3,~,~]=test_model(tmp_E2,length(tmp_E2));
%                     
%                     if tmp_pass3
%                         break
%                     end
%                 end
%             end
%             
%             x(j).sig=tmp_x_sig;
%         end
%         
%         minLength=inf;
% 
%         for j=1:numChannels
%             if length(x(j).sig) < minLength
%                 minLength=length(x(j).sig);
%             end
%         end
% 
%         x_sig=zeros(minLength,numChannels);
% 
%         for j=1:numChannels
%             x_sig(:,j)=x(j).sig(end-minLength+1:end);
%         end
% 
%         [tmp_mdl,tmp_E3,tmp_crit]=mvar(x_sig,config_crit);
%         [tmp_pass,tmp_h,tmp_pVal]=test_model(tmp_E3,length(tmp_E3));
% 
%         if tmp_pass
%             break
%         end
%     end
% end

% %% Multivariate processing (one signal at a time)
%  
% for i=1:maxOrder
%     %%
%     config_crit.orderRange=i;
%     config_crit.output=0;
%     
%     e_mdl=struct;
%     e_E=zeros(length(E_orig)-i,numChannels);
% %     tmp_h=nan(numChannels,1);
% %     tmp_pVal=nan(numChannels,1);
%     
%     for j=1:numChannels
%     	[e_mdl(j).mdl,e_E(:,j),e_crit(i,j)]=mvar(E_orig(:,j),config_crit);
% %         [~,tmp_h(j),tmp_pVal(j)]=test_model(e_E(:,j),size(E_orig,1));
%     end
%     
%     [tmp_pass,tmp_h,tmp_pVal]=test_model(e_E,size(E_orig,1));
% %     tmp_pass=all(tmp_h==0);
%     
%     fprintf('%d, Overall: Pass = %d\n',i,tmp_pass);
%     fprintf('%d, Ch. 1: h = %d, p = %.2f\n',i,tmp_h(1),tmp_pVal(1));
%     fprintf('%d, Ch. 2: h = %d, p = %.2f\n',i,tmp_h(2),tmp_pVal(2));
%     fprintf('%d, Ch. 3: h = %d, p = %.2f\n',i,tmp_h(3),tmp_pVal(3));
%     fprintf('%d, Ch. 4: h = %d, p = %.2f\n',i,tmp_h(4),tmp_pVal(4));
%     
%     x_sig=zeros(length(x_sig_orig)-i,numChannels);
%     
%     for j=1:numChannels
%         x_sig(:,j)=filter_signal(x_sig_orig(:,j),e_mdl(j).mdl.AR);
%     end
%     
%     config_crit.orderRange=1:50;
%     config_crit.output=1;
%     
%     [tmp_mdl,tmp_E,tmp_crit]=mvar(x_sig,config_crit);
%     [tmp_pass,tmp_h,tmp_pVal]=test_model(tmp_E,length(tmp_E));
%     fprintf('%d, Overall: Pass = %d\n',i,tmp_pass);
%     fprintf('%d, Ch. 1: h = %d, p = %.2f\n',i,tmp_h(1),tmp_pVal(1));
%     fprintf('%d, Ch. 2: h = %d, p = %.2f\n',i,tmp_h(2),tmp_pVal(2));
%     fprintf('%d, Ch. 3: h = %d, p = %.2f\n',i,tmp_h(3),tmp_pVal(3));
%     fprintf('%d, Ch. 4: h = %d, p = %.2f\n',i,tmp_h(4),tmp_pVal(4));
%     
%     t=(0:length(x_sig_orig)-1)/fs;
%     m=tmp_mdl_orig.order;
%     
%     figure(hFig1);
%     subplot(311);
%     plot(t,x_sig_orig(:,1),'b',t(1+i:end),x_sig(:,1),'r'); title('Signal'); xlim([0 t(end)])
%     subplot(312);
%     plot(t(1+m:end),E_orig(:,1),'b',t(1+i+tmp_mdl.order:end),tmp_E(:,1),'r'); title('Residuals'); xlim([0 t(end)])
%     ax=subplot(313); cla;
%     plot_acf(tmp_E(:,1),[],[],ax); 
%     
%     figure(hFig2);
%     subplot(311);
%     plot(t,x_sig_orig(:,2),'b',t(1+i:end),x_sig(:,2),'r'); title('Signal'); xlim([0 t(end)])
%     subplot(312);
%     plot(t(1+m:end),E_orig(:,2),'b',t(1+i+tmp_mdl.order:end),tmp_E(:,2),'r'); title('Residuals'); xlim([0 t(end)])
%     ax=subplot(313); cla;
%     plot_acf(tmp_E(:,2),[],[],ax); 
%     
%     figure(hFig3);
%     subplot(311);
%     plot(t,x_sig_orig(:,3),'b',t(1+i:end),x_sig(:,3),'r'); title('Signal'); xlim([0 t(end)])
%     subplot(312);
%     plot(t(1+m:end),E_orig(:,3),'b',t(1+i+tmp_mdl.order:end),tmp_E(:,3),'r'); title('Residuals'); xlim([0 t(end)])
%     ax=subplot(313); cla;
%     plot_acf(tmp_E(:,3),[],[],ax); 
%     
%     figure(hFig4);
%     subplot(311);
%     plot(t,x_sig_orig(:,4),'b',t(1+i:end),x_sig(:,4),'r'); title('Signal'); xlim([0 t(end)])
%     subplot(312);
%     plot(t(1+m:end),E_orig(:,4),'b',t(1+i+tmp_mdl.order:end),tmp_E(:,4),'r'); title('Residuals'); xlim([0 t(end)])
%     ax=subplot(313); cla;
%     plot_acf(tmp_E(:,4),[],[],ax); 
%     waitforbuttonpress;
% end

% 
% %% Multivariate processing (all at once); not working
%  
% for i=1:maxOrder
%     %%
%     config_crit.orderRange=i;
% % errorModelOrder=5;
% % config_crit.orderRange=errorModelOrder;
%     config_crit.output=0;
%     
%     [e_mdl,e_E,e_crit(i)]=mvar(E_orig,config_crit);
% %     [e_mdl,e_E,e_crit]=mvar(E_orig,config_crit);
%     
%     [tmp_pass,tmp_h,tmp_pVal]=test_model(e_E,size(E_orig,1));
%     fprintf('%d, Overall: Pass = %d\n',i,tmp_pass);
%     fprintf('%d, Ch. 1: h = %d, p = %.2f\n',i,tmp_h(1),tmp_pVal(1));
%     fprintf('%d, Ch. 2: h = %d, p = %.2f\n',i,tmp_h(2),tmp_pVal(2));
%     fprintf('%d, Ch. 3: h = %d, p = %.2f\n',i,tmp_h(3),tmp_pVal(3));
%     fprintf('%d, Ch. 4: h = %d, p = %.2f\n',i,tmp_h(4),tmp_pVal(4));
%     
%     x_sig=filter_signal(x_sig_orig,e_mdl.AR);
%     
%     config_crit.orderRange=1:50;
%     config_crit.output=1;
%     
%     [tmp_mdl,tmp_E,tmp_crit]=mvar(x_sig,config_crit);
%     [tmp_pass,tmp_h,tmp_pVal]=test_model(tmp_E,length(tmp_E));
%     fprintf('%d, Overall: Pass = %d\n',i,tmp_pass);
%     fprintf('%d, Ch. 1: h = %d, p = %.2f\n',i,tmp_h(1),tmp_pVal(1));
%     fprintf('%d, Ch. 2: h = %d, p = %.2f\n',i,tmp_h(2),tmp_pVal(2));
%     fprintf('%d, Ch. 3: h = %d, p = %.2f\n',i,tmp_h(3),tmp_pVal(3));
%     fprintf('%d, Ch. 4: h = %d, p = %.2f\n',i,tmp_h(4),tmp_pVal(4));
%     
%     t=(0:length(x_sig_orig)-1)/fs;
%     m=tmp_mdl_orig.order;
% %     t2=(errorModelOrder:length(x_sig_orig)-1)/fs;
% %     t2_e=(errorModelOrder:size(e_E,1)-1+errorModelOrder)/fs;
% %     t2=(i:length(x_sig_orig)-1)/fs;
% %     t2_e=(i:size(e_E,1)-1+i)/fs;
% %     t2_e=(i:size(tmp_E,1)-1+i)/fs;
%     
%%
    figure(hFig1);
    subplot(311);
    plot(t,x_sig_orig(:,1),'b',t(1+i:end),x_sig(:,1),'r'); title('Signal'); xlim([0 t(end)])
    subplot(312);
%     plot(t1_e,E_orig(:,1),'b',t2_e,e_mdl.x_hat(:,1),'r'); title('Residuals'); xlim([0 t1_e(end)])
    plot(t(1+m:end),E_orig(:,1),'b',t(1+i+tmp_mdl.order:end),tmp_E(:,1),'r'); title('Residuals'); xlim([0 t(end)])
    ax=subplot(313); cla;
    plot_acf(tmp_E(:,1),[],[],ax); 
    
    figure(hFig2);
    subplot(311);
    plot(t,x_sig_orig(:,2),'b',t(1+i:end),x_sig(:,2),'r'); title('Signal'); xlim([0 t(end)])
    subplot(312);
%     plot(t1_e,E_orig(:,2),'b',t2_e,e_mdl.x_hat(:,2),'r'); title('Residuals'); xlim([0 t1_e(end)])
    plot(t(1+m:end),E_orig(:,2),'b',t(1+i+tmp_mdl.order:end),tmp_E(:,2),'r'); title('Residuals'); xlim([0 t(end)])
    ax=subplot(313); cla;
    plot_acf(tmp_E(:,2),[],[],ax); 
    
    figure(hFig3);
    subplot(311);
    plot(t,x_sig_orig(:,3),'b',t(1+i:end),x_sig(:,3),'r'); title('Signal'); xlim([0 t(end)])
    subplot(312);
%     plot(t1_e,E_orig(:,3),'b',t2_e,e_mdl.x_hat(:,3),'r'); title('Residuals'); xlim([0 t1_e(end)])
    plot(t(1+m:end),E_orig(:,3),'b',t(1+i+tmp_mdl.order:end),tmp_E(:,3),'r'); title('Residuals'); xlim([0 t(end)])
    ax=subplot(313); cla;
    plot_acf(tmp_E(:,3),[],[],ax); 
    
    figure(hFig4);
    subplot(311);
    plot(t,x_sig_orig(:,4),'b',t(1+i:end),x_sig(:,4),'r'); title('Signal'); xlim([0 t(end)])
    subplot(312);
%     plot(t1_e,E_orig(:,4),'b',t2_e,e_mdl.x_hat(:,4),'r'); title('Residuals'); xlim([0 t1_e(end)])
    plot(t(1+m:end),E_orig(:,4),'b',t(1+i+tmp_mdl.order:end),tmp_E(:,4),'r'); title('Residuals'); xlim([0 t(end)])
    ax=subplot(313); cla;
    plot_acf(tmp_E(:,4),[],[],ax); 
%     waitforbuttonpress;
% end

% %% Setting up for p=1, iteratively
% 
% x_sig_orig=x.Rest(:,:,currTrial);
% 
% config_crit.orderRange=ar.Rest(currTrial).mdl.order;
% config_crit.output=0;
% tmp_mdl_orig=ar.Rest(currTrial).mdl;
% E_orig=res.Rest(currTrial).E;
% 
% maxIterations=100;
% p=1;
% epsilon=2e-2;
% 
% E=E_orig;
% x_sig=x_sig_orig;
% 
% %% Main loop
% 
% for i=1:maxIterations
%     phi_hat=calculate_serial_coefficients(E,p);
%     
%     if all(all(abs(phi_hat) < epsilon))
%         fprintf('\nEstimation of phi converged\n');
%         break
%     end
%     
%     x_sig=filter_signal(x_sig,phi_hat);
%     [tmp_mdl,E,~]=mvar(x_sig,config_crit);
%     fprintf('%d,',i);
% end

%% Plot the original and new signal to compare



% t1_y=(0:size(tmp_mdl_orig.x_hat,1)-1)/fs;
% t2_y=(0:size(tmp_mdl.x_hat,1)-1)/fs;
% 
% t2_e=(0:size(E,1)-1)/fs;
% 
% for i=1:numChannels
%     figure;
%     subplot(311);
%     plot(t1,x_sig_orig(:,i),'b',t2,x_sig(:,i),'r'); legend('x','x\_new'); drawnow;
%     subplot(312);
%     plot(t1_y,tmp_mdl_orig.x_hat(:,i),'b',t2_y,tmp_mdl.x_hat(:,i),'r'); legend('y','y\_new'); 
%     subplot(313);
%     plot(t1_e,E_orig(:,i),'b',t2_e,E(:,i),'r'); legend('e','e\_new'); 
% end

%% Plot the ACF of the residuals before and after

% plot_acf(E_orig(:,1)); title('Original E - Ch1')
% plot_acf(E_orig(:,2)); title('Original E - Ch2')
% plot_acf(E_orig(:,3)); title('Original E - Ch3')
% plot_acf(E_orig(:,4)); title('Original E - Ch4')
% plot_acf(E(:,1)); title('New E - Ch1')
% plot_acf(E(:,2)); title('New E - Ch2')
% plot_acf(E(:,3)); title('New E - Ch3')
% plot_acf(E(:,4)); title('New E - Ch4')



