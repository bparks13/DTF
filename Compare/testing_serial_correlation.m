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
% currChannel=3;
numChannels=size(x.Rest,2);

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

%% Testing different model orders for estimating AR coefficients for serially correlated errors

config_crit.output=0;
maxOrder=15;
e_crit=zeros(maxOrder,1);
% x_sig_orig=x.Rest(:,currChannel,currTrial);
x_sig_orig=x.Rest(:,:,currTrial);
tmp_mdl_orig=ar.Rest(currTrial).mdl;
% E_orig=res.Rest(currTrial).E(:,currChannel);
E_orig=res.Rest(currTrial).E;

t1=(0:size(x_sig_orig,1)-1)/fs;
t1_e=(0:size(E_orig,1)-1)/fs;

figure;
ax=subplot(221);
plot_acf(E_orig(:,1),[],[],ax); title('Ch. 1')
ax=subplot(222);
plot_acf(E_orig(:,2),[],[],ax); title('Ch. 2')
ax=subplot(223);
plot_acf(E_orig(:,3),[],[],ax); title('Ch. 3')
ax=subplot(224);
plot_acf(E_orig(:,4),[],[],ax); title('Ch. 4')

hFig1=figure('Name','Ch. 1');
hFig2=figure('Name','Ch. 2');
hFig3=figure('Name','Ch. 3');
hFig4=figure('Name','Ch. 4');

%%
 
for i=1:maxOrder
    %%
    config_crit.orderRange=i;
% errorModelOrder=5;
% config_crit.orderRange=errorModelOrder;
    config_crit.output=0;
    
    [e_mdl,e_E,e_crit(i)]=mvar(E_orig,config_crit);
%     [e_mdl,e_E,e_crit]=mvar(E_orig,config_crit);
    
    [tmp_pass,tmp_h,tmp_pVal]=test_model(e_E,size(E_orig,1));
    fprintf('%d, Overall: Pass = %d\n',i,tmp_pass);
    fprintf('%d, Ch. 1: h = %d, p = %.2f\n',i,tmp_h(1),tmp_pVal(1));
    fprintf('%d, Ch. 2: h = %d, p = %.2f\n',i,tmp_h(2),tmp_pVal(2));
    fprintf('%d, Ch. 3: h = %d, p = %.2f\n',i,tmp_h(3),tmp_pVal(3));
    fprintf('%d, Ch. 4: h = %d, p = %.2f\n',i,tmp_h(4),tmp_pVal(4));
    
    x_sig=filter_signal(x_sig_orig,e_mdl.AR);
    
    config_crit.orderRange=1:50;
    config_crit.output=1;
    
    [tmp_mdl,tmp_E,tmp_crit]=mvar(x_sig,config_crit);
    [tmp_pass,tmp_h,tmp_pVal]=test_model(tmp_E,length(tmp_E));
    fprintf('%d, Overall: Pass = %d\n',i,tmp_pass);
    fprintf('%d, Ch. 1: h = %d, p = %.2f\n',i,tmp_h(1),tmp_pVal(1));
    fprintf('%d, Ch. 2: h = %d, p = %.2f\n',i,tmp_h(2),tmp_pVal(2));
    fprintf('%d, Ch. 3: h = %d, p = %.2f\n',i,tmp_h(3),tmp_pVal(3));
    fprintf('%d, Ch. 4: h = %d, p = %.2f\n',i,tmp_h(4),tmp_pVal(4));
    
    t=(0:length(x_sig_orig)-1)/fs;
    m=tmp_mdl_orig.order;
%     t2=(errorModelOrder:length(x_sig_orig)-1)/fs;
%     t2_e=(errorModelOrder:size(e_E,1)-1+errorModelOrder)/fs;
%     t2=(i:length(x_sig_orig)-1)/fs;
%     t2_e=(i:size(e_E,1)-1+i)/fs;
%     t2_e=(i:size(tmp_E,1)-1+i)/fs;
    
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
    waitforbuttonpress;
end

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




