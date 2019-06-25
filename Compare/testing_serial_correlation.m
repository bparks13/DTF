%% testing_serial_correlation
%
%  This is a script to test the estimation of serially correlated errors
%
%   Based on The Econometric Analysis of Time Series, by Andrew Harvey
%

CCC;

%% Load

FILE='ET_CL_004__2018_06_20__run5__200Hz__Z_SCORE__BIC_(1).mat';
load(fullfile(get_root_path,'Files',FILE));
currTrial=1;
currChannel=1;
numChannels=size(x.Rest,2);

%% Setting up for p=1, iteratively

x_sig_orig=x.Rest(:,:,currTrial);

config_crit.orderRange=ar.Rest(currTrial).mdl.order;
config_crit.output=0;
% [tmp_mdl_orig,E_orig,~]=mvar(x_sig_orig,config_crit);
tmp_mdl_orig=ar.Rest(currTrial).mdl;
E_orig=res.Rest(currTrial).E;

maxIterations=100;
p=1;
epsilon=2e-2;

E=E_orig;
x_sig=x_sig_orig;

%% Main loop

for i=1:maxIterations
    phi_hat=calculate_serial_coefficients(E,p);
    
    if all(all(abs(phi_hat) < epsilon))
        fprintf('\nEstimation of phi converged\n');
        break
    end
    
    x_sig=filter_signal(x_sig,phi_hat);
    [tmp_mdl,E,~]=mvar(x_sig,config_crit);
    fprintf('%d,',i);
end

%% Plot the original and new signal to compare

t1=(0:length(x_sig_orig)-1)/fs;
t2=(0:length(x_sig)-1)/fs;
t1_y=(0:size(tmp_mdl_orig.x_hat,1)-1)/fs;
t2_y=(0:size(tmp_mdl.x_hat,1)-1)/fs;
t1_e=(0:size(E_orig,1)-1)/fs;
t2_e=(0:size(E,1)-1)/fs;

for i=1:numChannels
    figure;
    subplot(311);
    plot(t1,x_sig_orig(:,i),'b',t2,x_sig(:,i),'r'); legend('x','x\_new'); drawnow;
    subplot(312);
    plot(t1_y,tmp_mdl_orig.x_hat(:,i),'b',t2_y,tmp_mdl.x_hat(:,i),'r'); legend('y','y\_new'); 
    subplot(313);
    plot(t1_e,E_orig(:,i),'b',t2_e,E(:,i),'r'); legend('e','e\_new'); 
end

%% Plot the ACF of the residuals before and after

% plot_acf(E_orig(:,1)); title('Original E - Ch1')
% plot_acf(E_orig(:,2)); title('Original E - Ch2')
% plot_acf(E_orig(:,3)); title('Original E - Ch3')
% plot_acf(E_orig(:,4)); title('Original E - Ch4')
% plot_acf(E(:,1)); title('New E - Ch1')
% plot_acf(E(:,2)); title('New E - Ch2')
% plot_acf(E(:,3)); title('New E - Ch3')
% plot_acf(E(:,4)); title('New E - Ch4')




