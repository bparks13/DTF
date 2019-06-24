%% testing_serial_correlation
%
%  This is a script to test the estimation of serially correlated errors
%
%   Based on The Econometric Analysis of Time Series, by Andrew Harvey
%

CCC;

%% Load

FILE='ET_CL_001__2017_05_17__run12__PSD__Z_SCORE.mat';
load(fullfile(get_root_path,'Files',FILE));
currTrial=1;
currChannel=3;
numChannels=size(x.Rest,2);

%% Testing for p=1, iteratively

x_sig_orig=x.Rest(:,:,currTrial);

config_crit.orderRange=ar.Rest(currTrial).mdl.order;
config_crit.output=0;
[tmp_mdl_orig,E_orig,~]=mvar(x_sig_orig,config_crit);

maxIterations=100;
epsilon=1e-2;

E=E_orig;
x_sig=x_sig_orig;

for i=1:maxIterations
    phi_hat=calculate_serial_coefficients(E,2);
    
    if all(all(abs(phi_hat) < epsilon))
        fprintf('\nEstimation of phi converged\n');
        break
    end
    
    x_sig=filter_signal(x_sig,phi_hat);
    [tmp_mdl,E,~]=mvar(x_sig,config_crit);
    fprintf('%d,',i);
end

%% Plot the original and new signal to compare

t=(0:length(x_sig_orig)-1)/fs;
t_y=(0:size(tmp_mdl.x_hat,1)-1)/fs;
t_e=(0:size(E,1)-1)/fs;

for i=1:numChannels
    figure;
    subplot(311);
    plot(t,x_sig_orig(:,i),'b',t,x_sig(:,i),'r'); legend('x','x\_new'); drawnow;
    subplot(312);
    plot(t_y,tmp_mdl_orig.x_hat(:,i),'b',t_y,tmp_mdl.x_hat(:,i),'r'); legend('y','y\_new'); 
    subplot(313);
    plot(t_e,E_orig(:,i),'b',t_e,E(:,i),'r'); legend('e','e\_new'); 
end

%% Plot the ACF of the residuals before and after

plot_acf(E_orig(:,1)); title('Original E - Ch1')
plot_acf(E_orig(:,2)); title('Original E - Ch2')
plot_acf(E_orig(:,3)); title('Original E - Ch3')
plot_acf(E_orig(:,4)); title('Original E - Ch4')
plot_acf(E(:,1)); title('New E - Ch1')
plot_acf(E(:,2)); title('New E - Ch2')
plot_acf(E(:,3)); title('New E - Ch3')
plot_acf(E(:,4)); title('New E - Ch4')




