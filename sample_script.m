%% sample_script
%   
%  Loads one file from the server, and extracts a single trial to run in order to
%  highlight the workflow of the function calls
%

clear; close all; clc;

%% Load variables from the server

[x,fs]=load_data_sample();

%% Calculate MVAR model

[estMdl,E,BIC]=mvar(x);

%% Test whiteness

h=test_model(E,x);

if h
    disp('WARNING: Null hypothesis of uncorrelated errors rejected.')
end

%% Calculate DTF

freqRange=2:100;

gamma=dtf(estMdl,freqRange,fs);

%% Plot

labels={'Vim','Cortex'};
plot_connectivity(gamma,freqRange,labels)