CCC;

%% Load

load(fullfile(get_root_path,'Files','S01_D01_R01.mat'));
% load(fullfile(get_root_path,'Files','S01_D01_R04.mat'));

I=cell(6,1);
sig=cell(6,1);
alpha=[0.01,0.05,0.01,0.05,0.01,0.05]';

minIterationsPerRealization=[50,100];

%% Create multiple significance thresholds

I{1}=struct('Rest',1000,'MoveRight',1000,'MoveLeft',1000,'FeelRight',1000);
I{2}=struct('Rest',1000,'MoveRight',1000,'MoveLeft',1000,'FeelRight',1000);

% I{1}=struct('Rest',1000,'MoveRight',1000,'MoveLeft',1000,'ImagineRight',1000,'ImagineLeft',1000);
% I{2}=struct('Rest',1000,'MoveRight',1000,'MoveLeft',1000,'ImagineRight',1000,'ImagineLeft',1000);

sig{1}=calculate_significance_from_surrogate(surrogate_filt,alpha(1),'invariant');
sig{2}=calculate_significance_from_surrogate(surrogate_filt,alpha(2),'invariant');

config_surr.iterations=1000;
config_surr=check_number_of_iterations(config_surr,x_filt,minIterationsPerRealization(1));
[tmp_surr1,~]=surrogate_analysis(x_filt,fs,freqForAnalysis,config_mvar,config_surr); % I=Ind.(50)
sig{3}=calculate_significance_from_surrogate(tmp_surr1,alpha(3),'invariant');
sig{4}=calculate_significance_from_surrogate(tmp_surr1,alpha(4),'invariant');

I{3}=config_surr.iterations;
I{4}=config_surr.iterations;

config_surr.iterations=1000;
config_surr=check_number_of_iterations(config_surr,x_filt,minIterationsPerRealization(2));
[tmp_surr2,~]=surrogate_analysis(x_filt,fs,freqForAnalysis,config_mvar,config_surr); % I=Ind.(100)
sig{5}=calculate_significance_from_surrogate(tmp_surr2,alpha(5),'invariant');
sig{6}=calculate_significance_from_surrogate(tmp_surr2,alpha(6),'invariant');

I{5}=config_surr.iterations;
I{6}=config_surr.iterations;

%% Plot in different figures for the different conditions

hFig=struct;

fields=fieldnames(x_filt);
numFields=length(fields);

numChannels=4;

colors=linspecer(6);

for j=1:numFields
    hFig.(fields{j})=figure('Name',fields{j});

    for k=1:numChannels
        for l=1:numChannels
            if k ~= l
                subplot(numChannels,numChannels,(k-1)*numChannels+l);
                legendNames=cell(6,1);

                for i=1:6
                    plot([0,1],[sig{i}.(fields{j})(k,l),sig{i}.(fields{j})(k,l)],'Color',colors(i,:),'LineWidth',2); 
                    hold on;
                    legendNames{i}=sprintf('%d,%.2f',I{i}.(fields{j}),alpha(i));
                end
                
                legend(legendNames); ylim([0 1]);
            end
        end
    end
    
    drawnow
end



