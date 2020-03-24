function [surrogate,distribution,pxx]=surrogate_analysis(x,fs,freqRange,config_crit,config)
%% [results,distribution,pxx]=surrogate_analysis(x,fs,freqRange,config_crit,config)
%
%  Performs surrogate analysis on the given data, creating an empirical distribution
%  of samples and then calculating if the values found for connectivity are significant
%  based on this distribution
%
%   Inputs:
%    - x: Struct containing the original time series data. The first fields are the
%       conditions from the runs, and each field is a matrix of size is [n x c x r], where
%       n is the number of samples, c is the number of channels, and r is the number of
%       realizations.
%    - fs: Sampling frequency in Hz
%    - freqRange: Vector of the range of frequencies over which the connectivity is
%       measured 
%    - config_crit: Struct containing additional parameters, as defined in mvar.m
%    - config: Struct containing additional parameters
%       cond: Optional, a cell array containing one or more conditions to focus on, running
%           surrogate on each individual trial with no combining of trials
%       method: String indicating which method to use. All surrogate analysis uses
%           trial-shuffling, but this can be constrained to a single condition ['single'],
%           can shuffle trials from all conditions together ['combine'],
%           randomly shuffle individual time points from within each trial ['sample',
%           (default)], or shift the whole signal in time ['shift']
%       iterations: Int defining the number of times to run through the dataset
%       numSamples: Int defining how many consecutive samples to randomly select if the
%           method is defined as 'sample'. Default is 1
%
%   Outputs:
%    - surrogate: Struct containing all of the values created by the surrogate analysis
%    - distribution: Struct containing the specific trials used in each condition, to
%       check that the surrogate analysis has a uniform distribution. 
%    - pxx: Optional output, returns a struct containing all of the PSD values for each
%       iteration, as well as the averages and standard deviations
%  
% 'sample' based shuffling is done by randomly choosing a trial, and then randomly
% shuffling the time points, independently, for that channel, with a different
% trial for each channel
%
% 'shift' based shuffling is different for bivariate vs. multivariate signals. For
% bivariate cases, one channel is kept constant, and the other channel is systematically
% shifted one sample at a time. The number of permutations is then equal to the number of
% samples. For multivariate csaes, all channels are independently shifted, and the number
% of iterations is used here.
%
%  See also: mvar, dtf
%

%  Based on the trial-shuffling procedure for surrogate found in
%  10.1016/j.neuroimage.2004.09.036 

%% Load file

fields=fieldnames(x);
numFields=length(fields);
numIterations=0;
numSamplesUsed=1;
method='sample';

if isstruct(config)
    if isfield(config,'cond') && ~isempty(config.cond)
        fields=config.cond;
    end
    
    if isfield(config,'method') && ~isempty(config.method)
        method=config.method;
    end
    
    if isfield(config,'iterations') && ~isempty(config.iterations)
        if isstruct(config.iterations)
            numIterations=cell(numFields,1);
            
            for i=1:numFields
                numIterations{i}=config.iterations.(fields{i});
            end
        else
            numIterations=cell(numFields,1);
            
            for i=1:numFields
                numIterations{i}=config.iterations;
            end
        end
    end
    
    if isfield(config,'numSamples') && ~isempty(config.numSamples)
        if strcmp(method,'sample')
            numSamplesUsed=config.numSamples;
        else
            warning('Cannot define the number of samples for any method other than ''sample''');
        end
    end
end

if ~iscell(numIterations)
    numIterations=cell(numFields,1);
    
    for i=1:numFields
        numIterations{i}=1000;
    end
end

config_crit.output=0;

if isfield(config_crit,'modelOrder')
    config_crit=rmfield(config_crit,'modelOrder');
end

surrogate=struct;
distribution=struct;

if nargout == 3
    pxx=struct;
end
    
%% Surrogate analysis

if strcmp(method,'single')
    totalOperations=numFields*numIterations;

    for k=1:numFields
        currCond=fields{k};

        numSamples=size(x.(currCond),1);
        numChannels=size(x.(currCond),2);
        numTrials=size(x.(currCond),3);

        tmp_x=zeros(numSamples,numChannels);
        gamma_dist=zeros(numChannels,numChannels,length(freqRange),numIterations);

        distribution.(currCond)=nan(numIterations,numChannels);
        
        if nargout == 3
            pxx.(currCond)=struct('all',nan(length(freqRange),numChannels,numIterations),...
                'avg',nan(length(freqRange),numChannels),'std',nan(length(freqRange),numChannels));
        end

        for i=1:numIterations
            randomTrials=randperm(numTrials,numChannels);

            for j=1:numChannels
                tmp_x(:,j)=x.(currCond)(:,j,randomTrials(j));
            end

            [tmp_mdl,~,~]=mvar(tmp_x,config_crit);
            gamma_dist(:,:,:,i)=dtf(tmp_mdl,freqRange,fs);

            distribution.(currCond)(i,:)=randomTrials;
            
            if nargout == 3
                pxx.(currCond).all(:,:,i)=periodogram(tmp_x,[],freqRange,fs);
            end

            if mod((k-1)*numIterations+i,floor(totalOperations/100)) == 0
                fprintf('%d%%\n',floor(((k-1)*numIterations+i)/totalOperations*100));
            end
        end

        surrogate.(currCond)=gamma_dist;
        
        if nargout == 3
            pxx.(currCond).avg=mean(pxx.(currCond).all,3);
            pxx.(currCond).std=std(pxx.(currCond).all,[],3);
        end
    end
elseif strcmp(method,'combine')
    minLength=inf;
    numTrials=0;
    numChannels=size(x.(fields{1}),2);
    
    for i=1:length(fields)
        if minLength > size(x.(fields{i}),1)
            minLength=size(x.(fields{i}),1);
        end
        
        numTrials=numTrials+size(x.(fields{i}),3);
    end
    
    x_all=nan(minLength,numChannels,numTrials);
    names=cell(numTrials,1);
    count=1;
    
    for i=1:length(fields)
        for j=1:size(x.(fields{i}),3)
            x_all(:,:,count)=x.(fields{i})(1:minLength,:,j);
            names{count}=sprintf('%s%d',fields{i},j);
            count=count+1;
        end
    end
    
    tmp_x=zeros(minLength,numChannels);
    gamma_dist=zeros(numChannels,numChannels,length(freqRange),numIterations);

    tmp_dist=cell(numIterations,numChannels);

    if nargout == 3
        pxx=struct('all',nan(length(freqRange),numChannels,numIterations),...
            'avg',nan(length(freqRange),numChannels),'std',nan(length(freqRange),numChannels));
    end

    for i=1:numIterations
        randomTrials=randperm(numTrials,numChannels);
        
        for j=1:numChannels
            tmp_x(:,j)=x_all(:,j,randomTrials(j));
            tmp_dist{i,j}=names{randomTrials(j)};
        end
        
        [tmp_mdl,~,~]=mvar(tmp_x,config_crit);
        gamma_dist(:,:,:,i)=dtf(tmp_mdl,freqRange,fs);

        if nargout == 3
            pxx.all(:,:,i)=periodogram(tmp_x,[],freqRange,fs);
        end

        if mod(i,floor(numIterations/100)) == 0
            fprintf('%d%%\n',floor(i/numIterations*100));
        end
    end
    
    for i=1:length(fields)
        surrogate.(fields{i})=gamma_dist;
        distribution.(fields{i})=tmp_dist;
        
        if nargout == 3 
            pxx.(fields{i}).avg=mean(pxx.all,3);
            pxx.(fields{i}).std=std(pxx.all,[],3);
        end
    end
    
    pxx=rmfield(pxx,'all');
elseif strcmp(method,'sample')
    x_cell=struct2cell(x);
    dist_cell=cell(numFields,1);
    surr_cell=cell(numFields,1);
    
    parfor k=1:numFields
%     for k=1:numFields
        numSamples=size(x_cell{k},1);
        numChannels=size(x_cell{k},2);
        numTrials=size(x_cell{k},3);
        
        gamma_dist=zeros(numChannels,numChannels,length(freqRange),numIterations{k});
        dist_cell{k}=nan(numIterations{k},numChannels);

        for i=1:numIterations{k}
            randomTrials=randperm(numTrials,numChannels);

            tmp_x_shuffled=nan(numSamples,numChannels);
                
            numSampleIterations=floor(numSamples/numSamplesUsed);
            
            if numSamplesUsed == 1
                for m=1:numChannels
                    tmp_x=x_cell{k}(:,m,randomTrials(m));
                    randPerm=randperm(size(tmp_x,1));
                    
                    tmp_x_shuffled(:,m)=tmp_x(randPerm); % Randomly shuffle the time domain data
                end
                
                tmp_x=tmp_x_shuffled;
            else
                for m=1:numChannels
                    tmp_x=x_cell{k}(:,m,randomTrials(m));
                    for j=1:numSampleIterations
                        randomSample=randi(length(tmp_x)-numSamplesUsed+1,1);
                        tmp_x_shuffled((j-1)*numSamplesUsed+1:j*numSamplesUsed,m)=tmp_x(randomSample:randomSample+numSamplesUsed-1);
                        tmp_x(randomSample:randomSample+numSamplesUsed-1)=[];
                    end
                    
                    if numSampleIterations*numSamplesUsed < numSamples
                        tmp_x_shuffled(j*numSamplesUsed+1:end,m)=tmp_x;
                    end
                end
                
                tmp_x=tmp_x_shuffled;
            end

            [tmp_mdl,~,~]=mvar(tmp_x,config_crit);
            gamma_dist(:,:,:,i)=dtf(tmp_mdl,freqRange,fs);

            dist_cell{k}(i,:)=randomTrials;
            
            if mod(i,floor(numIterations{k}/20)) == 0
                fprintf('%s - %d%%\n',fields{k},floor(i/numIterations{k}*100));
            end
        end
        
        surr_cell{k}=gamma_dist;
    end
    
    distribution=cell2struct(dist_cell,fields,1);
    surrogate=cell2struct(surr_cell,fields,1);
elseif strcmp(method,'shift')
    for k=1:numFields
        currCond=fields{k};

        numSamples=size(x.(currCond),1);
        numChannels=size(x.(currCond),2);
        numTrials=size(x.(currCond),3);
        
        maxShift=numSamples-1;

        totalOperations=numFields*numIterations; 
        
        gamma_dist=zeros(numChannels,numChannels,length(freqRange),numIterations);
        distribution.(currCond)=nan(numIterations,numChannels);

        if nargout == 3
            pxx.(currCond)=struct('all',nan(length(freqRange),numChannels,numIterations),...
                'avg',nan(length(freqRange),numChannels),'std',nan(length(freqRange),numChannels));
        end

        if numChannels == 2 % Bivariate case is different
            for i=1:numIterations
                randomTrial=randi(numTrials,1);
                randomShift=randi(maxShift,1);
                
                tmp_x=nan(numSamples,numChannels);
                tmp_x(:,1)=x.(currCond)(:,1,randomTrial);
                tmp_x(:,2)=wshift('1',x.(currCond)(:,2,randomTrial),randomShift);
                
                [tmp_mdl,~,~]=mvar(tmp_x,config_crit);
                gamma_dist(:,:,:,i)=dtf(tmp_mdl,freqRange,fs);

                distribution.(currCond)(i,:)=[randomTrial,randomShift];

                if nargout == 3
                    pxx.(currCond).all(:,:,i)=periodogram(tmp_x,[],freqRange,fs);
                end

                if mod((k-1)*numIterations+i,floor(totalOperations/100)) == 0
                    fprintf('%d%%\n',floor(((k-1)*numIterations+i)/totalOperations*100));
                end
            end
            
            if nargout == 3
                pxx.(currCond).avg=mean(pxx.(currCond).all,3);
                pxx.(currCond).std=std(pxx.(currCond).all,[],3);
            end

            surrogate.(currCond)=gamma_dist;
        else
            disp('Multivariate time shifting is not implemented yet')
            return
        end
    end
end

end