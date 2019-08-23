function [surrogate,distribution,avg_pxx]=surrogate_analysis(file,config)
%% [results,distribution,avg_pxx]=surrogate_analysis(file,config)
%
%  Performs surrogate analysis on a pre-existing file, creating an empirical distribution
%  of samples and then calculating if the values found for connectivity are significant
%  based on this distribution
%
%   Inputs:
%    - file: String containing the filename to be opened. If this is empty, uigetfile is
%       called so the user can pick the file to open. Note that this is not the absolute
%       path, only the filename itself
%    - config: Optional struct containing additional parameters
%       cond: A cell array containing one or more conditions to focus on, running
%           surrogate on each individual trial with no combining of trials
%       method: String indicating which method to use. All surrogate analysis uses
%           trial-shuffling, but this can be constrained to a single condition ['single']
%           (default), can shuffle trials from all conditions together ['combine'],
%           randomly shuffle individual time points from within each trial ['sample'], or
%           shift the whole signal in time ['shift']
%       signal: String defining if the original signal is used ['original', default], or the
%           decorrelated signal ['decorr']
%       iterations: Int defining the number of times to run through the dataset
%       numSamples: Int defining how many consecutive samples to randomly select if the
%           method is defined as 'sample'
%
%   Outputs:
%    - surrogate: Struct containing all of the values created by the surrogate analysis
%    - distribution: Struct containing the specific trials used in each condition, to
%       check that the surrogate analysis has a uniform distribution. 
%        * For the bivariate use case of 'shift', distribution is comprised of the random
%          trial chosen in the first column, followed by the random shift of the channel
%          in the second column
%    - avg_pxx: Optional output, returns the average PSD values of the randomly shuffled
%       trials
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

method='single';
bool_original=true;
bool_structGiven=false;

%% Load file

if nargin==0 || isempty(file)
    file=uigetfile(fullfile(get_root_path,'Files','*.mat'));
    
    if file == 0
        surrogate=[];
        return
    end
elseif isstruct(file)
    bool_structGiven=true;
end

if ~bool_structGiven
    data=load(fullfile(get_root_path,'Files',file));
else
    data=file;
    
    [~,name,ext] = fileparts(data.newFile);
    file=[name,ext];
end

fs=data.fs;
freqForAnalysis=data.config_plot.freqLims;
freqRange=data.freqRange;
fields=fieldnames(data.x);
numIterations=1000;
numSamplesUsed=nan;

if nargin == 2 && isstruct(config)
    if isfield(config,'cond') && ~isempty(config.cond)
        fields=config.cond;
    end
    
    if isfield(config,'method') && ~isempty(config.method)
        method=config.method;
    end
    
    if isfield(config,'signal') && ~isempty(config.signal)
        if strcmp(config.signal,'original')
            bool_original=true;
        elseif strcmp(config.signal,'decorr') && isfield(data,'x_filt')
            bool_original=false;
        end
    end
    
    if isfield(config,'iterations')
        numIterations=config.iterations;
    end
    
    if isfield(config,'numSamples') && ~isempty(config.numSamples)
        if strcmp(method,'sample')
            numSamplesUsed=config.numSamples;
        else
            warning('Cannot define the number of samples for any method other than ''sample''');
        end
    end
end

config_crit=struct(...
    'orderSelection',data.config_crit.orderSelection,...
    'crit',data.config_crit.crit,...
    'orderRange',[],...
    'method',data.config_crit.method,...
    'fs',fs,...
    'freqRange',freqForAnalysis,...
    'output',0);

surrogate=struct;
distribution=struct;

if bool_original
    xVar='x';
    arVar='ar';
else
    xVar='x_filt';
    arVar='ar_filt';
end

if nargout==3
    sum_pxx=[];
    avg_pxx=struct;
end
    
%% Surrogate analysis

if strcmp(method,'single')
    numFields=length(fields);
    totalOperations=numFields*numIterations;

    for k=1:numFields
        currCond=fields{k};

        x=data.(xVar).(currCond);
        ar=data.(arVar).(currCond);

        config_crit.orderRange=round(summarize_model_orders(ar));
        numSamples=size(x,1);
        numChannels=size(x,2);
        numTrials=size(x,3);

        tmp_x=zeros(numSamples,numChannels);
        gamma_dist=zeros(numChannels,numChannels,length(freqRange),numIterations);

        distribution.(currCond)=nan(numIterations,numChannels);

        for i=1:numIterations
            randomTrials=randperm(numTrials,numChannels);

            for j=1:numChannels
                tmp_x(:,j)=x(:,j,randomTrials(j));
            end

            [tmp_mdl,~,~]=mvar(tmp_x,config_crit);
            gamma_dist(:,:,:,i)=dtf(tmp_mdl,freqRange,fs);

            distribution.(currCond)(i,:)=randomTrials;

            if mod((k-1)*numIterations+i,floor(totalOperations/100)) == 0
                fprintf('%d%%\n',floor(((k-1)*numIterations+i)/totalOperations*100));
            end
        end

        surrogate.(currCond)=gamma_dist;
    end
elseif strcmp(method,'combine')
    minLength=inf;
    numTrials=0;
    numChannels=size(data.(xVar).(fields{1}),2);
    
    for i=1:length(fields)
        if minLength > size(data.(xVar).(fields{i}),1)
            minLength=size(data.(xVar).(fields{i}),1);
        end
        
        numTrials=numTrials+size(data.(xVar).(fields{i}),3);
    end
    
    x=nan(minLength,numChannels,numTrials);
    names=cell(numTrials,1);
    count=1;
    
    for i=1:length(fields)
        for j=1:size(data.(xVar).(fields{i}),3)
            x(:,:,count)=data.(xVar).(fields{i})(1:minLength,:,j);
            names{count}=sprintf('%s%d',fields{i},j);
            count=count+1;
        end
    end
    
    tmp_x=zeros(minLength,numChannels);
    gamma_dist=zeros(numChannels,numChannels,length(freqRange),numIterations);

    tmp_dist=cell(numIterations,numChannels);
    config_crit.orderRange=round(summarize_model_orders(data.(arVar)));

    for i=1:numIterations
        randomTrials=randperm(numTrials,numChannels);
        
        for j=1:numChannels
            tmp_x(:,j)=x(:,j,randomTrials(j));
            tmp_dist{i,j}=names{randomTrials(j)};
        end
        
        [tmp_mdl,~,~]=mvar(tmp_x,config_crit);
        gamma_dist(:,:,:,i)=dtf(tmp_mdl,freqRange,fs);

        if mod(i,floor(numIterations/100)) == 0
            fprintf('%d%%\n',floor(i/numIterations*100));
        end
    end
    
    for i=1:length(fields)
        surrogate.(fields{i})=gamma_dist;
        distribution.(fields{i})=tmp_dist;
    end
elseif strcmp(method,'sample')
    numFields=length(fields);
    totalOperations=numFields*numIterations;
    
    for k=1:numFields
        currCond=fields{k};

        x=data.(xVar).(currCond);
        ar=data.(arVar).(currCond);

        config_crit.orderRange=round(summarize_model_orders(ar));
        numSamples=size(x,1);
        numChannels=size(x,2);
        numTrials=size(x,3);
        
        if nargout == 3
            sum_pxx=zeros(length(freqRange),numChannels);
        end

        gamma_dist=zeros(numChannels,numChannels,length(freqRange),numIterations);
        distribution.(currCond)=nan(numIterations,numChannels);

        for i=1:numIterations
            randomTrials=randperm(numTrials,numChannels);

            tmp_x_shuffled=nan(numSamples,numChannels);
                
            numSampleIterations=floor(numSamples/numSamplesUsed);
            
            if numSamplesUsed == 1
                for l=1:numChannels
                    tmp_x=x(:,l,randomTrials(l));
                    for j=1:numSampleIterations
                        randomSample=randi(length(tmp_x),1);
                        tmp_x_shuffled(j,l)=tmp_x(randomSample);
                        tmp_x(randomSample)=[];
                    end
                end
                
                tmp_x=tmp_x_shuffled;
            else
                for l=1:numChannels
                    tmp_x=x(:,l,randomTrials(l));
                    for j=1:numSampleIterations
                        randomSample=randi(length(tmp_x)-numSamplesUsed+1,1);
                        tmp_x_shuffled((j-1)*numSamplesUsed+1:j*numSamplesUsed,l)=tmp_x(randomSample:randomSample+numSamplesUsed-1);
                        tmp_x(randomSample:randomSample+numSamplesUsed-1)=[];
                    end
                    
                    if numSampleIterations*numSamplesUsed < numSamples
                        tmp_x_shuffled(j*numSamplesUsed+1:end,l)=tmp_x;
                    end
                end
                
                tmp_x=tmp_x_shuffled;
            end

            [tmp_mdl,~,~]=mvar(tmp_x,config_crit);
            gamma_dist(:,:,:,i)=dtf(tmp_mdl,freqRange,fs);

            distribution.(currCond)(i,:)=randomTrials;
            
            if nargout == 3
                sum_pxx=sum_pxx+periodogram(tmp_x,[],freqRange,fs);
            end

            if mod((k-1)*numIterations+i,floor(totalOperations/100)) == 0
                fprintf('%d%%\n',floor(((k-1)*numIterations+i)/totalOperations*100));
            end
        end
        
        if nargout == 3
            avg_pxx.(currCond)=sum_pxx/numIterations;
        end

        surrogate.(currCond)=gamma_dist;
    end
elseif strcmp(method,'shift')
    numFields=length(fields);
    
    for k=1:numFields
        currCond=fields{k};

        x=data.(xVar).(currCond);
        ar=data.(arVar).(currCond);

        config_crit.orderRange=round(summarize_model_orders(ar));
        numSamples=size(x,1);
        numChannels=size(x,2);
        numTrials=size(x,3);
        
        maxShift=numSamples-1;

        totalOperations=numFields*numIterations; 
        
        gamma_dist=zeros(numChannels,numChannels,length(freqRange),numIterations);
        distribution.(currCond)=nan(numIterations,numChannels);

        if numChannels == 2 % Bivariate case is different
            for i=1:numIterations
                randomTrial=randi(numTrials,1);
                randomShift=randi(maxShift,1);
                
                tmp_x=nan(numSamples,numChannels);
                tmp_x(:,1)=x(:,1,randomTrial);
                tmp_x(:,2)=wshift('1',x(:,2,randomTrial),randomShift);
                
                [tmp_mdl,~,~]=mvar(tmp_x,config_crit);
                gamma_dist(:,:,:,i)=dtf(tmp_mdl,freqRange,fs);

                distribution.(currCond)(i,:)=[randomTrial,randomShift];

                if mod((k-1)*numIterations+i,floor(totalOperations/100)) == 0
                    fprintf('%d%%\n',floor(((k-1)*numIterations+i)/totalOperations*100));
                end
            end

            surrogate.(currCond)=gamma_dist;
        end
    end
end

%% Either save the file, or return the struct

if nargout==0
    if bool_original
        save(fullfile(get_root_path,'Files',file),'-append','surrogate','distribution')
        clear surrogate distribution
    else
        surrogate_filt=surrogate; %#ok<*NASGU>
        distribution_filt=distribution;
        save(fullfile(get_root_path,'Files',file),'-append','surrogate_filt','distribution_filt')
        clear surrogate_filt surrogate distribution distribution_filt
    end
end

end