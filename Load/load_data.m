function [x,fs,x_all]=load_data(file,channels,condition,filtering)
%% [x,fs,x_all]=load_data(file,channels,condition,filtering)
%
%  Given a filename and the specific condition to take data from, returns the signal from
%  all the trials matching that condition, and the sampling frequency used for this file.
%  Filters the data before it is separated into its constitutent trials
%
%   Inputs:
%    - file: Filename specifying the full file path to the file
%    - condition: Assumed to be the condition in visual_stim, corresponding to different
%       conditions depending on the specific file opened. Ensure that the value given
%       matches the file given
%    - channels: Vector of ints defining which channels to extract. If it is a matrix, the
%       bipolar combination of channels is taken, with the second column channel subtracted
%       from the first column channel. Size is [c x 2], where c is the number of channels
%    - filtering: Struct containing optional additional filtering parameters
%       hpf: If defined, should be a struct containing num and den for a high pass filter. 
%           Default is to use a 4th order butterworth filter, with a cutoff of 1 Hz
%       comb: If defined, should be a struct containing num and den for a comb filter. 
%           Default is to use a 60 Hz comb filter, with a qFactor of 35. If this field is
%           empty, the default comb is used
%       lpf: If defined, should be a struct containing num and den for a low pass filter. 
%           No default
%       notch: If defined, should be a struct containing num and den for a notch filter. 
%           No default. Can contain more than one set of num and den for multiple notches
%       ma: If defined, should be an integer defining the number of samples to run
%           through the moving average
%       normalize: Separate from all other filtering techniques, the signal can be
%           normalized. String containing the method to normalize; 'none' performs no
%           additional normalization, 'z-score' normalizes the signals by dividing by the
%           standard deviation of the respective condition extracted (for x) or the
%           standard deviation of the entire signal
%       downsample: Separate from everything else, if defined, the signal will be
%           downsampled to match the sampling frequency given in this field. Note that
%           this value should be an even multiple of the original sampling frequency.
%           Additionally, it is recommended that this step is only performed if there is
%           no low-pass filtering done
%       realizations: Sub-struct with field 'length' defining the length of realizations
%           to split the trials into. Should be given as the number of samples to be taken
%           for each realization. If realizations is empty, default length is 1 second
%           times the sampling frequency. Note that this is the initial sampling
%           frequency, but the time of each realization in seconds will remain the same is
%           'downsample' is specified
%
%   Outputs:
%    - x: Matrix of values for all channels and all trials matching a particular
%       condition, with size [n x c x t], where n is the number of samples found in the
%       shortest length trial, c is the number of channels, and t is the number of trials
%       for that particular condition
%    - fs: Sampling frequency in Hz
%    - x_all: Optional output, which is the entire timeline of the signal for the entire
%       run
%
%  See also: varm, estimate, mvar, dtf, test_model, plot_connectivity
%

% example conditions: for ET_CL_04, 2018_06_20, run 5: 1 (rest), 2 (cue right), 3 (cue 
% left), 4 (move right), 5 (move left)

data=load(file);

fs=data.datastorage.src.LFP.Fs;

numChannels=size(channels,1);

x_all=zeros(size(data.datastorage.src.LFP.data,1)-1,numChannels);

if size(channels,2) > 1
    for i=1:numChannels
        x_all(:,i)=data.datastorage.src.LFP.data(1:end-1,channels(i,1))-data.datastorage.src.LFP.data(1:end-1,channels(i,2));
    end
else
    for i=1:numChannels
        x_all(:,i)=data.datastorage.src.LFP.data(1:end-1,channels(i));
    end
end

% Flexible filtering parameters from the filtering struct 

bool_normalize=false;
bool_realizations=false;
realizationLength=nan;

if nargin > 3 && isstruct(filtering)
    if isfield(filtering,'hpf')
        for i=1:numChannels
            x_all(:,i)=filtfilt(filtering.hpf.num,filtering.hpf.den,x_all(:,i));
        end
    else
        order_hp=4;
        cutoff_hp=1;
        [num_hp,den_hp]=CreateHPF_butter(fs,order_hp,cutoff_hp);

        for i=1:numChannels
            x_all(:,i)=filtfilt(num_hp,den_hp,x_all(:,i));
        end
    end
    
    if isfield(filtering,'comb')
        if isempty(filtering.comb)
            f0=60;
            w0=f0/(fs/2);
            qFactor=35;
            bw=w0/qFactor;
            [num_comb,den_comb]=iircomb(fs/f0,bw,'notch');

            for i=1:numChannels
                x_all(:,i)=filtfilt(num_comb,den_comb,x_all(:,i));
            end
        else
            num_comb=filtering.comb.num;
            den_comb=filtering.comb.den;

            for i=1:numChannels
                x_all(:,i)=filtfilt(num_comb,den_comb,x_all(:,i));
            end            
        end
    end

    if isfield(filtering,'lpf')
        for i=1:numChannels
            x_all(:,i)=filtfilt(filtering.lpf.num,filtering.lpf.den,x_all(:,i));
        end
    end

    if isfield(filtering,'notch')
        for i=1:length(filtering.notch)
            for j=1:numChannels
                x_all(:,j)=filtfilt(filtering.notch(i).num,filtering.notch(i).den,x_all(:,j));
            end
        end
    end

    if isfield(filtering,'ma')
        a=1;
        b=ones(filtering.ma,1)/filtering.ma;

        for j=1:numChannels
            x_all(:,j)=filter(b,a,x_all(:,j));
        end
    end
    
    if isfield(filtering,'normalize')
        if ~strcmp(filtering.normalize,'none')
            bool_normalize=true;
        end
    end
    
    if isfield(filtering,'realizations')
        bool_realizations=true;
        
        if isempty(filtering.realizations)
            realizationLength=1*fs; % 1 seconds worth of samples
        elseif isfield(filtering.realizations,'length')
            realizationLength=filtering.realizations.length;
        end
    end
end

% If no condition is given (condition is empty), return the whole signal

if isempty(condition) || ~any(condition)
    x=x_all;
    x_all=[];
    
    if isfield(filtering,'downsample')
        indToSkip=round(fs/filtering.downsample);
        x=downsample(x,indToSkip);
        fs=filtering.downsample;
    end

    if bool_normalize
        if strcmp(filtering.normalize,'z-score')
            for i=1:numChannels
                x(:,i) = x(:,i) / std(x(:,i));
            end
        end
    end
    
    return
end

% Extract all trials matching 'condition'

instruct=data.datastorage.src.visual_stim.data;

ind_change_all=find(diff(instruct) ~= 0)+1;
ind_change_curr=find(diff(instruct) ~= 0 & instruct(2:end) == condition)+1;

numTrials=length(ind_change_curr);

ind_curr_start=zeros(numTrials,1);
ind_curr_end=zeros(numTrials,1);

currTrial=1;

for i=1:length(ind_change_all)
    if instruct(ind_change_all(i)) == condition
        if i ~= length(ind_change_all)
            ind_curr_start(currTrial)=ind_change_all(i);
            ind_curr_end(currTrial)=ind_change_all(i+1);
        end
        currTrial=currTrial+1;
    end
end

if ind_curr_start(end) == 0
    ind_curr_start(end)=[];
    ind_curr_end(end)=[];
    numTrials=numTrials-1;
end

minLength=min(ind_curr_end-ind_curr_start);

if ~bool_realizations
    x=zeros(minLength,numChannels,numTrials);

    for i=1:numTrials
        x(:,:,i)=x_all(ind_curr_start(i):ind_curr_start(i)+minLength-1,:);
    end
    
    % Manually drop certain trials based on artifacts

    trialsToDrop=drop_trials(file,condition);
    x(:,:,trialsToDrop)=[];
    numTrials=numTrials-length(trialsToDrop);
else
    numRealizations=floor(minLength/realizationLength);
    
    x=zeros(realizationLength,numChannels,numTrials * numRealizations);
    
    currNum=1;
    
    for i=1:numTrials
        for j=1:numRealizations
            currStart=ind_curr_start(i)+(j-1)*realizationLength;
            currEnd=currStart+realizationLength-1;
            x(:,:,currNum)=x_all(currStart:currEnd,:);
            currNum=currNum+1;
        end
    end
    
    % Manually drop certain trials based on artifacts

    trialsToDrop=drop_trials(file,condition,numRealizations);
    x(:,:,trialsToDrop)=[];
    numTrials=numTrials-length(trialsToDrop);
end

% If filtering.normalize is defined, normalize the individual trials by the average std of
% the conditions, and normalize the overall signal by the std of the entire signal

if bool_normalize
    if strcmp(filtering.normalize,'z-score')
        for i=1:numChannels
            x_all(:,i) = x_all(:,i) / std(x_all(:,i));
        end

        avgStd=zeros(numChannels,1);

        for i=1:numTrials
            for j=1:numChannels
                avgStd(j)=avgStd(j)+std(x(:,j,i));
            end
        end

        avgStd=avgStd / numTrials;

        for i=1:numChannels
            x(:,i,:)=x(:,i,:)/avgStd(i);
        end
    end
end
    
if isfield(filtering,'downsample')
    indToSkip=round(fs/filtering.downsample);
    x_all=downsample(x_all,indToSkip);
    tmp_x=zeros(ceil(size(x,1)/indToSkip),size(x,2),size(x,3));
    
    for i=1:size(x,2)
        for j=1:size(x,3)
            tmp_x(:,i,j)=downsample(x(:,i,j),indToSkip);
        end
    end
    
    x=tmp_x;
    
    fs=filtering.downsample;
end

end