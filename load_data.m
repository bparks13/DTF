function [x,fs]=load_data(file,channels,condition,filtering)
%% [x,fs]=load_data(file,channels,condition,filtering)
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
%       hpf: If defined, should be a struct containing num and den. Default is to use a
%           4th order butterworth filter, with a cutoff of 1 Hz
%       comb: If defined, should be a struct containing num and den. Default is to use a
%           60 Hz comb filter, with a qFactor of 35
%       lpf: If defined, should be a struct containing num and den. No default
%       notch: If defined, should be a struct containing num and den. No default. Can
%           contain more than one set of num and den for multiple notches
%
%   Outputs:
%    - x: Matrix of values for all channels and all trials matching a particular
%       condition, with size [n x c x t], where n is the number of samples found in the
%       shortest length trial, c is the number of channels, and t is the number of trials
%       for that particular condition
%    - fs: Sampling frequency in Hz
%
%  See also: varm, estimate, mvar, dtf, test_model, plot_connectivity
%

% example conditions: for ET_CL_04, 2018_06_20, run 5: 1 (rest), 2 (cue right), 3 (cue 
% left), 4 (move right), 5 (move left)

data=load(file);

fs=data.datastorage.src.LFP.Fs;

tmp_x=zeros(size(data.datastorage.src.LFP.data,1)-1,length(channels));

numChannels=length(channels);

if size(channels,2) > 1
    for i=1:numChannels
        tmp_x(:,i)=data.datastorage.src.LFP.data(1:end-1,channels(i,1))-data.datastorage.src.LFP.data(1:end-1,channels(i,2));
    end
else
    for i=1:numChannels
        tmp_x(:,i)=data.datastorage.src.LFP.data(1:end-1,channels(i));
    end
end

order_hp=4;
cutoff_hp=2;
[num_hp,den_hp]=CreateHPF_butter(fs,order_hp,cutoff_hp);

for i=1:numChannels
    tmp_x(:,i)=filtfilt(num_hp,den_hp,tmp_x(:,i));
end

f0=60;
w0=f0/(fs/2);
qFactor=35;
bw=w0/qFactor;
[num_comb,den_comb]=iircomb(fs/f0,bw,'notch');

for i=1:numChannels
    tmp_x(:,i)=filtfilt(num_comb,den_comb,tmp_x(:,i));
end

% Flexible filtering parameters from the filtering struct 

if nargin > 3 && isstruct(filtering)
    if isfield(filtering,'lpf')
        for i=1:numChannels
            tmp_x(:,i)=filtfilt(filtering.lpf.num,filtering.lpf.den,tmp_x(:,i));
        end
    end
    
    if isfield(filtering,'notch')
        for i=1:length(filtering.notch)
            for j=1:numChannels
                tmp_x(:,j)=filtfilt(filtering.notch(i).num,filtering.notch(i).den,tmp_x(:,j));
            end
        end
    end
end

% If no condition is given (condition is empty), return the whole signal

if isempty(condition) || ~any(condition)
    x=tmp_x;
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

x=zeros(minLength,numChannels,numTrials);

for i=1:numTrials
    x(:,:,i)=tmp_x(ind_curr_start(i):ind_curr_start(i)+minLength-1,:);
end

end