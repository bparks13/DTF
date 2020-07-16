function [x,instruct]=load_data(datastorage,meta,filtering,datastorage_postacq,condition)
%% [x,instruct]=load_data(datastorage,meta,filtering,data_postacq,condition)
%
%  Given a filename and the specific condition to take data from, returns the signal from
%  all the trials matching that condition, and the sampling frequency used for this file.
%  Filters the data before it is separated into its constitutent trials
%
%   Inputs:
%    - datastorage: The datastorage struct in its entirety
%    - meta: Struct containing all relevant information for this run
%    -- settings: Struct defining settings relevant to the DTF being performed
%    --- cues_only: Boolean defining if the cues should be used by themselves (true), or
%         if there exists a postacquisition file with better resolution of movement epochs
%    --- extrap_method: String defining the method to use for extrapolation. Only used for
%         PC+S data. Inputs are ['linear', default], ['pchip'], and ['spline']
%    -- vars: Struct that contains the variables used in this run
%    --- channels: Depending on the inputs given, is either a vector of monopolar channels,
%         or a matrix of bipolar channels, with size [c x 2], where c is the number of
%         bipolar channels, and the second channel is subtracted from the first channel
%    --- labels: Cell array consisting of strings corresponding to the channels given back.
%         If monopolar channels are given, it lists the structure and number (i.e. Vim0), if
%         bipolar channels are given, it lists the structure and the order of subtraction
%        (i.e. Vim (3-2)).
%    --- conditions: Vector containing all possible conditions available for this specific
%         set of inputs. Size is [1 x m], where m is the number of conditions returned
%    --- cond_labels: Cell array consisting of strings denoting what each condition is
%    --- visit_type: String defining whether the run is 'intraop' or 'closed-loop'
%    --- postacq_type: String used as a regular expression of the order/numbering of conditions,
%         where the string R_CR_CL_MR_ML means there is Rest, Cue Right/Left, Move
%         Right/Left, and the corresponding condition values are 1,2,3,4,5, respectively. 
%    - filtering: Struct containing filtering parameters
%    -- hpf: If defined, should be a struct containing num and den for a high pass filter. 
%        Default is to use a 4th order butterworth filter, with a cutoff of 1 Hz
%    -- comb: If defined, should be a struct containing num and den for a comb filter. 
%        Default is to use a 60 Hz comb filter, with a qFactor of 35. If this field is
%        empty, the default comb is used
%    -- lpf: If defined, should be a struct containing num and den for a low pass filter. 
%        No default
%    -- notch: If defined, should be a struct containing num and den for a notch filter. 
%        No default. Can contain more than one set of num and den for multiple notches
%    -- ma: If defined, should be an integer defining the number of samples to run
%        through the moving average
%    -- normalize: Separate from all other filtering techniques, the signal can be
%        normalized. String containing the method to normalize; 'none' performs no
%        additional normalization, 'z-score' normalizes the signals by dividing by the
%        standard deviation of the respective condition extracted (for x) or the
%        standard deviation of the entire signal
%    -- downsample: Separate from everything else, if defined, the signal will be
%        downsampled to match the sampling frequency given in this field. Note that
%        this value should be an even multiple of the original sampling frequency.
%        Additionally, it is recommended that this step is only performed if there is
%        no low-pass filtering done
%    -- realizations: Sub-struct with field 'length' defining the length of realizations
%        to split the trials into. Should be given as the number of samples to be taken
%        for each realization. If realizations is empty, default length is 1 second
%        times the sampling frequency. Note that this is the initial sampling
%        frequency, but the time of each realization in seconds will remain the same is
%        'downsample' is specified
%    - datastorage_postacq: Struct containing the datastorage_postacq struct,
%       for aligning data according to acceleration
%    - condition: Integer value defining the state matching a specific condition
%
%   Outputs:
%    - x: Matrix of values for all channels and all trials matching a particular
%       condition, with size [n x c x t], where n is the number of samples found in the
%       shortest length trial, c is the number of channels, and t is the number of trials
%       for that particular condition
%    - instruct: Variable that contains either the cues for each event (cues_only == true)
%       or the time points of when the hands moved based on acceleration data (false)
%
%  See also: varm, estimate, mvar, dtf, test_model, plot_connectivity, load_variables,
%  extract_visual_stim_for_intraop, extract_visual_stim_for_closed_loop
%

% example conditions: for ET_CL_04, 2018_06_20, run 5: 1 (rest), 2 (cue right), 3 (cue 
% left), 4 (move right), 5 (move left)

%% HERE/STILL REWRITING THE HELP AND FUNCTION CALL ITSELF

if isfield(meta.vars,'visit_type')
    visit_type=meta.vars.visit_type;
else
    visit_type='';
end

if isfield(meta.settings,'cues_only')
    cues_only=meta.settings.cues_only;
else
    cues_only=true;
end

if isfield(meta.settings,'extrap_method') && ~isempty(meta.settings.extrap_method)
    extrap_method=meta.settings.extrap_method;
else
    extrap_method='linear';
end

if isfield(meta.vars,'postacq_type')
    postacq_type=meta.vars.postacq_type;
else
    error('postacq_type is a required field in the meta struct')
end

fs=datastorage.src.LFP.Fs;

channels=meta.vars.channels;
numChannels=size(channels,1);

x_all=zeros(size(datastorage.src.LFP.data,1)-1,numChannels);

if size(channels,2) > 1
    for i=1:numChannels
        x_all(:,i)=datastorage.src.LFP.data(1:end-1,channels(i,1))-datastorage.src.LFP.data(1:end-1,channels(i,2));
    end
else
    for i=1:numChannels
        x_all(:,i)=datastorage.src.LFP.data(1:end-1,channels(i));
    end
end

if strcmp(visit_type,'closed-loop')
    x_all=x_all*1e-3;    % Medtronic recordings are in µV, convert to mV to be consistent
end

% Check that all filtering parameters are stable

if ~check_stability(filtering)
    error('A filter is unstable. Check warnings for specifics.');
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

    if isfield(filtering,'notch')
        for i=1:length(filtering.notch)
            for j=1:numChannels
                x_all(:,j)=filtfilt(filtering.notch(i).num,filtering.notch(i).den,x_all(:,j));
            end
        end
    end

    if isfield(filtering,'lpf')
        for i=1:numChannels
            x_all(:,i)=filtfilt(filtering.lpf.num,filtering.lpf.den,x_all(:,i));
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

% Load the visual_stim (in whatever form) in case it is being returned for x_all

if strcmp(visit_type,'intraop') || strcmp(visit_type,'')
    instruct=extract_visual_stim_for_intraop(datastorage,datastorage_postacq,cues_only,postacq_type,realizationLength);
elseif strcmp(visit_type,'closed-loop')
    instruct=extract_visual_stim_for_closed_loop(datastorage_postacq,cues_only);
else
    error('Invalid visit_type variable given')
end

% If no condition is given (condition is empty), return the whole signal

if isempty(condition) || ~any(condition)
    x=x_all;
    
    if bool_normalize
        if strcmp(filtering.normalize,'z-score')
            for i=1:numChannels
                x(:,i) = x(:,i) / std(x(:,i));
            end
        end
    end
    
    if isfield(filtering,'downsample')
        freqRatio=round(fs/filtering.downsample);
    
        if freqRatio * filtering.downsample ~= fs % Uneven frequency ratio, need to interpolate data first
            t=(0:length(x)-1)/fs;
            t_new=linspace(0,t(end),round(t(end) * freqRatio * filtering.downsample));
            tmp_x=nan(length(t_new),numChannels);
            
            for i=1:numChannels
                tmp_x(:,i)=interp1(t,x(:,i),t_new,extrap_method); % Test different interpolation methods here
            end
            
            x=tmp_x;
        end
    
        x=downsample(x,freqRatio);
    end
    
    return
end

% Extract all trials matching 'condition' from instruct

[ind_curr_start,ind_curr_end,numTrials]=extract_trials_from_instruct(instruct,condition);

% ind_change_all=find(diff(instruct) ~= 0)+1;
% ind_change_curr=find(diff(instruct) ~= 0 & instruct(2:end) == condition)+1;
% 
% numTrials=length(ind_change_curr);
% 
% ind_curr_start=zeros(numTrials,1);
% ind_curr_end=zeros(numTrials,1);
% 
% currTrial=1;
% 
% for i=1:length(ind_change_all)
%     if instruct(ind_change_all(i)) == condition
%         if i ~= length(ind_change_all)
%             ind_curr_start(currTrial)=ind_change_all(i);
%             ind_curr_end(currTrial)=ind_change_all(i+1);
%         end
%         currTrial=currTrial+1;
%     end
% end
% 
% if ind_curr_start(end) == 0
%     ind_curr_start(end)=[];
%     ind_curr_end(end)=[];
%     numTrials=numTrials-1;
% end

if ~bool_realizations
    minLength=min(ind_curr_end-ind_curr_start);

    x=zeros(minLength,numChannels,numTrials);

    for i=1:numTrials
        x(:,:,i)=x_all(ind_curr_start(i):ind_curr_start(i)+minLength-1,:);
    end
    
    % Manually drop certain trials based on artifacts

    trialsToDrop=drop_trials(datastorage,condition);
    x(:,:,trialsToDrop)=[];
    numTrials=numTrials-length(trialsToDrop);
else
    numRealizations=floor((ind_curr_end-ind_curr_start)./realizationLength);
    
    x=zeros(realizationLength,numChannels,sum(numRealizations));
    
    currNum=1;
    
    for i=1:numTrials
        for j=1:numRealizations(i)
            currStart=ind_curr_start(i)+(j-1)*realizationLength;
            currEnd=currStart+realizationLength-1;
            x(:,:,currNum)=x_all(currStart:currEnd,:);
            currNum=currNum+1;
        end
    end
    
    % Manually drop certain trials based on artifacts

    trialsToDrop=drop_trials(datastorage,condition,sum(numRealizations));
    x(:,:,trialsToDrop)=[];
    numTrials=sum(numRealizations)-length(trialsToDrop);
end

% If filtering.normalize is defined, normalize the individual trials by the std of
% that realization

if bool_normalize
    if strcmp(filtering.normalize,'z-score')
        for i=1:numTrials
            for j=1:numChannels
                x(:,j,i)=x(:,j,i) / std(x(:,j,i));
            end
        end
    end
end
    
if isfield(filtering,'downsample')
    freqRatio=round(fs/filtering.downsample);
    
    if freqRatio * filtering.downsample ~= fs % Uneven frequency ratio, need to interpolate data first
        t=(0:length(x)-1)/fs;
        t_new=linspace(0,t(end),round(t(end) * freqRatio * filtering.downsample));
        tmp_x=nan(length(t_new),numChannels,size(x,3));

        for i=1:numChannels
            for j=1:size(x,3)
                tmp_x(:,i,j)=interp1(t,x(:,i,j),t_new,extrap_method); % Test different interpolation methods here
            end
        end
        
        x=tmp_x;
    end

    x=downsample(x,freqRatio);
end

if nargout == 1
    clear instruct;
end

end