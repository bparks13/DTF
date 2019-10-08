function x_shuffled=shuffle_time_series(x,shuffle_all)
%% x_shuffled=shuffle_time_series(x,shuffle_all)
%
%  Given a time-series matrix of values, returns the same number of
%  channels/trials/samples, but all samples within a trial have been randomly shuffled to
%  destroy temporal/spectral components.
%
%   Inputs:
%    - x: Matrix containing all time-series data. Size is [n x c x r], where n
%       is the number of samples, c is the number of channels, and r is the number of
%       realizations
%    - shuffle_all: Optional - Boolean denoting whether to shuffle all of the trials and
%       return a matrix of the same size as the input [true, default], or to return a
%       randomly selected subset of realizations that matches the number of channels given
%
%   Outputs:
%    - x_shuffled: Matrix containing all sample shuffled time-series data. Size is 
%       [n x c x r] for when shuffle_all is true, where n is the number of samples, c is
%       the number of channels, and r is the number of realizations. If shuffle_all is
%       false, the size is [n x c x c], where the number of realizations returned is equal
%       to the number of channels given
%
%  Examples:
%   x_shuffled = shuffle_time_series(x);
%   x_shuffled = shuffle_time_series(x, true);
%       - These will both return a matrix of the same size as the input 'x'
%
%   x_shuffled = shuffle_time_series(x, false);
%    	- This will return a smaller matrix, but is useful for surrogate analysis because
%    	it randomizes which trials are used each time, and randomizes the order of the
%    	samples each time, so the same shuffled trials are not reused
%   
%

if nargin == 1
    shuffle_all = true;
end

numSamples=size(x,1);
numChannels=size(x,2);
numTrials=size(x,3);
        
if shuffle_all
    x_shuffled=nan(size(x));
    
    for i=1:numTrials
        for j=1:numChannels
            tmp_x=x(:,j,i);
            
            for k=1:numSamples
                randomSample=randi(length(tmp_x),1);
                x_shuffled(k,j,i)=tmp_x(randomSample);
                tmp_x(randomSample)=[];
            end
        end
    end
else
    randomTrials=randperm(numTrials,numChannels);
    
    x_shuffled=nan(numSamples,numChannels);
    
    for i=1:numChannels
        tmp_x=x(:,i,randomTrials(i));
        
        for j=1:numSamples
            randomSample=randi(length(tmp_x),1);
            x_shuffled(j,i)=tmp_x(randomSample);
            tmp_x(randomSample)=[];
        end
    end
end

end