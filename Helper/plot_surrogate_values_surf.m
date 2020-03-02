function plot_surrogate_values_surf(surrogate,freqRange)
%% plot_surrogate_values_surf(surrogate,freqRange)
%
%  Plot all the values from surrogate in a 3D plot (see Ding et al. 2001 for examples).
%
%   Inputs:
%    - surrogate: Struct or matrix containing the surrogate values from all channels,
%       frequencies, and iterations. Size of each matrix should be [c x c x f x i], where
%       c is the number of channels, f is the number of frequencies, and i is the number
%       of iterations
%    - freqRange: Vector containing the values of the frequencies used
%
%   Outputs:
%    Figure containing the 3D plot of all surrogate values
%

if size(freqRange,1) < size(freqRange,2)
    freqRange=freqRange';
end

numChannels=size(surrogate,1);
numFrequencies=size(surrogate,3);
numIterations=size(surrogate,4);

iterValues=1:numIterations;

tmp_surr=squeeze(surrogate(1,2,:,:));

tmp_surr=sort(tmp_surr,2)';

figure;
surf(freqRange,iterValues,tmp_surr);

end