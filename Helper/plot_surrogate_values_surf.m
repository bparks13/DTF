function plot_surrogate_values_surf(surrogate,alpha,freqRange)
%% plot_surrogate_values_surf(surrogate,alpha,freqRange)
%
%  Plot all the values from surrogate in a 3D plot (see Ding et al. 2001 for examples).
%
%   Inputs:
%    - surrogate: Struct or matrix containing the surrogate values from all channels,
%       frequencies, and iterations. Size of each matrix should be [c x c x f x i], where
%       c is the number of channels, f is the number of frequencies, and i is the number
%       of iterations
%    - alpha: Defines the value above which a connection would be considered significant;
%       given as a decimal corresponding to a percentage (i.e. 0.01 == 1%)
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

% iterValues=1:numIterations;

tmp_surr=squeeze(surrogate(1,2,:,:));
ordered=sort(tmp_surr(:));

values=[tmp_surr(:),repmat(freqRange,numIterations,1)];

dtf_values=0:0.005:round(ordered(round((1-alpha)*length(ordered)))*1.2,2);

centers={dtf_values, freqRange(1):1:freqRange(end)};

figure('renderer','opengl'); subplot(numChannels,numChannels,2);
hist3(values,centers); xlabel('DTF'); ylabel('Frequency [Hz]'); 
view(2);
xlim([dtf_values(1) dtf_values(end)]); ylim([freqRange(1), freqRange(end)]);
ax=gca; ax.Children(1).FaceColor='interp'; ax.Children(1).CDataMode='auto';

end