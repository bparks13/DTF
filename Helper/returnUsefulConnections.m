function useful=returnUsefulConnections(freqBands,gamma,significance,combinations)
%% useful=returnUsefulConnections(freqBands,gamma,significance,combinations)
%
%  Helper function to reduce the number of combinations and only plot the bar plots of
%  combinations that have any significance in any frequency band
%
%   Inputs:
%    - freqBand: Indices pertaining to the frequencies in gamma to average over. Should be
%       a struct containing all the conditions in question
%    - gamma: Struct of matrices containing all of the connectivity values calculated for this
%       condition. Size is [d x d x f x t], where d is the number of channels, f is the
%       frequencies, and t is the number of trials
%    - significance: Struct of matrices defining what connectivity value is considered to be
%       significant. For 'invariant', the size is [c x c], where c is the number of
%       channels. For 'dependent', size is [c x c x f], where f is the number of 
%       frequencies. Significance type is automatically detected from the matrix size
%    - combinations: Matrix of indices pertaining to the channel combinations used in
%       gamma; size is [2*k x 2], where k is number of combinations given by the number of
%       channels (k = N!/(K!(N-K)!)). If this is not given, then the combinations are
%       automatically calculated
%
%   Outputs:
%    - useful: Returns combinations, but only with the significant combinations. Size is
%       dependant on how significant the connections are in each channel combination
%
%  See also: plot_bar_with_error.m
%

numFrequencyBands=length(freqBands);
fields=fieldnames(gamma);
numFields=length(fields);

if nargin==3
    numChannels=size(gamma.(fields{1}),1);
    combinations=[nchoosek(1:numChannels,2);nchoosek(numChannels:-1:1,2)];
end

numCombinations=size(combinations,1);
isSignificant=false(numCombinations,1);

for i=1:numFrequencyBands
    for j=1:numCombinations
        if ~isSignificant(j)
            for k=1:numFields
                tmp_mean=mean(mean(squeeze(gamma.(fields{k})(combinations(j,1),combinations(j,2),freqBands{i},:))));
                
                if tmp_mean > significance.(fields{k})(combinations(j,1),combinations(j,2))
                    isSignificant(j)=true;
                    break;
                end
            end
        end
    end
end

useful=combinations(isSignificant,:);

end