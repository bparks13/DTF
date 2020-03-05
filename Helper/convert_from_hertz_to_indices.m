function freqBandInd=convert_from_hertz_to_indices(freqBands,freqForAnalysis)
%% freqBandInd=convert_from_hertz_to_indices(freqBands,freqForAnalysis)
%
%  Given the frequency bands of interest and the frequency values that were used for the
%  connectivity measurements, return the indices of the frequency values that correspond
%  to those bands of interest
%
%   Inputs:
%    - freqBands: 2 column matrix of the start and end values for the frequency bands of
%       interest. Size is [f x 2], where f is the number of frequency bands 
%    - freqForAnalysis: Vector of frequency values in Hz that were used to create the
%       connectivity values. Size is [1 x n], where n is the total number of frequencies
%       evaluated
%
%   Outputs:
%    - freqBandInd: Cell array of vectors that are the indices that relate freqBands to
%       freqForAnalysis.
%

numBands=size(freqBands,1);

freqBandInd=cell(numBands,1);

for i=1:numBands
    tmp_start=find(freqBands(i,1)==freqForAnalysis,1);
    tmp_end=find(freqBands(i,2)==freqForAnalysis,1);
    
    freqBandInd{i}=tmp_start:tmp_end;
end

end