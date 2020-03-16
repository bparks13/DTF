function error=ste(x)
%% error=ste(x)
%
%  Calculates the standard error of the given signal by using the standard deviation
%  divided by the square root of the number of realizations in the calculation
%
%   Inputs:
%    - x: Matrix of values to calculate the standard error of. The last dimension is
%       assumed to be the realization dimension, and will be used to average across as
%       well as for counting the number of realizations for the standard error normalizing
%       factor of the square root of realizations
%
%   Outputs:
%    - error: Will be a matrix with one less dimension than the input x, and will contain
%       one standard error in each position
%

numDimensions=length(size(x));

numRealizations=size(x,numDimensions);

error=std(x,0,numDimensions)/sqrt(numRealizations);

end