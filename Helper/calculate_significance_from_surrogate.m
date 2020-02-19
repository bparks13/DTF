function significance=calculate_significance_from_surrogate(surrogate,alpha,type)
%% significance=calculate_significance_from_surrogate(surrogate,alpha,type)
%
%  This is a helper function to calculate the significance threshold based on the
%  surrogate values.
%
%   Inputs:
%    - surrogate: Matrix containing all of the surrogate analysis values for a particular
%       condition. Size should be [c x c x f x t], where c is the number of channels, f is
%       the number of frequencies being analyzed, and t is the number of iterations the
%       surrogate analysis was run for 
%    - alpha: Float defining what percentage of the distribution is considered to be a
%       significant amount of connection. Default is 0.01
%    - type: String defining the type of significance to run. Possible inputs include
%       'invariant', which does not separate by frequency, and 'dependent', which returns
%       a different significance value for each frequency. Default is 'invariant'
%
%   Outputs:
%    - significance: Depending on the 'type' specified, the size of significance can
%       change. For 'invariant', the size is [c x c], where c is the number of channels.
%       For 'dependent', size is [c x c x f], where f is the number of frequencies
%

numChannels=size(surrogate,1);
numFrequencies=size(surrogate,3);

if nargin==1
    alpha=0.01;
    type='invariant';
elseif nargin==2
    type='invariant';
end

if isstruct(surrogate)
    fields=fieldnames(surrogate);
    significance=struct;
    
    for i=1:length(fields)
        significance.(fields{i})=calculate_significance_from_surrogate(surrogate.(fields{i}),alpha,type);
    end
    
    return
end

if strcmp(type,'dependent')
    significance=ones(numChannels,numChannels,numFrequencies);

    for k=1:numChannels
        for l=1:numChannels
            if k ~= l
                for m=1:numFrequencies
                    orderedConnections=sort(squeeze(surrogate(k,l,m,:)));
                    significance(k,l,m)=orderedConnections(ceil(length(orderedConnections)*(1-alpha)));
                end
            end
        end
    end
elseif strcmp(type,'invariant')
    significance=ones(numChannels,numChannels);

    for k=1:numChannels
        for l=1:numChannels
            if k ~= l
                thisChannel=squeeze(surrogate(k,l,:,:));
                orderedConnections=sort(thisChannel(:));
                significance(k,l)=orderedConnections(ceil(length(orderedConnections)*(1-alpha)));
            end
        end
    end
end

end