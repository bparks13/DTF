function [pxx]=plot_psd(in,fs,config)
%% [pxx]=plot_psd(in,fs,config)
%
%  Given either a real signal or the AR coefficients, plot the estimated Power Spectral
%  Density
%
%   Inputs:
%    - in: Either a vector of signals [n x 1] or a vector of AR coefficients [m x 1],
%       where n is the number of samples and m is the model order, respectively
%    - fs: Sampling frequency in Hz
%    - config: Optional struct that contains additional parameters
%       inputType: String defining which type of input is given, either the signal itself
%       'signal', or the AR coefficients 'coeff'. Default is 'signal'
%
%   Outputs:
%    - Figure containing the logarithmic PSD of the input
%    - pxx: Optional output, containing the values of the PSD
%
%   See also: dtf, pwelch, estimate_ar_coefficients
%

bool_isSignal=false;
bool_isCoeff=false;

freqRange=1:(fs/2);
nFreqs=length(freqRange);

x=[];
ar=[];

if isstruct(config)
    if isfield(config,'inputType')
        if strcmp(config.inputType,'signal')
            bool_isSignal=true;
            x=in;
        elseif strcmp(config.inputType,'coeff')
            bool_isCoeff=true;
            ar=in;
        end
    end
end

if isempty(x) && isempty(ar)
    x=in;
    bool_isSignal=true;
end

if bool_isSignal
    pxx=pwelch(x,fs,fs/2,1:(fs/2),fs); 
    figure; plot(10*log10(pxx));
elseif bool_isCoeff
    pxx=nan(nFreqs,1);
    
    for i=1:nFreqs
        
    end
end

end