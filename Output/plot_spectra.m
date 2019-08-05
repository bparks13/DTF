function plot_spectra(S,spectral_range,S_orig)
%% plot_spectra(S,spectral_range,S_orig)
%
%  Given the spectra calculated from the AR coefficients, plot the spectra of the AR
%  coefficients. If the original coefficients are known (simulated data), both can be
%  plotted on the same plot to compare
%
%   Inputs:
%    - S: Spectra calculated from the estimated AR coefficients. Size is [c x c x f],
%       where c is the number of channels and f is the number of frequencies that the
%       spectra is calculated over
%    - spectral_range: Vector defining the range of frequencies over which the spectrum is
%       calculated, for plotting purposes
%    - S_orig: Optional input. If the original AR coefficients are known, the original
%       spectra can also be plotted to compare the two spectrums
%
%   Outputs:
%    Figure containing subplots for each channel and its corresponding spectra
%

bool_SingleSpectra=true;
bool_reformatOrig=false;

if nargin==3
    bool_SingleSpectra=false;
    
    if size(S_orig,3)==1 % If the FFT output is given instead of the spectral output, reformat when plotting
        bool_reformatOrig=true;
    end
end

numChannels=size(S,1);

figure;

for i=1:numChannels
    subplot(numChannels,1,i);
    plot(spectral_range,abs(squeeze(S(i,i,:))));
    
    if ~bool_SingleSpectra
        hold on;
        
        if bool_reformatOrig
            plot(spectral_range,S_orig(:,i));
        else
            plot(spectral_range,abs(squeeze(S_orig(i,i,:))));
        end
        legend('S','S\_orig')
    else
        legend('S')
    end
end

end