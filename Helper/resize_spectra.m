function P=resize_spectra(S)
%% P=resize_spectra(S)
%
%  Given a matrix that contains both auto- and cross-spectra components, return only the
%  auto spectral values. Additionally, P is defined as |S|
%
%   Inputs:
%    - S: Matrix containing auto- and cross-spectral components. Size is [c x c x f],
%       where c is the number of channels, and f is the number of frequencies
%    
%   Outputs:
%    - P: Matrix containing only the auto-spectral components. Size is [f x c]
%
%  See also: mvar, plot_connectivity
%

P=nan(size(S,3),size(S,1));

for j=1:size(S,1)
    P(:,j)=squeeze(abs(S(j,j,:)));
end

end