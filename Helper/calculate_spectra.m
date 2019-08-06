function S=calculate_spectra(H,C)
%% S=calculate_spectra(H)
%
%  Helper function to take a transfer function and the covariance matrix, and produce the
%  spectral matrix, which is defined as S = H * C * H'
%
%   Inputs:
%    - H: Transfer function matrix. Size should be [c x c x f], where c is the number of
%       channels and f is the number of frequencies 
%    - C: Covariance matrix. Size is [c x c], where c is the number of channels
%
%   Outputs:
%    - S: Spectral matrix. Size is [c x c x f]

S=nan(size(H));

for i=1:size(H,3)
    S(:,:,i)=H(:,:,i) * C * H(:,:,i)';
end

end