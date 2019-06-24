function phi_hat=calculate_serial_coefficients(E,order)
%% phi_hat=calculate_serial_coefficients(E,order)
%
%   Inputs:
%    - E: Residuals from the AR model. Matrix with size [n x c], where n is the number of
%       samples and c is the number of channels
%    - order: Model order (number of lags) to calculate the coefficients for
%
%   Outputs:
%    - phi_hat: Serial correlated coefficients of the residuals. Matrix of size [o x c],
%       where o is the model order, and c is the number of channels
%

numChannels=size(E,2);

phi_hat=zeros(order,numChannels);

for i=1:order
    phi_hat(i,:) = diag((E((1+i):end,:)'*E(1:(end-i),:)) ./ (E((1+i):end,:)'*E((1+i):end,:)));
end

end