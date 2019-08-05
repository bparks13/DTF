% pr14_1.m 
% A test for identifying coefficients from a time series

% The program prints the coefficients we use (a and b)
% plus their estimates (vector A) obtained with the arfit function
% ['If needed INSTALL ARfit toolbox: http://climate-dynamics.org/software/#arfit']
% ['put the toolbox files in your current folder, or use "set path" in the ']
% ['File menu to include it in the MALAB paths']
% installed in path C:\Program Files\MATLAB71\toolbox\ar_tools
CCC;
% Set coefficients a and b
a=0.95;         
b=-.55;
% Parameters & initial values for the time series  
N=1000;        
e=randn(1,N+3);                         % GWN input
x(1)=0;e(1)=0;
x(2)=0;e(2)=0;

% create the time series using autoregression
for i=3:N+3
    x(i)=a*x(i-1)+b*x(i-2)+e(i);
end;
% Remove the 1st 2nd zero-valued-points from x
x=x(3:N+3); 
% Make e the same length
e=e(3:N+3);
% normalize x & e
x=x-mean(x);
x=x/std(x)^2;
e=e-mean(e);
e=e/std(e)^2;

% ARfit toolbox 
[w,A,C,SBC,FPE,th]=arfit(x',0,4);       % arfit function
% Show the results
('coefficients in aurtoregression time series')
[a b]
('estimated coefficients using arfit')
A