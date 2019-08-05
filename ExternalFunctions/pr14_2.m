%% pr14_2

CCC;

%% 
sr=1000;                    % sample rate is 1000 samples/s
dt=1/sr;                    % sample interval
Nyq=sr/2;                   % Nyquist frequency
T=10;                       % Epoch of 10s
F=0:1/T:Nyq;                % Frequency Scale
% the A matrices for the autoregressive process
A1=[.5 0 0;0 .8 0;0 0 .7];  % x1, x2, and x3 depend on their own k-1 values
A2=[-.8 0 0;0 0 0;0 0 0];   % only x1 depends on its own k-2 value
A4=[-.2 0 0;.4 0 0;0 0 0];  % x1 and x2 both depend on x1 k-4 value
A7=[0 0 0;.8 0 0;.4 0 0];   % x2 and x3 both depend on x1 k-7 value
% you may uncomment the following line to add a dependence x2 -> x3
%A7=[0 0 0;.8 0 0;.4 .9 0]; % x2 and x3 both depend on x1 k-7 value and 
                            % x3 also on x2 k-7 value
A13=[0 0 0;0 0 0;.9 0 0];   % x3 depends on x1 k-13 value
I=eye(3);
%time series
X=zeros(3,length(F));
% loop to simulate the time series note it starts at 14 since the max delay
% is 13 [as given by A13]
for k=14:length(F)
    E=randn(3,1);
    X(:,k)=A1*X(:,k-1)+A2*X(:,k-2)+A4*X(:,k-4)+A7*X(:,k-7)+A13*X(:,k-13)+E;
end

%%
k=0;
for f=0:1/T:Nyq                                     % LOOP to compute H
    k=k+1;
    % compute the z values for each delay and use the imaginary part
    % j*2*pi*f*dt to obtain the values for the frequency response
    z1=exp(-1i*2*pi*f*dt);
    z2=z1^2;
    z4=z1^4;
    z7=z1^7;
    z13=z1^13;
    
    H(:,:,k)=inv(I-A1*z1-A2*z2-A4*z4-A7*z7-A13*z13);  % Eq (14.5c)
    H2(:,:,k)=abs(H(:,:,k)).^2;
end

%%
% Loops to plot the signals & compute and plot DTF
figure;hold on;
xmax=max(F);
% first plot the three signals
subplot(4,3,1);
plot(X(1,:),'k');hold;axis([0 length(F) min(X(1,:)) max(X(1,:))]);
t='signal-1';title(t);
axis('off');
subplot(4,3,2);
plot(X(2,:),'k');hold;axis([0 length(F) min(X(2,:)) max(X(2,:))]);
t='signal-2';title(t);
axis('off');
subplot(4,3,3);
plot(X(3,:),'k');hold;axis([0 length(F) min(X(3,:)) max(X(3,:))]);
t='signal-3';title(t);
axis('off');
% Next compute and plot the Spectra and DTFs
for n=1:3                               
    for m=1:3
        if n==m
            ttl=['Spectrum ' num2str(m)];
            % put the plot variable [the Spectrum] in a temporary array
            for k=1:5001;tmp(k)=H2(n,m,k);end
            ymax=max(1.*tmp);
        else
            ttl=[' ' num2str(m) '->' num2str(n)];
            DTF(n,m,:)=H2(n,m,:)./sum(H2(n,:,:));       % Eq (14.7)
            % put the plot variable [the DTF] in a temporary array
            for k=1:5001;tmp(k)=DTF(n,m,k);end
            ymax=1.;
        end
        subplot(4,3,n*3+m);
        plot(F,tmp,'k');
        axis([0 xmax 0 ymax]);
        title(ttl);
    end
end

return

%% 

AR=zeros(3,3,13); %#ok<*UNRCH>
AR(:,:,1)=A1;
AR(:,:,2)=A2;
AR(:,:,4)=A4;
AR(:,:,7)=A7;
AR(:,:,13)=A13;
