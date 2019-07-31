% pr14_3
% cross-correlation
% run pr 14_2 first!!
msg=['run pr 14_2 first!! AND ZOOM IN to see the result']
c11=xcov(X(1,:),X(1,:),'coeff');
c12=xcov(X(1,:),X(2,:),'coeff');
c13=xcov(X(1,:),X(3,:),'coeff');
c23=xcov(X(2,:),X(3,:),'coeff');
figure;hold
plot(c11,'k.-');
plot(c12,'b');
plot(c13,'g');
plot(c23,'r');
title('Cross-Correlation: bl-autocor; blue-1,2; green-1,3; red-2,3')
ylabel('ZOOM IN!!')