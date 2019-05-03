%% create a "ground-truth" signal, feed it back to the varm/estimate functions
% This is to verify if the A(1,2) coefficient is 1?2 or 2?1

truth_model=ar(1).estMdl;

numperiods=length(x)-truth_model.P*2;
Y=forecast(truth_model,numperiods,x(:,:,1));