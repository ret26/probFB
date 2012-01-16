function [X,covX] = getAR2LDSOutput(Xfin,Pfin);

% function [X,covX] = getAR2LDSOutput(Xfin,Pfin);
%
% Get the output in the correct format from the AR2 filter bank Kalman
% Smoother
%
% INPUTS
% Xfin = mean of the latent variables, size [N,2D,T] (N=1)
% Pfin = covariance of the latent variables, size []
%
% OUTPUTS
% X = mean of the latent variables, size [D,T]
% covX = covariance of the latent variables, size [D,D,T]

[N,TwoD,T] = size(Xfin);  

D =  TwoD/2;

X = Xfin(1,[1:2:TwoD],:);

covX = zeros(D,D,T);

for d=1:D
%  X(d,:) = reshape(Xfin(1,2*(d-1)+1,:),[1,T]);
  covX(d,:,:) = Pfin(2*(d-1)+1,[1:2:TwoD],:);
end

