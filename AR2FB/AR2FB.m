function [X,covX] = AR2FB(y,Lam,Var,vary,varargin)
  
  % function [X,covX] = AR2FB(y,Lam,Var,vary)
  %
  % Second order auto-regressive filter bank
  %
  % With optional inputs and outputs:
  % function [X,covX] = AR2FB(y,Lam,Var,vary,verbose,KF)
  %
  % INPUTS
  % y = data [T,1]
  % Lam = dynamical parameters [D,2]
  % Var = dynamic noise variance [D,1]
  % vary = observation noise 
  % OPTIONAL INPUTS:
  % verbose = 1 => verbose output
  % KF = 1 => Kalman Filter (rather than smoothing)
  %  
  % OUTPUTS
  % X = AR2 process mean values, size [D,T]
  % covX = Covariances of these values, [D,D,T]
  %
  % I could modify this to return the sufficient statistics and
  % likelihood if requested by the user - see kalman.m
  
T = length(y);

if nargin<=4
  verbose = 0 ;
else
  verbose = varargin{1};
end

if nargin<=5
  KF = 0 ;
else
  KF = varargin{2};
end

% Get parameters for smoothing
[A,Q,C,R,x0,P0] =  ar2LDSParams(Lam,Var,vary);
  
% Kalman Smoothing
[lik,Xfin,Pfin] = kalman(A,C,Q,R,x0,P0,reshape(y,[1,1,T]),verbose,KF);

% Output
[X,covX] = getAR2LDSOutput(Xfin,Pfin);

