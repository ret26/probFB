function [Z,covZ] = probFB(y,lamx,varx,om,vary,varargin)
  
  % function [Z,covZ] = probFB(y,lamx,varx,om,vary,varargin)
  %
  % Probabilistic Filter Bank (see Turner 2010, Chapter 5 for details)
  %
  % With optional inputs and outputs:
  % function [Z,covZ] = probSTFT(y,lamx,varx,om,vary,verbose,KF)
  %
  % INPUTS
  % lamx = dynamical AR parameters [D,1]
  % varx = dynamical noise parameters [D,1]
  % om = mean frequencies of the sinusoids [D,1]
  % vary = oberservation noise
  % y = Data, size [T,1] 
  %
  % OPTIONAL INPUTS:
  % verbose = 1 => verbose output
  % KF = 1 => Kalman Filter (rather than smoothing)
  %  
  % OUTPUTS
  % Z = probabilistic filter bank process mean values (these are
  %     complex), size [D,T]
  % covZ = Covariances of these values, [2D,2D,T]
  %
  % I could modify this to return the sufficient statistics and
  % likelihood if requested by the user - see kalman.m
  
T = length(y);

if nargin<=5
  verbose = 0 ;
else
  verbose = varargin{1};
end

if nargin<=6
  KF = 0 ;
else
  KF = varargin{2};
end
  
% Kalman Smoothing
[lik,Xfin,Pfin] = kalmanSlowFB(lamx,varx,om, vary,y,verbose,KF);

% Output
[Z,covZ] = getFBLDSOutput(Xfin,Pfin);

