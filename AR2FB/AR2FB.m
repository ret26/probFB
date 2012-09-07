function [X,varargout] = AR2FB(y,Lam,Var,vary,varargin)
  
  % function X = AR2FB(y,Lam,Var,vary)
  %
  % Second order auto-regressive filter bank (see Turner 2010,
  % Chapter 5 for details)
  %
  % In the standard mode above the FFT is used to compute the AR
  % filter bank. This is fast. When the optional inputs/output below
  % are added/requested, the kalman mode is used which is slower.
  % 
  % NOTE: THAT THE FFT METHOD GIVES A SLIGHTLY DIFFERENT SOLUTION
  % FROM THE KALMAN BASED METHODS (SEE BELOW). SPECIFICALLY, AT THE
  % START AND END OF THE SIGNAL THERE ARE DISCREPANCIES DUE TO THE
  % CIRCULAR BOUNDARY CONDITIONS ASSUMED. IN GENERAL THOUGH, THE
  % TWO METHODS WILL BE EXTREMELY SIMILAR.
  %
  % With optional inputs and outputs:
  % function [X,covX] = AR2FB(y,Lam,Var,vary,verbose,KF)
  %
  % In this mode (or when the noise is non-stationary), the Kalman
  % smoother/filter is used to compute the AR filter bank. This is
  % rather slower, but it can handle non-stationary noise and it
  % returns the uncertainty values. 
  %
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
  % OPTIONAL OUTPUTS:
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

if KF==0 & nargout==1 & length(vary)==1

  X = AR2FBFFT(y,Lam,Var,vary);
  
else
  % Get parameters for smoothing
  [A,Q,C,R,x0,P0] =  ar2LDSParams(Lam,Var,vary);
  
  % Kalman Smoothing
  [lik,Xfin,Pfin] = kalman(A,C,Q,R,x0,P0,reshape(y,[1,1,T]),verbose,KF);

  % Output
  [X,covX] = getAR2LDSOutput(Xfin,Pfin);
end

if nargout==2
  varargout(1) = {covX};
else
  varargout(1) = {[]};
end