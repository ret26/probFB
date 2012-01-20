function fftCov = getGPSESpec(len,T);
  

% function fftCov = getGPSESpec(len,T);
%
% Returns the spectrum of a GP with a squared exponential kernel.
%
% <x(t)x(t')> = exp(-1/(2*len^2)*(t-t')^2)
%
% INPUTS
% len = length scale of the GP
% T = number of frequencies over which to evaluate the spectrum
%
% OUTPUTS
% fftCov = power spectrum of the GP, size [T,1]

    
  omega = 2*pi*[[0:floor(T/2)],[-floor(T/2)+1:1:-1]]'/T;

  fftCov = sqrt(2*pi*len^2)*exp(-omega.^2*len.^2/2);
  
%  tiny = 1e-6;
%  fftCov = fftCov+tiny;
  