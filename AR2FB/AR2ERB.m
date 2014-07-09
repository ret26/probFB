function [WERB,cf,H] = AR2ERB(lam,varx)

% function WERB = AR2ERB(lam,varx)
%
% computes the Equivalent Rectangular Bandwidth for AR2 filters
%
% INPUTS
% lams = AR(2) dynamical parameter [2,D]
% varx = AR(2) dynamical noise [1,D]
%
% OUTPUTS
% ERB = equivalent rectangular bandwidth (relative to sample rate
%       i.e. the max value of ERB is 1/2) [1,D]
% cf = centre frequencies [1,D]
% H = peak power at centre-frequency [1,D]

D = length(varx);
NumFreqs = 10000;
RngFreqs = [0,1/2];
WERB = zeros(1,D);
H = zeros(1,D);
cf = zeros(1,D);

for d=1:D
  [Freqs,SpecTau] = getSpecAR2(lam(:,d)',varx(d),NumFreqs,RngFreqs);
  dFreq = Freqs(2)-Freqs(1);
  I = sum(SpecTau)*dFreq;
  [H(d),pos] = max(SpecTau);
  cf(d) = Freqs(pos);
  WERB(d) = I/H(d); 
end


