function test_suite = test_AR2ERB
  initTestSuite;

% Tests: 
%
% function [WERB,cf,H] = AR2ERB(lam,varx)
%
%

function testPlot

clear;

dispFigs = 1;
lam = [1,-0.9]';
varx = 1;
tau = length(lam);

[WERB,cf,H] = AR2ERB(lam,varx);

NumFreqs = 50000;
RngFreqs = [0,1/2];
[Freqs,SpecTau] = getSpecAR2(lam(:)',varx,NumFreqs,RngFreqs);

I1 = sum(SpecTau)*(Freqs(2)-Freqs(1));
I2 = WERB*H;

if dispFigs==1
  figure
  hold on
  plot(Freqs,SpecTau,'.r')
  plot([cf-WERB/2,cf-WERB/2,cf+WERB/2,cf+WERB/2],[0,H,H,0])
  [I1,I2,I1-I2]
end

tol = 1e-4;


assertVectorsAlmostEqual(I1,I2,'absolute',tol,0)
