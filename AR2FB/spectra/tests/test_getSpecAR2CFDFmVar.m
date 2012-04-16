
function test_suite = test_getSpecAR2CFDFmVar
  initTestSuite;

  
function test_compareWithAlternativeSpectrumComputation

dispFigs = 0;

lam = [1,-0.5];
varx = 1;
tau = length(lam);

NumFreqs = 50;
RngFreqs = [0,1/2];
[Freqs1,Spec1] = getSpecAR2(lam,varx,NumFreqs,RngFreqs);

cosCF = lam(1)*(lam(2)-1)/(4*lam(2));
cosDF = sqrt(-(1+lam*lam')/lam(2)-4*cosCF^2-2);

mVar = getmVarAR2(lam,varx);

[Freqs2,Spec2] = getSpecAR2CFDFmVar(cosCF,cosDF,mVar,NumFreqs,RngFreqs);

if dispFigs==1
  figure
  hold on
  plot(Freqs1,Spec1,'.r')
  plot(Freqs2,Spec2,'-k')
end

tol = 1e-5;
assertVectorsAlmostEqual(Spec1,Spec2,'absolute',tol,0)

