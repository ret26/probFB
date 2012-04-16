
function test_suite = test_getSpecAR2cosFreq
  initTestSuite;

  
function test_compareWith_getSpecAR2

dispFigs = 0;
lam = [1,-0.5];
varx = 1;
tau = length(lam);

numAves = 1000;
NumFreqs = 50;
RngFreqs = [0,1/2];
[Freqs,Spec1] = getSpecAR2(lam,varx,NumFreqs,RngFreqs);


[CF,DF, mVar] = AR22freq(lam,varx);
[cosCF,cosDF] = CFDF2cosCFDF(CF,DF);
 
Spec2 = getSpecAR2cosFreq(cosCF,cosDF,mVar,Freqs);


if dispFigs==1
  figure
  hold on
  plot(Freqs,Spec1,'.r')
  plot(Freqs,Spec2,'-k')
end

tol = 1e-3;
assertVectorsAlmostEqual(Spec1,Spec2,'absolute',tol,0)

