
function test_suite = test_getVarSpec
  initTestSuite;

% tests function omSq=getVarSpec(X)
  
function test_sinusoid

T = 100;

om = 2*pi*[1/5,1/10];

X = cos([0:T]'*om);

freqSq1 = (om/(2*pi)).^2;

freqSq2=getVarSpec(X);

tol = 1e-3;
assertVectorsAlmostEqual(freqSq1,freqSq2,'absolute',tol,0)


function test_whiteNoise

T = 10000;

X = randn(T,1);

freqSq1 = mean(linspace(-1/2,1/2,T)'.^2);

freqSq2=getVarSpec(X);

[freqSq1,freqSq2]

tol = 1e-1;
assertVectorsAlmostEqual(freqSq1,freqSq2,'absolute',tol,0)


