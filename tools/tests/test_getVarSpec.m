
function test_suite = test_getVarSpec
  initTestSuite;

% tests function omSq=getVarSpec(X)
  
function test_sinusoid

T = 100;

om = 2*pi*[1/10,1/50];

X = cos([0:T]'*om+2*pi*rand);

freqSq1 = (om/(2*pi)).^2;

freqSq2=getVarSpec(X);

%[freqSq1;freqSq2]

tol = 1e-3;
assertVectorsAlmostEqual(freqSq1,freqSq2,'absolute',tol,0)


function test_whiteNoise

T = 10000;

X = randn(T,1);

freqSq1 = mean(linspace(-1/2,1/2,T)'.^2);

freqSq2=getVarSpec(X);

%[freqSq1,freqSq2]

tol = 1e-1;
assertVectorsAlmostEqual(freqSq1,freqSq2,'absolute',tol,0)


function test_multiSinusoids

T = 2000;

om1 = 2*pi*[1/10,1/5];
om2 = 2*pi*[1/7,1/3];

alp = 0.3;

X = alp*cos([0:T]'*om1+2*pi*rand)+(1-alp)*cos([0:T]'*om2+2*pi*rand);


freqSq1 = (alp^2*(om1/(2*pi)).^2+(1-alp)^2*(om2/(2*pi)).^2)/(alp^2+(1-alp)^2);

freqSq2=getVarSpec(X);

%[freqSq1;freqSq2]

tol = 1e-3;
assertVectorsAlmostEqual(freqSq1,freqSq2,'absolute',tol,0)
