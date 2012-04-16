
function test_suite = test_cosFreq2AR2
  initTestSuite;

% function [Lam,Var] = cosFreq2AR2(cosCF,cosDF,mVar);
  
function test_consistent

rand('state',6);
D = 5;
CF = rand(D,1)/2;
maxDF = min([1/2-CF,CF]')'*2;

mVar= exp(randn(D,1));
DF = maxDF.*rand(D,1);

[cosCF,cosDF] = CFDF2cosCFDF(CF,DF);

[Lam,Var] = cosFreq2AR2(cosCF,cosDF,mVar);

[CF2,DF2,mVar2] = AR22freq(Lam,Var);

tol = 1e-5;
assertVectorsAlmostEqual(CF,CF2,'absolute',tol,0)
assertVectorsAlmostEqual(DF,DF2,'absolute',tol,0)
assertVectorsAlmostEqual(mVar,mVar2,'absolute',tol,0)


