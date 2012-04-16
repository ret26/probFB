
function test_suite = test_getmVarAR2
  initTestSuite;

% function mVar = getmVarAR2(Lam,Var)
  
function test_AR1

dispFigs = 0;
D= 5;
Lam = [rand(D,1),zeros(D,1)];
Var = exp(randn(D,1));

mVar1 = getmVarAR2(Lam,Var);
mVar2 = Var./(1-Lam(:,1).^2);


tol = 1e-3;
assertVectorsAlmostEqual(mVar1,mVar2,'absolute',tol,0)


function test_AR2_specialCase

dispFigs = 0;
Lam = [1,-1/2];
Var = exp(randn);
mVar1 = getmVarAR2(Lam,Var);
mVar2 = 2.4*Var;


tol = 1e-3;
assertVectorsAlmostEqual(mVar1,mVar2,'absolute',tol,0)


