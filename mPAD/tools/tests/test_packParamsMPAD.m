function test_suite =  test_packParamsMPAD
  initTestSuite;

% test 
% 

function test_correctAutoCorrelationGPs

G = rand;
Lam1 = rand;
Var1 = rand;

Len2 = rand;
Var2 = rand;

Mu2 = rand;
vary = rand;

Params = packParamsMPAD(G,Lam1,Var1,Len2,Var2,Mu2,vary);

[G2,Lam12,Var12,Len22,Var22,Mu22,vary2] = unpackParamsMPAD(Params);

tol = 1e-5;
assertVectorsAlmostEqual(G,G2,'absolute',tol,0)
assertVectorsAlmostEqual(Lam1,Lam12,'absolute',tol,0)
assertVectorsAlmostEqual(Var1,Var12,'absolute',tol,0)
assertVectorsAlmostEqual(Len2,Len22,'absolute',tol,0)
assertVectorsAlmostEqual(Var2,Var22,'absolute',tol,0)
assertVectorsAlmostEqual(Mu2,Mu22,'absolute',tol,0)
assertVectorsAlmostEqual(vary,vary2,'absolute',tol,0)
