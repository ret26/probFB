function test_suite =  test_packDimsMPAD
  initTestSuite;

% test 
% 

function test_correctAutoCorrelationGPs

D = rand;
K = rand;
T = rand;

Dims = packDimsMPAD(D,K,T);

[D2,K2,T2] = unpackDimsMPAD(Dims);

tol = 1e-5;
assertVectorsAlmostEqual(D2,D,'absolute',tol,0)
assertVectorsAlmostEqual(K2,K,'absolute',tol,0)
assertVectorsAlmostEqual(T2,T,'absolute',tol,0)
