function test_suite =  test_getAR2LDSOutput
  initTestSuite;

% test function [X,covX] = getAR2LDSOutput(Xfin,Pfin);

function test_interfaceWithKalman

Lam = [1.5,-0.7;
       1,-0.8];

Var = [1,1.2];

vary = 0.03;

[A,Q,C,R,x0,P0] =  ar2LDSParams(Lam,Var,vary);

T = 10;
Y = randn(1,1,T);

[lik,Xfin,Pfin,varargout] = kalman(A,C,Q,R,x0,P0,Y);
 
[X,covX] = getAR2LDSOutput(Xfin,Pfin);
 
tol = 1e-5;
assertVectorsAlmostEqual(X,squeeze(Xfin(1,[1,3],:)),'absolute',tol,0)
assertVectorsAlmostEqual(covX(1,1,:),Pfin(1,1,:),'absolute',tol,0)
assertVectorsAlmostEqual(covX(2,2,:),Pfin(3,3,:),'absolute',tol,0)
assertVectorsAlmostEqual(covX(2,1,:),Pfin(1,3,:),'absolute',tol,0)
