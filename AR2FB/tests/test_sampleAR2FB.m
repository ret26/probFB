function test_suite = test_sampleAR2FB
  initTestSuite;

% Tests: 
%
% function [Y,X] = sampleAR2FBSlow(Lam,Var,vary,T);
%
function testStatisticsMatch

% Tests the auto correlation of the samples matches the
% analytic auto-correlation of the AR(2) processes

dispFigs=0;

randn('seed',2)

 T = 100000;
 D = 1;
 K = 2;

 Lam = [1.1,-0.9;
        1.5,-0.95];
 Var = [1,0.5]';

 vary = 0;
 
autoCor1A = getAutoCorARTau(Lam(1,:)',Var(1),2);
autoCor2A = getAutoCorARTau(Lam(2,:)',Var(2),2);

[Y,X] = sampleAR2FB(Lam,Var,vary,T);

autoCor1B = [var(X(:,1)); mean(X(2:T,1).*X(1:T-1,1))];
autoCor2B = [var(X(:,2)); mean(X(2:T,2).*X(1:T-1,2))];
 
tol =3e-1;
assertVectorsAlmostEqual(autoCor1A,autoCor1B,'absolute',tol,0)
assertVectorsAlmostEqual(autoCor2A,autoCor2B,'absolute',tol,0)

tol =1e-5;
assertVectorsAlmostEqual(sum(X,2),Y,'absolute',tol,0)

