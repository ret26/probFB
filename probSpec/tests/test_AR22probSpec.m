function test_suite = test_AR22probSpec
  initTestSuite;

% Tests: 
%
% function [om,lamx,varx] = AR22probSpec(Lam,Var) 
%

function testComputeByHand

% Two time-step example which can be computed by hand

D = 2;

lamx = [.9;.7];
varx = rand(length(lamx),1);
om = [1/10,1/20]';

[Lam,Var] = probSpec2AR2(om,lamx,varx);
[om2,lamx2,varx2] = AR22probSpec(Lam,Var);

tol = 1e-5;
assertVectorsAlmostEqual(om2,om,'absolute',tol,0)
assertVectorsAlmostEqual(lamx2,lamx,'absolute',tol,0)
assertVectorsAlmostEqual(varx2,varx,'absolute',tol,0)
