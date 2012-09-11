function test_suite = test_probSpec2freq
  initTestSuite;

% Tests: 
%
% function [fmax,df,varMa] = probSpec2freq(om,lamx,varx)
%

function testConsistent

dispFigs = 0;

fmax = [1/8,1/4]';
df = [1/32,1/80]';
mVar = [1,2]';

[om,lamx,varx] = freq2probSpec(fmax,df,mVar);
[fmax2,df2,mVar2] = probSpec2freq(om,lamx,varx);

tol = 1e-3;
assertVectorsAlmostEqual(fmax,fmax2,'absolute',tol,0)
assertVectorsAlmostEqual(df,df2,'absolute',tol,0)
assertVectorsAlmostEqual(mVar,mVar2,'absolute',tol,0)



function test_small_lambda

dispFigs = 1;


om = pi/10; lamx = 0.1; varx = 1;

[fmax2,df2,mVar2] = probSpec2freq(om,lamx,varx);

tol = 1e-5;
assertVectorsAlmostEqual(1/2,df2,'absolute',tol,0)
