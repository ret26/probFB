function test_suite = test_freq2probSpec
  initTestSuite;

% Tests: 
%
% function [om,lamx,varx] = freq2probSpec(fmax,df,varMa)
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


