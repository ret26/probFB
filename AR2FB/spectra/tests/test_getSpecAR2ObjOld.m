function test_suite = test_getSpecAR2ObjOld
  initTestSuite;

% Tests:
%
% function [Obj,dObj]=getSpecAR2ObjOld(theta,vary,specTar);
% 
  
function testRunWithoutError

% Tests the function runs at all - useful for debugging

N = 100;
D= 5;
theta = randn(3*D,1);
vary = rand;
specTar = abs(fft(randn(N,1))).^2;

[Obj,dObj]=getSpecAR2Obj(theta,vary,specTar);
Obj=getSpecAR2Obj(theta,vary,specTar);

function testGradientNumerically

% numerically evaluates the gradient and compares the output of the
% code to this

N = 100;
D= 5;
theta = randn(3*D,1);
vary = rand;
specTar = abs(fft(randn(N,1))).^2;

delta =1e-6;
vals = checkgrad('getSpecAR2Obj',theta,delta,vary,specTar);

tol = 1e-6;
assertElementsAlmostEqual(vals,0,'absolute',tol,0)
