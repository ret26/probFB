function test_suite = test_getSpecAR2Obj
  initTestSuite;

% Tests:
%
% function [Obj,dObj]=getSpecAR2Obj(theta,vary,specTar);
% 
  
function testRunWithoutError

% Tests the function runs at all - useful for debugging

N = 100;
D= 5;
theta = [randn(D,1);log(rand(2*D,1)/2)];
vary = rand;
specTar = abs(fft(randn(N,1))).^2;
limCF = cumsum(rand(D,2)'/4)';
limDF = cumsum(rand(D,2)'/4)';
minVar = exp(randn(D,1));
bet = rand;

[Obj,dObj]=getSpecAR2Obj(theta,vary,bet,specTar,minVar,limCF,limDF);
Obj=getSpecAR2Obj(theta,vary,bet,specTar,minVar,limCF,limDF);

function testGradientNumerically

% numerically evaluates the gradient and compares the output of the
% code to this

N = 100;
D= 5;
theta = [randn(D,1);log(rand(2*D,1)/2)];
vary = rand;
specTar = abs(fft(randn(N,1))).^2;

limCF = cumsum(rand(D,2)'/4)';
limDF = cumsum(rand(D,2)'/4)';
minVar = exp(randn(D,1));
bet = rand;

delta =1e-6;
vals = checkgrad('getSpecAR2Obj',theta,delta,vary,bet,specTar,minVar,limCF,limDF);

tol = 1e-6;
assertElementsAlmostEqual(vals,0,'absolute',tol,0)
