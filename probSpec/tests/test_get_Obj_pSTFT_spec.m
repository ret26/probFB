function test_suite = test_get_Obj_pSTFT_spec
  initTestSuite;

% Tests:
%
% function [Obj,dObj]=getSpecAR2Obj(theta,vary,specTar);
% 
  
function testRunWithoutError

% Tests the function runs at all - useful for debugging

N = 100;
D= 5;
theta = randn(3*D,1);
vary = rand;
specTar = abs(fft(randn(N,1))).^2;
limOm = cumsum(rand(D,2)'/4)';
limLam = cumsum(rand(D,2)'/2)';
minVar = exp(randn(D,1));
bet = rand;


[Obj]= get_Obj_pSTFT_spec(theta,vary,specTar,minVar,limOm, ...
				     limLam,bet);

[Obj,dObj]= get_Obj_pSTFT_spec(theta,vary,specTar,minVar,limOm, ...
				     limLam,bet);

function testGradientNumerically

% numerically evaluates the gradient and compares the output of the
% code to this

N = 100;
D= 5;
theta = randn(3*D,1);
vary = rand;
specTar = abs(fft(randn(N,1))).^2;
limOm = cumsum(rand(D,2)'/4)';
limLam = cumsum(rand(D,2)'/2)';
minVar = exp(randn(D,1));
bet = rand;

delta =1e-6;
vals = checkgrad('get_Obj_pSTFT_spec',theta,delta,vary,specTar,minVar,limOm, ...
				     limLam,bet);

 tol = 1e-6;
 assertElementsAlmostEqual(vals,0,'absolute',tol,0)
