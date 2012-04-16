function test_suite = test_getMAAutoCorObj
  initTestSuite;

% Tests:
%
% function [Obj,dObj] = getMAAutoCorObj(phi,autoCorTar,smth)
% 
  
function testRunWithoutError

% Tests the function runs at all - useful for debugging

T = 5;
autoCorTar = randn(T,1);
phi = randn(T,1);
smth = 1;
[Obj,dObj] = getMAAutoCorObj(phi,autoCorTar,smth);



function testTequal2

% numerically evaluates the gradient and compares the output of the
% code to this

T = 2;
autoCorTar = randn(T,1);
phi = randn(T,1);
delta= 1e-5;

[Obj,dObj] = getMAAutoCorObj(phi,autoCorTar,0);

delta1 = phi(1)^2+phi(2)^2-autoCorTar(1);
delta2 = phi(1)*phi(2)-autoCorTar(2);
Obj2 = 1/2*delta1^2+1/2*delta2^2;

dObj2(1,1) = 2*delta1*phi(1)+delta2*phi(2);
dObj2(2,1) = 2*delta1*phi(2)+delta2*phi(1);

tol = 1e-6;
assertElementsAlmostEqual(Obj2,Obj,'absolute',tol,0)
assertElementsAlmostEqual(dObj2,dObj,'absolute',tol,0)

function testGradientNumerically

% numerically evaluates the gradient and compares the output of the
% code to this

T = 16;
autoCorTar = randn(T,1);
phi = randn(T,1);
delta= 1e-5;
smth = 1;

vals = checkgrad('getMAAutoCorObj',phi,delta,autoCorTar,smth);

tol = 1e-6;
assertElementsAlmostEqual(vals,0,'absolute',tol,0);
