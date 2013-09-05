function test_suite = test_get_obj_recon_FB_mag
  initTestSuite;

% Tests:
%
% function [obj,dobj] = get_obj_recon_FB_mag(y,aTar,spec)
% 
  
function testRunWithoutError

% Tests the function runs at all - useful for debugging

T = 100;
D = 3;
spec = randn(T,D);
aTar = exp(randn(T,D));
y = randn(T,1);

[obj,dobj] = get_obj_recon_FB_mag(y,aTar,spec,'a');

function testGradientNumerically_a

% numerically evaluates the gradient and compares the output of the
% code to this
T = 100;
D = 3;
spec = randn(T,D);
aTar = exp(randn(T,D));
y = randn(T,1);

[obj,dobj] = get_obj_recon_FB_mag(y,aTar,spec,'a');

delta =1e-6;
vals = checkgrad('get_obj_recon_FB_mag',y,delta,aTar,spec,'a');

tol = 1e-6;
assertElementsAlmostEqual(vals,0,'absolute',tol,0)

function testGradientNumerically_loga

% numerically evaluates the gradient and compares the output of the
% code to this
T = 100;
D = 3;
spec = randn(T,D);
aTar = exp(randn(T,D));
y = randn(T,1);

[obj,dobj] = get_obj_recon_FB_mag(y,aTar,spec,'loga');

delta =1e-6;
vals = checkgrad('get_obj_recon_FB_mag',y,delta,aTar,spec,'loga');

tol = 1e-6;
assertElementsAlmostEqual(vals,0,'absolute',tol,0)


function testGradientNumerically_prior

% numerically evaluates the gradient and compares the output of the
% code to this
T = 100;
D = 3;
spec = randn(T,D);
aTar = exp(randn(T,D));
y = randn(T,1);
lam = 0.9;
vary = (1-lam^2)*1;

[obj,dobj] = get_obj_recon_FB_mag(y,aTar,spec,'loga',lam,vary);

delta =1e-6;
vals = checkgrad('get_obj_recon_FB_mag',y,delta,aTar,spec,'loga',lam,vary);

tol = 1e-6;
assertElementsAlmostEqual(vals,0,'absolute',tol,0)
