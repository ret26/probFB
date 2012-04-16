function test_suite = test_MAAutoCor2Phi
  initTestSuite;

% Tests:
%
% function function autoCor = getAutoCorMA(phi)
% 
  
function testAR2

fmax = 100/1000;
df = 50/1000;
varMa = 1;
T = 100;

% AR2 version of the model
[lamx,varx] = freq2AR2(fmax,df,varMa);

% Autocorrelation of the model
autoCor = getAutoCorARTau(lamx',varx,T);

% Moving average version of the model
numIts = 1000;
smth = 0;

phi = MAAutoCor2Phi(autoCor,numIts,smth);

autoCorMA = getAutoCorMA(phi);

tol = 1e-1;
assertElementsAlmostEqual(autoCorMA,autoCor,'absolute',tol,0)

