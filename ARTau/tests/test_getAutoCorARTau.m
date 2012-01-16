
function test_suite = test_getAutoCorARTau
  initTestSuite;

function test_compareWithSample

dispFigs = 1;
lam = [1,-0.5];
varx = 1;

T = 100000; N = 1; 
X = sampleARTau(lam,varx,T,N);
tau = length(lam);

TAuto = tau+10;

tau = length(lam);
AutoCorEst= zeros(TAuto,1);

for t=1:TAuto
  AutoCorEst(t) = mean(X(1+t-1:T).*X(1:T-t+1));
end

%tic
AutoCor = getAutoCorARTau(lam',varx,TAuto);
%toc

if dispFigs==1
  figure
  hold on
  plot(AutoCor,AutoCorEst,'.r')
  plot([-1,1],[-1,1],'-k')
end

tol = 1e-1;
assertVectorsAlmostEqual(AutoCor,AutoCorEst,'absolute',tol,0)

