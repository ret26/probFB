
function test_suite = test_sampleARTau
  initTestSuite;

function test_compareSampleWithAutoCorrelation

dispFigs = 0;
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

function test_compareFistTauSamplesWithAutoCorrelation

dispFigs = 0;
lam = [1,-0.4];
varx = 1;

T = 2; N = 50000; 
X = sampleARTau(lam,varx,T,N);
tau = length(lam);

tau = length(lam);
TAuto = tau;

autoCorEst(1,1)= mean(X(1,:).^2);
autoCorEst(2,1)= mean(X(1,:).*X(2,:));

%tic
autoCor = getAutoCorARTau(lam',varx,TAuto);
%toc

if dispFigs==1
  figure
  hold on
  plot(autoCor,autoCorEst,'.r')
  plot([-1,1],[-1,1],'-k')
end

tol = 1e-1;
assertVectorsAlmostEqual(autoCor,autoCorEst,'absolute',tol,0)

