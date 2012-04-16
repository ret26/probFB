function test_suite = test_getAutoCorMA
  initTestSuite;

% Tests:
%
% function autoCor = getAutoCorMA(phi)
% 
  
function testRunWithoutError

% Tests the function runs at all - useful for debugging

dispFig =0;
T = 5;
phi = randn(T,1);

autoCor = getAutoCorMA(phi);

K = 10000;
% noi = randn(K+T,1);
% x= zeros(K,1);
% for k=1:K
%   x(k) = phi'*noi(k+T-1:-1:k);
% end

x = sampleMA(phi,K);

for t=1:T
autoCorEmp(t) = mean(x(1:K-t+1).*x(t:K));
end

if dispFig==1
figure
hold on
plot(autoCor,'-k')
plot(autoCorEmp,'-r')
end

tol = 3e-1;
assertElementsAlmostEqual(autoCor,autoCorEmp','absolute',tol,0)

