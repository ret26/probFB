function test_suite = test_Z2_to_X2
  initTestSuite;

% Tests: 
%
% function [Obj,dObj] = getObj_mPAD_noG(x2,y,fftCov,Params);
%

function test_equivalent_GP_SE

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Forward model Settings
dispFigs=0;
T = 100000;
K = 2;
LAG = 100;
tol = 10;

Len2 = [5,3];
Var2 = [3,8];

Z2 = randn(T,K);
X2 = Z2_to_X2(Z2,Len2,Var2,tol);


% general parameters
cvarTheory(1,:) = Var2(1)*exp(-1/(2*Len2(1).^2)*(0:LAG-1).^2);
cvarTheory(2,:) = Var2(2)*exp(-1/(2*Len2(2).^2)*(0:LAG-1).^2);

for lag=1:LAG
  cvarExperi(1,lag) = sum(X2(1:end-lag+1,1).*X2(lag:1:end,1))/(T-lag+1);
  cvarExperi(2,lag) = sum(X2(1:end-lag+1,2).*X2(lag:1:end,2))/(T-lag+1);
end

if dispFigs==1
  
figure
subplot(2,1,1)
hold on
plot([0:LAG-1],cvarExperi(1,:),'-k','linewidth',2)
plot([0:LAG-1],cvarTheory(1,:),'-r')

subplot(2,1,2)
hold on
plot([0:LAG-1],cvarExperi(2,:),'-k','linewidth',2)
plot([0:LAG-1],cvarTheory(2,:),'-r')

end


delta1 = mean(abs(cvarTheory(1,:)-cvarExperi(1,:)));
delta2 = mean(abs(cvarTheory(2,:)-cvarExperi(2,:)));


tol = 1e-1;
assertVectorsAlmostEqual(delta1,0,'absolute',tol,0)
assertVectorsAlmostEqual(delta2,0,'absolute',tol,0)


