function test_suite = test_sampleAR2FB
  initTestSuite;

% Tests: 
%
% function [Y,X] = sampleAR2FBSlow(Lam,Var,vary,T);
%
function testStatisticsMatch

% Tests the auto correlation of the samples matches the
% analytic auto-correlation of the AR(2) processes

dispFigs=0;

randn('seed',2)

 T = 100000;
 D = 2;
 K = 2;

 Lam = [1.1,-0.9;
        1.5,-0.95];
 Var = [1,0.5]';

 vary = 0;
 
autoCor1A = getAutoCorARTau(Lam(1,:)',Var(1),2);
autoCor2A = getAutoCorARTau(Lam(2,:)',Var(2),2);

[Y,X] = sampleAR2FB(Lam,Var,vary,T);

autoCor1B = [var(X(:,1)); mean(X(2:T,1).*X(1:T-1,1))];
autoCor2B = [var(X(:,2)); mean(X(2:T,2).*X(1:T-1,2))];
 
tol =3e-1;
assertVectorsAlmostEqual(autoCor1A,autoCor1B,'absolute',tol,0)
assertVectorsAlmostEqual(autoCor2A,autoCor2B,'absolute',tol,0)

tol =1e-5;
assertVectorsAlmostEqual(sum(X,2),Y,'absolute',tol,0)

function testFirstAndLastSamplesIndependent

% Makes sure that there are no
% wrap around effects arising
% from the use of the fft

dispFigs=1;

randn('seed',2)

 T = 1000;
 D = 2;
 K=1000;

 Lam = [1.1,-0.9;
        1.5,-0.95];
 Var = [1,0.5]';

 vary = 0;

% To test whether the first and the last sample are correlated:

Xs = zeros(T,D,K);
for k=1:K
  [Y,X] = sampleAR2FB(Lam,Var,vary,T);
  Xs(:,:,k) =  X;
end


autoCorA = getAutoCorARTau(Lam(1,:)',Var(1),T);
autoCorB = getAutoCorARTau(Lam(2,:)',Var(2),T);
autoCorEmpA = zeros(T,1);
autoCorEmpB = zeros(T,1);
for t=1:T
  autoCorEmpA(t) = mean(Xs(1,1,:).*Xs(t,1,:));
  autoCorEmpB(t) = mean(Xs(1,2,:).*Xs(t,2,:));
end



if dispFigs==1
  figure
  subplot(2,1,1)
  title('should be independent - i.e. no correlation')
  plot(squeeze(Xs(1,1,:)),squeeze(Xs(T,1,:)),'.k')
  subplot(2,1,2)
  plot(squeeze(Xs(1,2,:)),squeeze(Xs(T,2,:)),'.k')

  figure
  subplot(2,1,1)
  plot(autoCorEmpA,'-k','linewidth',2)
%  plot(autoCorA,'-r')

  subplot(2,1,2)
  plot(autoCorEmpB,'-k','linewidth',2)
%  plot(autoCorB,'-r')

end

%keyboard
%tol =3e-1;
%assertVectorsAlmostEqual(xCorLastFirstA,0,'absolute',tol,0)
%assertVectorsAlmostEqual(xCorLastFirstB,0,'absolute',tol,0)
