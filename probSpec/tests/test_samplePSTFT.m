function test_suite = test_samplePSTFT
  initTestSuite;

% Tests: 
%
% function [Y,X] = samplePSTFT(Lam,Var,Om,vary,T);
%
function testRealMatchAR2

% Tests the auto correlation of the real/imaginary parts of the samples matches
% those computed using the AR method.

dispFigs=1;

randn('seed',2)

 T = 100000;
 D = 2;
 K = 2;
 Lam1  = [0.9,0.7];
 Var1  = (1-Lam1.^2);
 Om = [pi/10,pi/20]';
 
 Lam2 = [2*Lam1'.*cos(Om),-Lam1'.^2];
 Var2 = (1+Lam1.^2).*Var1;

 vary = 0;
 
[Ya,Za] = samplePSTFTslow(Lam1,Var1,Om,vary,T);

X1 = real(Za);
tau = 50;

autoCor1A = zeros(tau,1);
autoCor2A = zeros(tau,1);

for t=1:tau
  autoCor1A(t) = mean(X1(t:T,1).*X1(1:T-t+1,1));
  autoCor2A(t) = mean(X1(t:T,2).*X1(1:T-t+1,2));
end

X2 = imag(Za);

autoCor1B = zeros(tau,1);
autoCor2B = zeros(tau,1);

for t=1:tau
  autoCor1B(t) = mean(X2(t:T,1).*X2(1:T-t+1,1));
  autoCor2B(t) = mean(X2(t:T,2).*X2(1:T-t+1,2));
end

%%%%%%

[Yb,Zb] = samplePSTFT(Lam1,Var1,Om,vary,T);

X1 = real(Zb);

autoCor1C = zeros(tau,1);
autoCor2C = zeros(tau,1);

for t=1:tau
  autoCor1C(t) = mean(X1(t:T,1).*X1(1:T-t+1,1));
  autoCor2C(t) = mean(X1(t:T,2).*X1(1:T-t+1,2));
end


X2 = imag(Zb);

autoCor1D = zeros(tau,1);
autoCor2D = zeros(tau,1);

for t=1:tau
  autoCor1D(t) = mean(X2(t:T,1).*X2(1:T-t+1,1));
  autoCor2D(t) = mean(X2(t:T,2).*X2(1:T-t+1,2));
end


if dispFigs==1
  figure
  subplot(2,1,1)
  hold on
  plot(autoCor1A,'-k','linewidth',2)
  plot(autoCor1B,'-r')
  plot(autoCor1C,'-c','linewidth',1/2)
  plot(autoCor1D,'-y','linewidth',1/2)
 
  subplot(2,1,2)
  hold on
  plot(autoCor2A,'-k','linewidth',2)
  plot(autoCor2B,'-r')
  plot(autoCor2C,'-c','linewidth',1/2)
  plot(autoCor2D,'-y','linewidth',1/2)

end

tol =3e-1;
assertVectorsAlmostEqual(autoCor1A,autoCor1B,'absolute',tol,0)
assertVectorsAlmostEqual(autoCor1A,autoCor1C,'absolute',tol,0)
assertVectorsAlmostEqual(autoCor2A,autoCor2B,'absolute',tol,0)
assertVectorsAlmostEqual(autoCor2A,autoCor2C,'absolute',tol,0)

tol =1e-5;
assertVectorsAlmostEqual(sum(real(Zb.*exp(i*[1:T]'*Om')),2),Yb,'absolute',tol,0)

function testRealMatchAR2CloseOm

% Tests the auto correlation of the real/imaginary parts of the samples matches
% those computed using the AR method. Testing where the frequencies
% are low, where the spectra from the negative/positive regions overlap

dispFigs=1;

randn('seed',2)

 T = 100000;
 D = 2;
 K = 2;
 Lam1  = [0.9,0.8];
 Var1  = (1-Lam1.^2);
 Om = [pi/100,pi/50]';
 
 Lam2 = [2*Lam1'.*cos(Om),-Lam1'.^2];
 Var2 = (1+Lam1.^2).*Var1;

 vary = 0;
 
[Ya,Za] = samplePSTFTslow(Lam1,Var1,Om,vary,T);

X1 = real(Za);
tau = 50;

autoCor1A = zeros(tau,1);
autoCor2A = zeros(tau,1);

for t=1:tau
  autoCor1A(t) = mean(X1(t:T,1).*X1(1:T-t+1,1));
  autoCor2A(t) = mean(X1(t:T,2).*X1(1:T-t+1,2));
end

X2 = imag(Za);

autoCor1B = zeros(tau,1);
autoCor2B = zeros(tau,1);

for t=1:tau
  autoCor1B(t) = mean(X2(t:T,1).*X2(1:T-t+1,1));
  autoCor2B(t) = mean(X2(t:T,2).*X2(1:T-t+1,2));
end

%%%%%%

[Yb,Zb] = samplePSTFT(Lam1,Var1,Om,vary,T);

X1 = real(Zb);

autoCor1C = zeros(tau,1);
autoCor2C = zeros(tau,1);

for t=1:tau
  autoCor1C(t) = mean(X1(t:T,1).*X1(1:T-t+1,1));
  autoCor2C(t) = mean(X1(t:T,2).*X1(1:T-t+1,2));
end


X2 = imag(Zb);

autoCor1D = zeros(tau,1);
autoCor2D = zeros(tau,1);

for t=1:tau
  autoCor1D(t) = mean(X2(t:T,1).*X2(1:T-t+1,1));
  autoCor2D(t) = mean(X2(t:T,2).*X2(1:T-t+1,2));
end


if dispFigs==1
  figure
  subplot(2,1,1)
  hold on
  plot(autoCor1A,'-k','linewidth',2)
  plot(autoCor1B,'-r')
  plot(autoCor1C,'-c','linewidth',1/2)
  plot(autoCor1D,'-y','linewidth',1/2)
 
  subplot(2,1,2)
  hold on
  plot(autoCor2A,'-k','linewidth',2)
  plot(autoCor2B,'-r')
  plot(autoCor2C,'-c','linewidth',1/2)
  plot(autoCor2D,'-y','linewidth',1/2)

end

tol =3e-1;
assertVectorsAlmostEqual(autoCor1A,autoCor1B,'absolute',tol,0)
assertVectorsAlmostEqual(autoCor1A,autoCor1C,'absolute',tol,0)
assertVectorsAlmostEqual(autoCor2A,autoCor2B,'absolute',tol,0)
assertVectorsAlmostEqual(autoCor2A,autoCor2C,'absolute',tol,0)

tol =1e-5;
assertVectorsAlmostEqual(sum(real(Zb.*exp(i*[1:T]'*Om')),2),Yb,'absolute',tol,0)



% function testFirstAndLastSamplesIndependent

% % Makes sure that there are no
% % wrap around effects arising
% % from the use of the fft

% dispFigs=1;

% randn('seed',2)

%  T = 1000;
%  D = 2;
%  K=1000;

%  Lam = [1.1,-0.9;
%         1.5,-0.95];
%  Var = [1,0.5]';

%  vary = 0;

% % To test whether the first and the last sample are correlated:

% Xs = zeros(T,D,K);
% for k=1:K
%   [Y,X] = sampleAR2FB(Lam,Var,vary,T);
%   Xs(:,:,k) =  X;
% end


% autoCorA = getAutoCorARTau(Lam(1,:)',Var(1),T);
% autoCorB = getAutoCorARTau(Lam(2,:)',Var(2),T);
% autoCorEmpA = zeros(T,1);
% autoCorEmpB = zeros(T,1);
% for t=1:T
%   autoCorEmpA(t) = mean(Xs(1,1,:).*Xs(t,1,:));
%   autoCorEmpB(t) = mean(Xs(1,2,:).*Xs(t,2,:));
% end



% if dispFigs==1
%   figure
%   subplot(2,1,1)
%   title('should be independent - i.e. no correlation')
%   plot(squeeze(Xs(1,1,:)),squeeze(Xs(T,1,:)),'.k')
%   subplot(2,1,2)
%   plot(squeeze(Xs(1,2,:)),squeeze(Xs(T,2,:)),'.k')

%   figure
%   subplot(2,1,1)
%   plot(autoCorEmpA,'-k','linewidth',2)
% %  plot(autoCorA,'-r')

%   subplot(2,1,2)
%   plot(autoCorEmpB,'-k','linewidth',2)
% %  plot(autoCorB,'-r')

% end

% %keyboard
% %tol =3e-1;
% %assertVectorsAlmostEqual(xCorLastFirstA,0,'absolute',tol,0)
% %assertVectorsAlmostEqual(xCorLastFirstB,0,'absolute',tol,0)
