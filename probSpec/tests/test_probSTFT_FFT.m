function test_suite = test_probSTFT_FFT
  initTestSuite;

% Tests: 
%
% function S = probSTFT_FFT(y,lamx,varx,om,vary)
%

function testCompare_FFT_And_Kalman

% Compare to the Kalman filtering method

T = 200;
D = 2;
dispFigs=1;

lamx = [.9;.7];
varx = rand(length(lamx),1);
om = [1/10,1/20]';

vary = rand;

[y,Z] = samplePSTFT(lamx,varx,om,vary,T);
%y = randn(T,1);

[Z1,covZ1] = probSTFT(y,lamx,varx,om,vary);
Z2 = probSTFT_FFT(y,lamx,varx,om,vary);

if dispFigs==1
  figure
  subplot(2,1,1)
  hold on
  title('real')
  plot(real(Z1)','-k','linewidth',2)
  plot(real(Z2)','-r')
  legend('kalman','fft')
  
  subplot(2,1,2)
  hold on
  title('imaginary')
  plot(imag(Z1)','-k','linewidth',2)
  plot(imag(Z2)','-r')
  
  error=mean((abs(Z1(:)-Z2(:))).^2)
end

tol = 1e-3;
tEdge = 60;
assertVectorsAlmostEqual(Z1(:,tEdge:T-tEdge),Z2(:,tEdge:T-tEdge),'absolute',tol,0)

