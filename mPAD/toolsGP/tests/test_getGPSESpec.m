function test_suite =  test_getGPSESpec
  initTestSuite;

% test 
% function fftCov = getGPSESpec(len,T);

function test_correctAutoCorrelation

dispFigs=1;
randn('state',1)
var = 1.5;
len = 5.32;

T = 1000;

fftCov = getGPSESpec(len,T);

autoCor1 = var*ifft(fftCov);

dts = [0:T/2,[T/2-1:-1:1]];

autoCor2 = var*exp(-1/(2*len^2)*(dts).^2);


if dispFigs==1
  figure;
  hold on
  plot(autoCor2,'-k','linewidth',2)
  plot(autoCor1,'-r','linewidth',1)
end

tol = 1e-5;
assertVectorsAlmostEqual(autoCor1,autoCor1,'absolute',tol,0)
