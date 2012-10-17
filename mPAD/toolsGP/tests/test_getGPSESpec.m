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

dts = [0:T/2,[T/2-1:-1:1]]';

autoCor2 = var*exp(-1/(2*len^2)*(dts).^2);


if dispFigs==1
  figure;
  hold on
  plot(autoCor2,'-k','linewidth',2)
  plot(autoCor1,'-r','linewidth',1)
end

tol = 1e-5;
assertVectorsAlmostEqual(autoCor1,autoCor2,'absolute',tol,0)


function test_correctAutoCorrelation_2D

dispFigs=1;
randn('state',1)
var = 1.5;
len = [5.32;4.12];

T = 1000;

fftCov = getGPSESpec(len,T);

autoCor1A = var*ifft(fftCov(:,1));
autoCor1B = var*ifft(fftCov(:,2));

dts = [0:T/2,[T/2-1:-1:1]]';

autoCor2A = var*exp(-1/(2*len(1)^2)*(dts).^2);
autoCor2B = var*exp(-1/(2*len(2)^2)*(dts).^2);


if dispFigs==1
  figure;
  subplot(2,1,1)
  hold on
    
  plot(autoCor2A,'-k','linewidth',2)
  plot(autoCor1A,'-r','linewidth',1)

  subplot(2,1,2)
  hold on
  plot(autoCor2B,'-k','linewidth',2)
  plot(autoCor1B,'-r','linewidth',1)

end

tol = 1e-5;
assertVectorsAlmostEqual(autoCor1A,autoCor2A,'absolute',tol,0)
assertVectorsAlmostEqual(autoCor1B,autoCor2B,'absolute',tol,0)
