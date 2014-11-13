%function test_suite =  test_getGPSEisoSpec
function test_getGPSEisoSpec
  %initTestSuite;
  test_correctSpectrum
  test_correctDerivatives
  test_correctAutoCorrelation
end

% test 
% function fftCov = getGPSESpec(len,T);

function test_correctSpectrum

dispFigs=1;
randn('state',1)
var = 1.5;
len = 5.32;

T = 1000;

fftCov1 = getGPSEisoSpec(var,len,T);



dts = [0:T/2,[T/2-1:-1:1]]';
autoCor = var*exp(-1/(2*len^2)*(dts).^2);
fftCov2 = fft(autoCor);

if dispFigs==1
  figure;
  hold on
  plot(abs(fftCov2),'-k','linewidth',2)
  plot(fftCov1,'-r','linewidth',1)
end

tol = 1e-5;
assertVectorsAlmostEqual(fftCov1,fftCov2,'absolute',tol,0)
end

function test_correctDerivatives

dispFigs=1;
randn('state',1)
var = 1.5;
len = 5.32;

T = 1001;

[fftCov1,dfftCov1] = getGPSEisoSpec(var,len,T);
delta = 1e-6;
fftCov21 = getGPSEisoSpec(exp(log(var)+delta),len,T);
fftCov22 = getGPSEisoSpec(var,exp(log(len)+delta),T);

dfftCov22 = (fftCov22-fftCov1)/delta;
dfftCov21 = (fftCov21-fftCov1)/delta;


if dispFigs==1
  figure;
  hold on
  plot(dfftCov21,'-k','linewidth',2)
  plot(dfftCov1(:,1),'+k','linewidth',1)
  plot(dfftCov22,'-b','linewidth',2)
  plot(dfftCov1(:,2),'+b','linewidth',1)
end

tol = 1e-1;
assertVectorsAlmostEqual(dfftCov1(:,2),dfftCov22,'absolute',tol,0)
assertVectorsAlmostEqual(dfftCov1(:,1),dfftCov21,'absolute',tol,0)
end

function test_correctAutoCorrelation

dispFigs=1;
randn('state',1)
var = 1.5;
len = 5.32;

T = 1000;

fftCov = getGPSEisoSpec(var,len,T);

autoCor1 = ifft(fftCov);

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
end


