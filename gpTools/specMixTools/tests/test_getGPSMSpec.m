function test_suite = test_getGPSMSpec
initTestSuite
%test_correctSpectrum
%test_correctAutoCorrelation
end

function test_correctSpectrum
dispFigs=1;
randn('state',1)
vars = [1.5 0.2 2.1];
lens = [5.32 80.1 3.5];
freqs = [0.05 0.25 0.37];

T = 100000;
fftCov1 = getGPSMSpec(vars,lens,freqs,T);


dts = [0:T/2,[T/2-1:-1:1]]';
%dts = [0:T-1]';
K = 3;
autoCor = zeros(length(dts),1);
for k = 1:K
    autoCor = autoCor + vars(k)*exp(-1/(2*lens(k)^2)*(dts).^2).*cos(2*pi*freqs(k)*dts);
end
fftCov2 = fft(autoCor);

if dispFigs==1
  figure;
  hold on
  plot(abs(fftCov2),'-k','linewidth',2)
  plot(fftCov1,'-r','linewidth',1)
end

tol = 1e-3;
assertVectorsAlmostEqual(fftCov1,fftCov2,'absolute',tol,0)
end

function test_correctAutoCorrelation
dispFigs=1;
randn('state',1)

vars = [1.5 0.2 2.1];
lens = [5.32 80.1 3.5];
freqs = [0.05 0.25 0.37];

T = 100000;
fftCov1 = getGPSMSpec(vars,lens,freqs,T);
autoCor1 = ifft(fftCov1);


dts = [0:T/2,[T/2-1:-1:1]]';
%dts = [0:T-1]';
K = 3;
autoCor2 = zeros(length(dts),1);
for k = 1:K
    autoCor2 = autoCor2 + vars(k)*exp(-1/(2*lens(k)^2)*(dts).^2).*cos(2*pi*freqs(k)*dts);
end

if dispFigs==1
  figure;
  hold on
  plot(autoCor2,'-k','linewidth',2)
  plot(autoCor1,'-r','linewidth',1)
end
%max(autoCor1-autoCor2)
tol = 1e-5;
assertVectorsAlmostEqual(autoCor1,autoCor2,'absolute',tol,0)
end
