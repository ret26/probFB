
function test_suite = test_get_pSTFT_spec
  initTestSuite;

 
  
function test_compareWithSample

dispFigs = 0;
lam = 0.8;
varx = 1/(1-lam.^2);
om = pi/5;

numAves = 1000;
NumFreqs = 100;
freqs = linspace(-1/2,1/2,NumFreqs);
spec = get_pSTFT_spec(freqs,lam,varx,om);

specY = zeros(NumFreqs,1);
for n=1:numAves

  T = NumFreqs; N = 1; 
  x1 = sampleARTau(lam,varx,T,N);
  x2 = sampleARTau(lam,varx,T,N);

  z = x1 + i*x2; 
  y = real(exp(i*om*[1:T]').*z);
  
  specYCur =(abs(fft(y)).^2);
  specY = specY+specYCur;
end

specNorm1 = specY/sqrt(sum(specY.^2));
specNorm1 = [specNorm1(T/2+2:end);specNorm1(1:T/2+1)];
specNorm2 = spec/sqrt(sum(spec.^2));

if dispFigs==1
  figure
  hold on
  plot(freqs,specNorm2,'.k')
  plot(freqs,specNorm1,'.r')
  legend('parametric','samples')
end

tol = 1e-1;
assertVectorsAlmostEqual(max(abs(specNorm1-specNorm2')),0,'absolute',tol,0)



function test_use_to_sample

lam = 0.8;
varx = (1-lam.^2);
om = pi/5;

NumFreqs = 10000;
freqs = linspace(0,1/2,NumFreqs);
spec = get_pSTFT_spec(freqs,lam,varx,om);
    
spec = [spec,spec(end-1:-1:2)]';
T = 2*(NumFreqs-1);
y =  ifft(sqrt(spec).*fft(randn(T,1)));

tol = 1e-1;
assertVectorsAlmostEqual(var(y),1,'absolute',tol,0)
assertVectorsAlmostEqual(mean(y),0,'absolute',tol,0)

