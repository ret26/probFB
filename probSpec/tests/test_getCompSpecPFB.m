
function test_suite = test_getCompSpecPFB
  initTestSuite;

  % Tests: function [Freqs,Spec] =
  % getCompSpecPFB(lams,varx,om,NumFreqs,RngFreqs)
  
  % Could do with testing higher order processes here i.e. beyond
  % second order
  %
  % I haven't tested here that the spectrum has the correct scale
  % (the tests focus on shape). However, the function sampleAR2FB
  % uses the spectrum to sample from an AR process and the tests
  % verify that the scale is correct.
  
function test_compareWithSample

dispFigs = 1;
lam = 0.8;
varx = 1/(1-lam.^2);
om = pi/10;

numAves = 1000;
NumFreqs = 100;
RngFreqs = [-1/2,1/2];
[Freqs,Spec] = getCompSpecPFB(lam,varx,om,NumFreqs,RngFreqs);

specY = zeros(NumFreqs,1);
for n=1:numAves

  T = NumFreqs; N = 1; 
  x1 = sampleARTau(lam,varx,T,N);
  x2 = sampleARTau(lam,varx,T,N);

  z = x1 + i*x2; 
  y = exp(i*om*[1:T]').*z;
  
  specYCur =(abs(fft(y)).^2);
  specY = specY+specYCur;
end

specNorm1 = specY/sqrt(sum(specY.^2));
specNorm1 = [specNorm1(T/2+2:end);specNorm1(1:T/2+1)];
specNorm2 = Spec/sqrt(sum(Spec.^2));

if dispFigs==1
  figure
  hold on
  plot(Freqs,specNorm2,'.k')
  plot(Freqs,specNorm1,'.r')
  legend('parametric','samples')
end

tol = 1e-1;
assertVectorsAlmostEqual(max(abs(specNorm1-specNorm2')),0,'absolute',tol,0)

