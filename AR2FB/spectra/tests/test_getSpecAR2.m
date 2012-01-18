
function test_suite = test_getSpecAR2
  initTestSuite;

  % Could do with testing higher order processes here i.e. beyond
  % second order
  %
  % I haven't tested here that the spectrum has the correct scale
  % (the tests focus on shape). However, the function sampleAR2FB
  % uses the spectrum to sample from an AR process and the tests
  % verify that the scale is correct.
  
function test_compareWithSample

dispFigs = 0;
lam = [1,-0.5];
varx = 1;
tau = length(lam);

numAves = 1000;
NumFreqs = 50;
RngFreqs = [0,1/2];
[Freqs,SpecTau] = getSpecAR2(lam,varx,NumFreqs,RngFreqs);

specX = zeros(NumFreqs,1);
for n=1:numAves

  T = NumFreqs*2; N = 1; 
  X = sampleARTau(lam,varx,T,N);

  specXCur =(abs(fft(X)).^2);
  specX = specX+specXCur(1:NumFreqs);
end

specNorm1 = specX/sqrt(sum(specX.^2));
specNorm2 = SpecTau/sqrt(sum(SpecTau.^2));

if dispFigs==1
  figure
  hold on
  plot(Freqs,specNorm1,'.r')
  plot(Freqs,specNorm2,'-k')
end

tol = 1e-1;
assertVectorsAlmostEqual(max(abs(specNorm1-specNorm2')),0,'absolute',tol,0)

function test_compareWithARTau

dispFigs = 0;
lam = [1,-0.5];
varx = 1;
tau = length(lam);


NumFreqs = 50;
RngFreqs = [0,1/2];
[Freqs,SpecTau1] = getSpecAR2(lam,varx,NumFreqs,RngFreqs);
[Freqs,SpecTau2] = getSpecARTau(lam,varx,NumFreqs,RngFreqs);


if dispFigs==1
  figure
  hold on
  plot(Freqs,SpecTau1,'.r')
  plot(Freqs,SpecTau2,'-k')
end

tol = 1e-5;
assertVectorsAlmostEqual(SpecTau1,SpecTau2,'absolute',tol,0)

