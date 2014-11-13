function test_suite = test_getObjVFESMFFT
%initTestSuite;
test_checkgrad1
test_checkgrad2
test_checkgrad3

end


function test_checkgrad1
T = 20;
y = randn(T,1);
specy = abs(y).^2;

param = randn(4,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TEST

delta = 1e-5;

d=checkgrad('getObjVFESMFFT',param,delta,specy,20);

tol = 1e-5;
assertVectorsAlmostEqual(d,0,'absolute',tol,0)

end

function test_checkgrad2
T = 100;
y = randn(T,1);
specy = abs(y).^2;

param = randn(10,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TEST

delta = 1e-5;

d=checkgrad('getObjVFESMFFT',param,delta,specy,30);

tol = 1e-5;
assertVectorsAlmostEqual(d,0,'absolute',tol,0)
end

function test_checkgrad3
varx = 0.4;
lenx = 500;
freqx = 0.1;
vary = 0.05;
T = 40000;
y = sampleGPSM(varx,lenx,freqx,T);
figure, plot(y)
specy = abs(fft(y)).^2;
delta = 1e-5;
params = log([varx lenx freqx vary]');

d=checkgrad('getObjVFESMFFT',params,delta,specy,1000);

tol = 1e-5;
assertVectorsAlmostEqual(d,0,'absolute',tol,0)
end