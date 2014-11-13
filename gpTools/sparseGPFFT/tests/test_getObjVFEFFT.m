function test_suite = test_getObjSparseGP
%initTestSuite;
test_checkgrad1
test_checkgrad2
test_checkgrad3

end


function test_checkgrad1
T = 20;
y = randn(T,1);
specy = abs(y).^2;

kernel.name = 'SM';
kernel.K = 1;
param = randn(4,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TEST

delta = 1e-5;

d=checkgrad('getObjVFEFFT',param,delta,specy,kernel,10);

tol = 1e-6;
assertVectorsAlmostEqual(d,0,'absolute',tol,0)

end

function test_checkgrad2
T = 100;
y = randn(T,1);
specy = abs(y).^2;

kernel.name = 'SM';
kernel.K = 3;
param = randn(10,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TEST

delta = 1e-5;

d=checkgrad('getObjVFEFFT',param,delta,specy,kernel,30);

tol = 1e-6;
assertVectorsAlmostEqual(d,0,'absolute',tol,0)
end

function test_checkgrad3
varx = 0.4;
lenx = 500;
freqx = 0.2;
vary = 0.05;
T = 40000;
y = sampleGPSM(varx,lenx,freqx,T);
specy = abs(fft(y)).^2;
delta = 1e-5;
params = log([varx lenx freqx vary]');
kernel.name = 'SM';
kernel.K = 1;
d=checkgrad('getObjVFEFFT',params,delta,specy,kernel,500);

tol = 1e-6;
assertVectorsAlmostEqual(d,0,'absolute',tol,0)
end