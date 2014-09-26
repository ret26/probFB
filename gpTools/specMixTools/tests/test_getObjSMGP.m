function test_suite = test_getObjSMGP
initTestSuite;
%test_checkgrad1
%test_checkgrad2
%test_checkgrad3
%test_ObjVsLenx
%test_ObjVsLenxVarx
end
% Tests:
%
% [f,df] = getObjSMGP(params,specy)

function test_checkgrad1
T = 20;
y = randn(T,1);
specy = abs(y).^2;

param = randn(4,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TEST

delta = 1e-5;

d=checkgrad('getObjSMGP',param,delta,specy);

tol = 1e-6;
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

d=checkgrad('getObjSMGP',param,delta,specy);

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
d=checkgrad('getObjSMGP',params,delta,specy);

tol = 1e-6;
assertVectorsAlmostEqual(d,0,'absolute',tol,0)
end
%{
function test_ObjVsLenx
varx = 0.4;
lenx = 500;
freqx = 0.2;
vary = 0.05;
T = 20000;
y = sampleGPSM(varx,lenx,freqx,T);
specy = abs(fft(y)).^2;

lens = [1:9 10:5:1000];
L = length(lens);
vals = zeros(L,1);
for l = 1:L
    params = log([varx lens(l) freqx vary]');
    vals(l) = getObjSMGP(params,specy);
end
keyboard
figure,
plot(lens,vals,'*')
xlabel('lengthscale'), ylabel('Obj value');

end
%}
%{
function test_ObjVsLenxVarx
varx = 0.4;
lenx = 500;
freqx = 0.2;
vary = 0.05;
T = 20000;
y = sampleGPSM(varx,lenx,freqx,T);
specy = abs(fft(y)).^2;

lens = [1:9 10:5:1000];
vars = [0.05:0.01:0.6];
L1 = length(lens);
L2 = length(vars);
vals = zeros(L1,L2);
for l1 = 1:L1
    for l2 = 1:L2
        fprintf('progress %d/%d, %d/%d\n',l1,L1,l2,L2);
        params = log([vars(l2) lens(l1) freqx vary]');
        vals(l1,l2) = getObjSMGP(params,specy);
    end
end
keyboard
figure,
plot3(lens,vars,vals,'*')
xlabel('lengthscale'), ylabel('varx'), zlabel('Obj value');

end
%}