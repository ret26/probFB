function test_suite = test_getObj_GTFtNMF_FB_hier
  initTestSuite;

% Tests: 
%
% 
% function [Obj,dObj] =
% getObj_GTFtNMF_FB_hier(zlogw,y,lamv,varv,omv,vary,lenx,mux,varx,preZ,tol);
%

function test_check_gradient

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Forward model Settings

D = 4;
K = 3;
T = 20;
y = randn(T,1);
Z = randn(T,K);
logW = randn(K,D);
lamv = [0.8,0.3,0.6,0.9];
varv = exp(randn(1,D));
omv = rand(1,D)*pi;
vary = 0.001;
lenx = exp(randn(1,K));
mux = randn(1,K);
varx = exp(randn(1,K));
tol = 7;
preZ = 1.2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TEST

delta = 1e-5;

zlogw = [Z(:);logW(:)];

[Obj,dObj] =  getObj_GTFtNMF_FB_hier(zlogw,y,lamv,varv,omv,vary, ...
					   lenx,mux,varx,preZ,tol);

d=checkgrad('getObj_GTFtNMF_FB_hier',zlogw,delta,y,lamv,varv,omv, ...
	    vary,lenx,mux,varx,preZ,tol)

tol = 1e-6;
assertVectorsAlmostEqual(d,0,'absolute',tol,0)


