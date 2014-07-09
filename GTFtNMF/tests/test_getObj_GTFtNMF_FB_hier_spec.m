function test_suite = test_getObj_GTFtNMF_FB_hier_spec
  initTestSuite;

% Tests: 
%
% 
% function [Obj,dObj] = getObj_GTFtNMF_FB_hier_spec(zlogwspec,y, ...
% 						   vary,lenx,mux,varx,preZ,tol);

function test_check_gradient

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Forward model Settings

D = 40;
K = 40;
T = 20;
y = randn(T,1);
Z = randn(T,K);
logW = randn(K,D);
lamv = rand(1,D);
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

zlogwspec = [Z(:);logW(:);log(varv)';log(lamv./(1-lamv))';log((cos(omv')+1)./(1-cos(omv')))];


[Obj,dObj] =  getObj_GTFtNMF_FB_hier_spec(zlogwspec,y,vary,lenx,mux,varx,preZ,tol);

d=checkgrad('getObj_GTFtNMF_FB_hier_spec',zlogwspec,delta,y,vary, ...
	    lenx,mux,varx,preZ,tol)

tol = 1e-6;
assertVectorsAlmostEqual(d,0,'absolute',tol,0)


