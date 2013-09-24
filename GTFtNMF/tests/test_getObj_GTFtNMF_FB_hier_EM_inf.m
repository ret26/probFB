function test_suite = test_getObj_GTFtNMF_FB_hier_EM_inf
  initTestSuite;

% Tests: 
%
% 
% function [Obj,dObj] =
% getObj_GTFtNMF_FB_hier_EM_inf(z,W,VV,RVV,lamv,varv,lenx,mux,varx,preZ,tol);
%

function test_check_gradient

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Forward model Settings

D = 4;
K = 3;
T = 20;
VV = rand(T,D);
RVV = rand(T-1,D);
Z = randn(T,K);
W = exp(randn(K,D));
lamv = [0.8,0.3,0.6,0.9];
varv = exp(randn(1,D));

lenx = exp(randn(1,K));
mux = randn(1,K);
varx = exp(randn(1,K));
preZ = rand;
tol = 7;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TEST

delta = 1e-5;

[Obj,dObj] = getObj_GTFtNMF_FB_hier_EM_inf(Z(:),W,VV,RVV,lamv,varv,lenx,mux,varx,preZ,tol);

d=checkgrad('getObj_GTFtNMF_FB_hier_EM_inf',Z(:),delta,W,VV,RVV,lamv,varv,lenx,mux,varx,preZ,tol)

tol = 1e-6;
assertVectorsAlmostEqual(d,0,'absolute',tol,0)


