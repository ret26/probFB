function test_suite = test_GTFNMF_EM
  initTestSuite;

% Tests: 
%
% [W,H,mnV,covV,info] = GTFNMF_EM(y,W,H,lamv,varv,omv,vary,varargin)
%

function test_likelihood_increases


T = 100;
D = 5;
K = 3;

mux = [1,1,1];
varx = [3,3,3];
lenx = [10,20,30];

W = [1/3,0,1/3,0,1/3;
     0,1/3,0,2/3,0;
     1/3,1/3,1/3,0,0];

lamv  = ones(1,D)*0.995;
varv  = (1-lamv.^2);
omv = linspace(pi/4,pi/30,D);

vary = 1e-3; % only used in inference and not for generation
tol = 11;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% generate 
[y,V,H,Z] = randGTFtNMF(W,lamv,varv,omv,0,lenx,mux,varx,tol,T);

A = sqrt(1/2*H*W);

HInit = ones(T,K);%exp(randn(T,K));
WInit = rand(K,D);

[WEst1,HEst1,mnV,covV,info1] = GTFNMF_EM(y,WInit,HInit,lamv,varv,omv,vary);

tol = 1e-4;
assertVectorsAlmostEqual(sum(diff(info1.lik)<0),0,'absolute',tol,0)

