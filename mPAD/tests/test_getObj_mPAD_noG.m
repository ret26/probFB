function test_suite = test_getObj_mPAD_noG
  initTestSuite;

% Tests: 
%
% function [Obj,dObj] = getObj_mPAD_noG(x2,y,fftCov,Params);
%

function test_check_gradient

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Forward model Settings

T = 20;
FSamp = 8000;
D = 3;
K = 2;
Dims = packDimsMPAD(D,K,T);

% Carrier Processes
mVar1 = ones(1,D);
CF = [400,750,1000]/FSamp;
DF = [40,50,80]/FSamp;

% Modulator Processes
Len2 = [5,10]; % time constant in samples of modulators
Mu2 = randn(D,1);
Var2 = 9*ones(K,1);

G = [1,0;
     0,1;
     1,1];

vary = 0.001;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GENERATE DATA

[om,Lam1,Var1] = freq2probSpec(CF,DF,mVar1);

Params = packParamsMPAD(G,Lam1,Var1,Len2,Var2,Mu2,vary,om);
[y,X1,X2,Amp] = sampleMPAD(Params,Dims);

tau = 5*max(Len2); % offset added to avoid wrap-around effects 
Tx = 2^ceil(log2(T+tau)); % make the duration a power of 2 so fft is faster
fftCov = getGPSESpec(Len2,Tx);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TEST

delta = 1e-5;

X2init = randn(Tx,K);

x2 = X2init(:);

GScale = rand(K,1);
d=checkgrad('getObj_mPAD_noG',x2,delta,y,fftCov,Params)

tol = 1e-6;
assertVectorsAlmostEqual(d,0,'absolute',tol,0)


