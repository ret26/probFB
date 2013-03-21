function test_suite = test_getObj_mPAD_fixG_basis_func
  initTestSuite;

% Tests: 
%
% function  [Obj,dObj] = getObj_mPAD_fixG_basis_func(z2G,y,Bas,Params,GScale,tol);
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TEST

delta = 1e-5;

Z2init = randn(T,K);

z2g = [Z2init(:);Params.G(:)];

GScale = rand(K,1);
d=checkgrad('getObj_mPAD_fixG_basis_func',z2g,delta,y,Params,GScale,tol)

tol = 1e-6;
assertVectorsAlmostEqual(d,0,'absolute',tol,0)


