clear;

pat = genpath('~/Synchronised/probFB');
addpath(pat);
addpath ~/Synchronised/Optimisation/

T = 10000;

% D = 3;
% CF = [1/32;1/20;1/10];
% DF = [1/128;1/50;1/100];
% varMa = [1;4;2];

D = 4;
CF = [1/32;1/20;1/10;1/15];
DF = [1/128;1/100;1/100;1/500];
varMa = [1;4;2;2];

vary = 0;

[LamTr,VarTr] = freq2AR2(CF,DF,varMa);

[y,X] = sampleAR2FB(LamTr,VarTr,vary,T);

profile on

[LamEst,VarEst,Info] = fitAR2FB(y,D);

profile off
profile viewer

[CFEst,dFEst,mVar] = AR22freq(LamEst,VarEst);

