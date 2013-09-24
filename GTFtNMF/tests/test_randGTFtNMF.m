function test_suite = test_randGTFtNMF
  initTestSuite;

% Tests: 
%
% function [y,V,Z] = randGTFtNMF(W,lamv,varv,omv,vary,lenx,mux,varx,tol,T);
%

function test_statistics

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Forward model Settings

T = 100000;
D = 4;
K = 3;

W = exp(randn(K,D));

mux = [1/2,3/4,1/3];
varx = [8,16,4];
lenx = [5,75,100];

lamv  = [0.99,0.995,0.95,0.99];
varv  = (1-lamv.^2);
omv = [pi/10,pi/20,pi/7,pi/5]';

vary = 0.0;
tol = 10;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TEST
[y,V,H,Z] = randGTFtNMF(W,lamv,varv,omv,vary,lenx,mux,varx,tol,T);

A = sqrt(1/2*H*W);

figure
plot(real(V(:,1))/max(real(V(:,1)))); 
hold on; 
plot(A(:,1)/max(A(:,1)),'-r')

logH = log(H);

muxEst = mean(logH);
varxEst = var(logH);
[mux;muxEst]
[varx;varxEst]
maxLag = 500;
acorr1Est = xcorr(logH(:,1)-muxEst(1),logH(:,1)-muxEst(1),maxLag,'unbiased')/varxEst(1);

dt = [-maxLag:1:maxLag];
acorr1 = exp(-1/(2*lenx(1)^2)*(dt).^2);

figure
hold on;
plot(dt,acorr1,'-k','linewidth',2)
plot(dt,acorr1Est,'-r')

%tol = 1e-1;
%assertVectorsAlmostEqual(muinf,mnH,'absolute',tol,0)



% mnH = mean(H,1);
% varH = var(H,1);
% meanlagH = (mean(H(2:end,:).*H(1:end-1,:),1)-mnH.^2)./varH;

% [muinf;mnH]
% [varinf;varH]
% [lam;meanlagH]

% tol = 1e-1;
% assertVectorsAlmostEqual(muinf,mnH,'absolute',tol,0)
% assertVectorsAlmostEqual(varinf,varH,'absolute',tol,0)
% assertVectorsAlmostEqual(lam,meanlagH,'absolute',tol,0)


