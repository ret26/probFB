function test_suite = test_GTFtNMF_EM_inf
  initTestSuite;

% Tests: 
%
% [W,H,mnV,covV,info] = GTFtNMF_EM_inf(y,W,H,lamv,varv,omv,vary,lenx,mux,varx,varargin)
%

function test_recover_true

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Forward model Settings
dispFigs = 1;

T = 500;
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

vary = 1e-3;
tol = 11;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TEST
[y,V,H,Z] = randGTFtNMF(W,lamv,varv,omv,vary,lenx,mux,varx,tol,T);

A = sqrt(1/2*H*W);

HInit = exp(randn(T,K));
WInit = rand(K,D);


[HEst,mnV,covV,info] = GTFtNMF_EM_inf(y,W,HInit,lamv,varv,omv,vary, ...
					lenx,mux,varx);

if dispFigs==1
figure
for k=1:K
  subplot(K,1,k)
  hold on
  plot(H(:,k),'-k')
  plot(HEst(:,k),'-r')
  set(gca,'yscale','log')
end

end

assertTrue(sum(abs(HInit(:)-H(:)))>sum(abs(H(:)-HEst(:))),'temporal basis functions not closer to true')

