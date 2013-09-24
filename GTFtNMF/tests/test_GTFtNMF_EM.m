function test_suite = test_GTFtNMF_EM
  initTestSuite;

% Tests: 
%
% [W,H,mnV,covV,info] = GTFtNMF_EM(y,W,H,lamv,varv,omv,vary,lenx,mux,varx,varargin)
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

[WEst1,HEst1,mnV,covV,info1] = GTFtNMF_EM(y,WInit,HInit,lamv,varv,omv,vary, ...
				       lenx,mux,varx);

tol = 1e-4;
assertVectorsAlmostEqual(sum(diff(info1.lik)<0),0,'absolute',tol,0)

% AEst1 = sqrt(1/2*HEst1*WEst1);

% for d=1:D
%   figure
%   hold on
%   plot(real(V(:,d)),'-k')
%   plot(A(:,d),'-r','linewidth',2)
%   plot(AEst1(:,d),'-b','linewidth',2)
%   legend('true filter','true amplitude', 'EM estimate','joint optimisation')
% end

% envelope_error = [mean((A(:)-AEst1(:)).^2),mean((A(:)-AEst2(:)).^2)]

% figure
% hold on
% title('spectrum of y')
% specy = abs(fft(y));
% plot(specy(1:T/2))

% figure
% hold on
% plot(info1.lik,'-b')
% plot(info2.lik,'-g')
% legend('EM','Joint')

% figure
% for k=1:K
%   subplot(K,1,k)
%   hold on
%   plot(H(:,k),'-r')
%   plot(HEst1(:,k),'-b')
%   plot(HEst2(:,k),'-g')
%   set(gca,'yscale','log')
% end

% % following only has meaning if the parameters are ordered correctly:
% %tempral_basis_error = [mean((H(:)-HEst1(:)).^2),mean((H(:)-HEst2(:)).^2)]

% figure
% subplot(3,1,1)
% imagesc(W)

% subplot(3,1,2)
% imagesc(WEst1)

% subplot(3,1,3)
% imagesc(WEst2)

% % following only has meaning if the parameters are ordered correctly:
% %spectral_basis_error = [mean((W(:)-WEst1(:)).^2),mean((W(:)-WEst2(:)).^2)]


