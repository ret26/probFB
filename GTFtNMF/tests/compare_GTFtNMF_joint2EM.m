% compares EM to joint optimisation for learning the spectral and
% temporal basis functions in GTFtEM

clear;

% Forward model Settings
T = 10000;
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
% TEST DATA
[y,V,H,Z] = randGTFtNMF(W,lamv,varv,omv,0,lenx,mux,varx,tol,T);

A = sqrt(1/2*H*W);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Probilistic filter bank to get envelopes

Z = probFB(y,lamv,varv,omv,0);
Asq = abs(Z)';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NMF to initialise features
HInit = exp(randn(T,K));
ks = ceil(T*rand(K,1));
WInit = Asq(ks,:);
varyNMF = zeros(T,D);

% initialise with regular NMF
Opts.restarts = 50;
Opts.numIts = 1000;
[WEst1,HEst1,info1] = nmf_fp(Asq,WInit,HInit,varyNMF,Opts);
AEst1 = sqrt(1/2*HEst1*WEst1);
HEst1 = HEst1+1e-5;
fastness = mean(diff(log(HEst1)).^2)./var(log(HEst1));

% now apply GTFNMF
[WEst1b,HEst1b,mnV1b,covV1b,info1b] = GTFNMF(y,WEst1,HEst1,lamv,varv,omv,vary);

% reorder in terms of temporal speed and remove zeros
fastness = mean(diff(log(HEst1b)).^2)./var(log(HEst1b));
[val,ind] = sort(fastness,'descend');
HEst1b = HEst1b(:,ind); 
WEst1b = WEst1b(ind,:);


% now compare
[WEst2,HEst2,mnV2,covV2,info2] = GTFtNMF_EM(y,WEst1b,HEst1b,lamv, ...
					    varv,omv,vary, lenx,mux,varx);

OptsJoint.numIts = 50;
OptsJoint.progress_chunk = 5;
[WEst3,HEst3,mnV3,covV3,info3] = GTFtNMF(y,WEst1b,HEst1b,lamv,varv,omv,vary, ...
				       lenx,mux,varx,OptsJoint);

AEst2 = sqrt(1/2*HEst2*WEst2);
AEst3 = sqrt(1/2*HEst3*WEst3);

for d=1:D
  figure
  hold on
  plot(real(V(:,d)),'-k')
  plot(A(:,d),'-r','linewidth',2)
  plot(AEst1(:,d),'-m','linewidth',2)
  plot(AEst2(:,d),'-b','linewidth',2)
  plot(AEst3(:,d),'-g','linewidth',2)
  legend('true filter','true amplitude','NMF', 'EM estimate','joint optimisation')
end

snr1 = 10*log10(mean(A(:).^2))-10*log10(mean((A(:)-AEst2(:)).^2));
snr2 = 10*log10(mean(A(:).^2))-10*log10(mean((A(:)-AEst3(:)).^2));

envelope_error = [snr1,snr2]

figure
hold on
title('spectrum of y')
specy = abs(fft(y));
plot(specy(1:T/2))

figure
hold on
plot(cumsum(info2.tim)/60,info2.lik(1:end-1),'-b')
plot(cumsum(info3.tim)/60,info3.lik,'-g')
legend('EM','Joint')
xlabel('time /mins')
ylabel('log-lik per data-point')

figure
for k=1:K
  subplot(K,1,k)
  hold on
  plot(H(:,k),'-r')
  plot(HEst2(:,k),'-b')
  plot(HEst3(:,k),'-g')
  set(gca,'yscale','log')
end

figure
for k=1:K
  subplot(K,1,k)
  hold on
  plot(H(:,k),'-r')
  plot(HEst1b(:,k),'-b')
  plot(HEst1(:,k),'-g')
  set(gca,'yscale','log')
end


% following only has meaning if the parameters are ordered correctly:
%tempral_basis_error = [mean((H(:)-HEst1(:)).^2),mean((H(:)-HEst2(:)).^2)]

figure
subplot(5,1,1)
imagesc(W)

subplot(5,1,2)
imagesc(WEst1)

subplot(5,1,3)
imagesc(WEst1b)

subplot(5,1,4)
imagesc(WEst2)

subplot(5,1,5)
imagesc(WEst3)


% following only has meaning if the parameters are ordered correctly:
%spectral_basis_error = [mean((W(:)-WEst1(:)).^2),mean((W(:)-WEst2(:)).^2)]

