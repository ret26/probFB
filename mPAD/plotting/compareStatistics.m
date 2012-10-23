function compareStatistics(y1,y2,X1a,X1b,X2a,X2b,Aa,Ab,fs)

% function compareStatistics(X1a,X1b,Aa,Ab)
%
% Compares the statistics of two sounds. E.g. one might be taken from
% a training sound and the other might be produced by sampling from a
% generative model


figure
subplot(2,3,1)
plot(var(real(X1a)),'-r','linewidth',2); hold on; plot(var(real(X1b)));
legend('sampled','data')
title('variance of X1')

subplot(2,3,2)
plot(var(X2a),'-r','linewidth',2); hold on; plot(var(X2b))
title('variance of X2')

subplot(2,3,3)
plot(var(Aa),'-r','linewidth',2); hold on; plot(var(Ab))
title('variance of A')

subplot(2,3,4)
covAa = cov(Aa);
covAb = cov(Ab);
hold on; plot(covAa(:),'-r','linewidth',2); plot(covAb(:));
title('covariance of A - comodulation')

subplot(2,3,5)
cov_logAa = cov(log(Aa));
cov_logAb = cov(log(Ab));
hold on; plot(cov_logAa(:),'-r','linewidth',2); plot(cov_logAb(:));
title('covariance of log(A) - comodulation')



NumFreq = 1000;
[pg1,varpg] = welchMethod(y1,NumFreq,20); 
[pg2,varpg] = welchMethod(y2,NumFreq,20); 

freqs = linspace(0,fs,NumFreq);

figure
hold on
plot(freqs,pg1,'k')
plot(freqs,pg2,'-r')
set(gca,'yscale','log')
legend('y1','y2','location','northeast')
xlabel('frequency /Hz')
ylabel('PSD')