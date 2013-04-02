function compareStatistics(y1,y2,X1a,X1b,X2a,X2b,Aa,Ab,fs)

% function compareStatistics(X1a,X1b,Aa,Ab)
%
% Compares the statistics of two sounds. E.g. one might be taken from
% a training sound and the other might be produced by sampling from a
% generative model


[T,D] = size(Aa); 

figure
subplot(2,3,1)
plot(var(real(X1a)),'-r','linewidth',2); hold on; plot(var(real(X1b)));
legend('data','sampled')
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


figure

ax_lim = [1/2,D+1/2];
subplot(2,2,1)
hold on
title('covariance data amp')
imagesc(covAa-diag(diag(covAa)))
set(gca,'xlim',ax_lim,'ylim',ax_lim)

subplot(2,2,2)
hold on
title('covariance sampled amp')
imagesc(covAb-diag(diag(covAb)))
set(gca,'xlim',ax_lim,'ylim',ax_lim)

subplot(2,2,3)
hold on
title('covariance data log-amp')
imagesc(cov_logAa-diag(diag(cov_logAa)))
set(gca,'xlim',ax_lim,'ylim',ax_lim)

subplot(2,2,4)
hold on
title('covariance sampled log-amp')
imagesc(cov_logAb-diag(diag(cov_logAb)))
set(gca,'xlim',ax_lim,'ylim',ax_lim)



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

figure
subplot(2,1,1)
imagesc(log(Aa)'-mean(log(Aa(:))))

subplot(2,1,2)
imagesc(log(Ab)'-mean(log(Ab(:))))


[T,D] = size(Aa);

figure;

maxLag = 500;

rootD = ceil(sqrt(D));



for d=1:D
  
  logAa = log(Aa(:,d));
  logAa = logAa - mean(logAa);

  logAb = log(Ab(:,d));
  logAb = logAb - mean(logAb);

  aCorr= xcorr(logAa,maxLag);
  bCorr= xcorr(logAb,maxLag);

%  aCorr = aCorr/max(aCorr);
%  bCorr = bCorr/max(bCorr);
  
  subplot(rootD,rootD,d)
  hold on
  plot(aCorr,'-k','linewidth',2)
  plot(bCorr,'-r')
  

end


