function plot_statistics(statistics,varargin)

% function plot_statistics(statistics)
%    or
% function plot_statistics(statistics1,statistics2)
% 
% Plot the summary statistics returned by the function
% get_stastistics. If two statistics structures are input to the
% function they get displayed side by side for comparison.
%
% INPUT
% statistics = structure of statistics
%     pg = periodogram [NumFreq,1]
%     covA = covariance between Hilbert envelopes in different channels [n_filts,n_filts]
%     covlogA = covariance between log-Hilbert envelopes in
%               different channels [n_filts,n_filts]
%     ACorr = autocorrelation of the Hilbert envelopes in each
%             channel [D,2*maxLag+1]
%     logACorr = autocorrelation of the log Hilbert envelopes in
%                each channel [D,2*maxLag+1]
%     ekurtA = excess kurtosis of Hilbert envelopes in each channel [D,1]
%     ekurtlogA = excess kurtosis of log Hilbert envelopes in each channel [D,1]
% optional:
% statistics2 = a second structure containing fields identical to
%     the ones listed above

if nargin>1
  two_statistics = true;
  statistics2 = varargin{1};
else
  two_statistics = false;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% periodograms

figure
hold on
plot(statistics.pg_freqs,statistics.pg,'k')

if two_statistics
 plot(statistics2.pg_freqs,statistics2.pg,'-r')
 legend('y1','y2','location','northeast')
end
set(gca,'yscale','log')

title('periodograms')
xlabel('frequency /Hz')
ylabel('PSD')
drawnow;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% power in each channel
figure

subplot(2,1,1)
hold on
plot(statistics.channel_freqs,diag(statistics.covA),'-k')

if two_statistics
 plot(statistics.channel_freqs,diag(statistics2.covA),'-r')
 legend('y1','y2','location','northeast')
end
xlabel('channel centre freq /Hz')
ylabel('envelope power')


subplot(2,1,2)
hold on
plot(statistics.channel_freqs,diag(statistics.cov_logA),'-k')

if two_statistics
 plot(statistics.channel_freqs,diag(statistics2.cov_logA),'-r')
 legend('y1','y2','location','northeast')
end
xlabel('channel centre freq /Hz')
ylabel('log-envelope power')
drawnow;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% covariances

figure

if two_statistics
  subplot(2,2,1)
else
  subplot(2,1,1)
end

hold on
title('y1: covariance amp (diagonals removed)')
covplot = statistics.covA-diag(diag(statistics.covA));
surf(statistics.channel_freqs,statistics.channel_freqs,covplot,'edgecolor','none')
axlim = [min(statistics.channel_freqs),max(statistics.channel_freqs)];
set(gca,'xlim',axlim,'ylim',axlim)
set(gca,'xscale','log','yscale','log')
colorbar;

if two_statistics
  subplot(2,2,2)
else
  subplot(2,1,2)
end

hold on
title('y1: covariance log amp (diagonals removed)')
covplot = statistics.cov_logA-diag(diag(statistics.cov_logA));
surf(statistics.channel_freqs,statistics.channel_freqs,covplot,'edgecolor','none')
axlim = [min(statistics.channel_freqs),max(statistics.channel_freqs)];
set(gca,'xlim',axlim,'ylim',axlim)
set(gca,'xscale','log','yscale','log')
colorbar;

if two_statistics
  
  subplot(2,2,3)
  hold on
  title('y2: covariance amp (diagonals removed)')
  covplot = statistics2.covA-diag(diag(statistics2.covA));
  surf(statistics2.channel_freqs,statistics2.channel_freqs,covplot,'edgecolor','none')
  axlim = [min(statistics2.channel_freqs),max(statistics2.channel_freqs)];
  set(gca,'xlim',axlim,'ylim',axlim)
  set(gca,'xscale','log','yscale','log')
  colorbar;

  subplot(2,2,4)
  hold on
  title('y2: covariance log amp (diagonals removed)')
  covplot = statistics2.cov_logA-diag(diag(statistics2.cov_logA));
  surf(statistics2.channel_freqs,statistics2.channel_freqs,covplot,'edgecolor','none')
  axlim = [min(statistics2.channel_freqs),max(statistics2.channel_freqs)];
  set(gca,'xlim',axlim,'ylim',axlim)
  set(gca,'xscale','log','yscale','log')
  colorbar;
end
drawnow;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% kurtoses

figure
subplot(2,1,1)
hold on
plot(statistics.channel_freqs,statistics.ekurtA,'-k')

if two_statistics
 plot(statistics.channel_freqs,statistics2.ekurtA,'-r')
 legend('y1','y2','location','northeast')
end
xlabel('channel centre freq /Hz')
ylabel('envelope excess kurtosis')

subplot(2,1,2)
hold on
plot(statistics.channel_freqs,statistics.ekurtlogA,'-k')

if two_statistics
 plot(statistics.channel_freqs,statistics2.ekurtlogA,'-r')
 legend('y1','y2','location','northeast')
end
xlabel('channel centre freq /Hz')
ylabel('log-envelope excess kurtosis')
drawnow;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% autocorrelation

figure

n_filt = length(statistics.channel_freqs);

k = ceil(sqrt(n_filt));

for n=1:n_filt
  subplot(k,k,n)
  hold on
  plot(statistics.ACorr(n,:),'-k','linewidth',2)
  if two_statistics
    plot(statistics2.ACorr(n,:),'-r')
  end

  if n==1
    
  end
end
drawnow;
