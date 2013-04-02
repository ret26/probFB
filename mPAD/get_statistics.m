function statistics = get_statistics(y1,varargin)
  
% function get_statistics(y1)
% or
% function get_statistics(y1,y2,fs,n_filts,flim,df,maxLag)
%
% returns the following statistics for the signal y1
% 
% INPUTS
% y1 = signal 1 [T,1]
% y2 = signal 2 [T,1]
% optional: 
% fs = sampling rate
% n_filts = number of filters
% flim = [min filter centre frequency, max filter centre frequency]
% df = filter bandwidth expressed as fraction of centre frequency
% maxLag = lag in samples to compute the auto correlation over 
%
% OUTPUTS
% statistics = structure of statistics
%     pg = periodogram [NumFreq,1]
%     pg_freqs = frequencies for the periodogram [NumFreq,1]
%     channel_freqs = centre-frequencies of the channels [n_filts,1]
%     covA = covariance between Hilbert envelopes in different channels [n_filts,n_filts]
%     cov_logA = covariance between log-Hilbert envelopes in
%               different channels [n_filts,n_filts]
%     ACorr = autocorrelation of the Hilbert envelopes in each
%             channel [D,2*maxLag+1]
%     logACorr = autocorrelation of the log Hilbert envelopes in
%                each channel [D,2*maxLag+1]
%     ekurtA = excess kurtosis of Hilbert envelopes in each channel [D,1]
%     ekurtlogA = excess kurtosis of log Hilbert envelopes in each channel [D,1]


% set the sampling rate
if nargin>2 
  fs = varargin{1};
else
  fs = 1;
end

% number of filters
if nargin>3 
  n_filts = varargin{2};
else
  n_filts = 30;
end

% min and max frequencies of the filters
if nargin>4 
  FLim = varargin{3};
else
  FLim = [2/100,10/22]; 
end

% width of the processes wrt centre-frequency
if nargin>5 
  dfFrac = varargin{4};
else
  dfFrac = 1/10; 
end

% set the lag for the auto correlation
if nargin>6 
  maxLag = varargin{5};
else
  maxLag = 500;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
% compare spectra (peridograms)
NumFreq = 1000;
freqs = linspace(0,fs/2,NumFreq);
[pg1,varpg] = welchMethod(y1,NumFreq,20); 

statistics.pg = pg1;
statistics.pg_freqs = freqs;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% envelope comparisons

% define filters for envelope comparisons

mVar = ones(n_filts,1);
fmax = logspace(log10(FLim(1)),log10(FLim(2)),n_filts)';
[om,lamx] = freq2probSpec(fmax,fmax*dfFrac,1);
varx = mVar.*(1-lamx.^2);
statistics.channel_freqs = fmax;

% apply filters and extract Hilbert envelopes
Z1 = probFB(y1,lamx,varx,om,0);
A1 = abs(Z1');

% statistics of amplitudes

% covariance of the amplitudes and log amplitudes
statistics.covA = cov(A1);
statistics.cov_logA = cov(log(A1));


% auto-correlation of the envelopes and log-envelopes
logaCorr1 = zeros(n_filts,maxLag*2+1);
aCorr1 = zeros(n_filts,maxLag*2+1);

for d=1:n_filts
  A1cur = A1(:,d)-mean(A1(:,d));
  aCorr1(d,:)= xcorr(A1cur,maxLag,'unbiased');
  
  logA1 = log(A1(:,d));
  logA1 = logA1 - mean(logA1);
  logaCorr1(d,:)= xcorr(logA1,maxLag,'unbiased');
end

statistics.logACorr = logaCorr1;
statistics.ACorr = aCorr1;

% kurtosis
ekurtA = zeros(n_filts,1);
ekurtlogA = zeros(n_filts,1);
for d=1:n_filts
  A4 = mean((A1(:,d)-mean(A1(:,d))).^4);
  A2 = mean((A1(:,d)-mean(A1(:,d))).^2);
  ekurtA(d) = A4./(A2^2)-3;

  logA1 = log(A1(:,d));
  logA4 = mean((logA1-mean(logA1)).^4);
  logA2 = mean((logA1-mean(logA1)).^2);
  ekurtlogA(d) = logA4./(logA2^2)-3;
end

statistics.ekurtA = ekurtA;
statistics.ekurtlogA = ekurtlogA;
