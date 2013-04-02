function errors = statistics_error(stats1,stats2,varargin)
  
% errors = statistics_error(stats1,stats2)
%   or
% errors = statistics_error(stats1,stats2,verb)
% 
% computes errors in the statistics and returns
% 
% INPUTS
% stats1 and stats2 = both structures of statistics
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
% verb = verbose output printing errors to screen if verb=true

fnames = fieldnames(stats1);

n_fields = length(fnames);
errors = zeros(n_fields,1);

for n=1:n_fields
  s1 = getfield(stats1,fnames{n});
  s2 = getfield(stats2,fnames{n});
  errors(n) = mean((s1(:)-s2(:)).^2);  
end

if nargin>2
  if varargin{1}
    disp('------- errors -------')
    disp(['spectra: ',num2str(errors(1))])
    disp(['covariance amplitudes: ',num2str(errors(4))])
    disp(['covariance log amplitudes: ',num2str(errors(5))])
    disp(['autocorrelation amplitudes: ',num2str(errors(7))])
    disp(['autocorrelation log amplitudes: ',num2str(errors(6))])
    disp(['excess kurtosis amplitudes: ',num2str(errors(8))])
    disp(['excess kurtosis log amplitudes: ',num2str(errors(9))])
  end
end