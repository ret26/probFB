function params = initSMParams(s,K,wsize,ovlapPc,varargin)
%initSECosParams initialises spectral mixture parameters
% First, a smoothed spectrum of the signal is obtained to remove 
% noisy peaks. K peaks are chosen to initialise the params
%  Inputs:
%    s: signal
%    K: number of components in the kernel
%    wsize: window size in pwelch method [rectangular window]
%    ovlapPc: overlapping percentage
%  Output:
%    params: hyperparameters of the kernel
%  Optional input:
%    varargin{1} 1 to plot the spectrum
%
% Written by Thang Bui and Richard Turner
% Last modified: 9/2014


window  = ones(wsize,1); % rectangular window for now
overlap = round(ovlapPc*wsize);
[pxx,f] = pwelch(s,window,overlap,[],[],'onesided','power');
[pks,locs] = findpeaks(pxx,'SortStr','descend');


% the spectral density can be over-smoothed, there may be less than K peaks
% randomly choose some more
if length(pks) < K
    remaining   = setdiff(1:length(pxx),locs)';
    needed      = K-length(pks);
    chosenIdx   = randperm(length(remaining),needed);
    locs        = [locs; remaining(chosenIdx)];
end
locs    = locs(1:K);
fpeaks  = f(locs);
pks     = pxx(locs);
pks     = pks/min(pks)+1;
% TODO: may need to get subband data and choose lengthscales accordingly

sigstd = std(s);
params = zeros(K*3+1,1);
params(2:3:3*K)     = log(50);      % lengthscales=50 for now
params(1:3:3*K)     = log(pks/2);     % varx
params(3:3:3*K)     = log(fpeaks);  % frequencies
params(3*K+1)       = log(1/5*sigstd);  % vary

if length(varargin) == 1
    opt = varargin{1};
    if opt == 1 % plot
        pxxlog = 10*log10(pxx);
        figure(1), plot(f,pxxlog,'-b',fpeaks,pxxlog(locs),'+r')
        xlabel('freq'), ylabel('power (dB)');
    end
end

end