function [f,df] = getObjSMGP(params,specy)
%getObjSMGP returns the neg log likelihood of GP model, computed 
% in frequency domain
%  Inputs:
%    params: hyperparameters of the spectral mixture kernel and noise variance
%    specy: signal spectrum
%  Outputs:
%    f: objective function value
%    df: vector of derivatives wrt hypers and noise variance
% See also trainSMGP_freq.m
%
% Written by Richard Turner and Thang Bui
% Last modified: 9/2014


specy = specy(:);
L = length(params);
K = floor(L/3); % number of components in the mixtures

T = length(specy);

% params of the covariance function
sig2s = exp(params(1:3:3*K));
lengs = exp(params(2:3:3*K));
freqs = exp(params(3:3:3*K));

% noise variance
vary = exp(params(end));

% obtain the fft of the covariance function
[fftCov,dfftCov] = getGPSMSpec(sig2s,lengs,freqs,T);
% objective function in frequency domain
f = 1/2*sum(log(fftCov+vary)) + 1/(2*T)*sum(specy./(fftCov+vary));
% 
dobjdspec = 1./(2*(fftCov+vary)) - 1/(2*T)*specy./((fftCov+vary).^2);

df = zeros(L,1);
for i = 1:3*K
    df(i) = sum(dobjdspec.*dfftCov(:,i));
end
df(L) = sum(dobjdspec)*vary; 

% normalise
f = f/T;
df = df/T;

debug = 0;
if (debug)
    fftCovTrue = getGPSMSpec(0.2,500,0.2,T);
    figure(1), plot(1:T,specy/T,1:T,fftCov,1:T,fftCovTrue), 
    legend('specy','specCov','true'), %pause(0.01)
    xlim(T*freqs + [-100 100]) 
    drawnow
end

end