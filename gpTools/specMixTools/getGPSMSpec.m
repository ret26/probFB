function [fftCov,varargout] = getGPSMSpec(vars,lens,freqs,T)
%getGPSMSpec gets spectrum of a spectral mixture kernel
%  Inputs:
%    vars: variances of components in the kernel
%    lens: lengthscales of components in the kernel
%    freqs: frequencies of components in the kernel
%    T: length of data
%  Outputs:
%    fftCov: spectrum of the kernel, evaluated at [0:T/2 -T/2:1]
%    dfftCov: optional, derivatives of the kernel wrt hypers
%
%  k(f) = sum_k sqrt(pi/2)*sig2k*lengk ...
%        *(exp(-2*pi^2*lengk^2*(f-freqk).^2) ...
%        + exp(-2*pi^2*lengk^2*(f+freqk).^2));
%
% See also getGPSESpec.m
% 
% Written by Richard Turner and Thang Bui
% Last modified: 9/2014


f = [0:ceil(T/2) -floor(T/2)+1:1:-1]'/T;

% deal with only one dimensional case for now
K = length(lens);
fftCovMat = zeros(T,K);
for k = 1:K
    sig2k = vars(k);
    lengk = lens(k);
    freqk = freqs(k);
    fftCovMat(:,k) = sqrt(pi/2)*sig2k*lengk ...
        *(exp(-2*pi^2*lengk^2*(f-freqk).^2) ...
        + exp(-2*pi^2*lengk^2*(f+freqk).^2));
end

fftCov = sum(fftCovMat,2);
% derivatives with respect to log of vars, lens and freqs
if nargout>1
    dfftCov = zeros(T,3*K);
    for k = 1:K
        sig2k = vars(k);
        lengk = lens(k);
        freqk = freqs(k);
        dfftCov(:,3*(k-1)+1) = fftCovMat(:,k); % d/dlogvarx_k
        M1 = -4*pi^2*lengk^2 * sqrt(pi/2)*sig2k*lengk*exp(-2*pi^2*lengk^2*(f-freqk).^2);
        M2 = -4*pi^2*lengk^2 * sqrt(pi/2)*sig2k*lengk*exp(-2*pi^2*lengk^2*(f+freqk).^2);
        dfftCov(:,3*(k-1)+2) = fftCovMat(:,k) ...
            + (f-freqk).^2.*M1 + (f+freqk).^2.*M2;
        dfftCov(:,3*(k-1)+3) = -freqk*(f-freqk).*M1 ...
            + freqk*(f+freqk).*M2;
    end
    varargout{1} = dfftCov;
end

end