function y = sampleGPSM(varx,lenx,freqx,T)
%sampleGPSM returns a function drawn from a GP with a spectral mixture
% kernel
%  Inputs: 
%    varx, lenx, freqx: hyperparameters of the kernel
%    T: length of signal to generate
%  Output:
%    y: signal
%
% Written by Richard Turner
% Mdified by Thang Bui
% Last modified: 9/2014

% TODO: make tau an optional param, 

tau = 0*max(lenx); % offset
Tx = 2^ceil(log2(T+tau));
fftCov = getGPSMSpec(varx,lenx,freqx,Tx);
wn = randn(Tx,1);

y = real(ifft(sqrt(fftCov).*fft(wn)));
y = y(1:T);

end