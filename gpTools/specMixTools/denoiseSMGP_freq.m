function [mf,vf] = denoiseSMGP_freq(vars,lens,freqs,vary,y)

T = length(y);

fy = fft(y);
fftCov1 = getGPSMSpec(vars,lens,freqs,T);
fftCov2 =  fftCov1 + vary;

mf = real(ifft(fftCov1./fftCov2.*fy));
vf = 1/T*sum(fftCov1 - fftCov1./fftCov2.*fftCov1);
end