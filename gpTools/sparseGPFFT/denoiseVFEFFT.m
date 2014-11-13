function [mf, vf] = denoiseVFEFFT(y,oriSpec,vary,M)
    sparseSpec = getVFEFFTSpec(oriSpec,M);
    specy = fft(y);
    T = length(y);
    mf = real(ifft(sparseSpec./(sparseSpec+vary).*specy));
    vf = 1/T*sum(oriSpec - sparseSpec./(sparseSpec+vary).*sparseSpec); ...
end