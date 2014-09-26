function [mf, vf] = denoiseSparseGP(y,oriSpec,vary,M)
    sparseSpec = getSparseGPSpec(oriSpec,M);
    specy = fft(y);
    T = length(y);
    mf = real(ifft(sparseSpec./(sparseSpec+vary).*specy));
    vf = 1/T*sum(oriSpec - sparseSpec./(sparseSpec+vary).*sparseSpec); ...
end