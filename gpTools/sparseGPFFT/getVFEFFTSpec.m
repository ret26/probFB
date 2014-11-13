function spec = getSparseGPSpec(oriSpec,M)
    T = length(oriSpec);
    if (M>=T)
        error('Number of inducing points M must be smaller or equal T')
    end
    halfT = floor(T/2);
    halfM = floor(M/2);
    spec = zeros(T,1);
    spec(1:halfM) = oriSpec(1:halfM);
    spec(end-halfM+1:end) = oriSpec(end-halfM+1:end);
end