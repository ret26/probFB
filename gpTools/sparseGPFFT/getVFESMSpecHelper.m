function [fftCovD,dfftCovD,fftCovY,dfftCovY] = getVFESMSpecHelper(hypers,T,M)
	% number of components in the mixture
	D = length(hypers)/3;
	% params of the covariance function
	sig2s = exp(hypers(1:3:3*D));
	lengs = exp(hypers(2:3:3*D));
	freqs = exp(hypers(3:3:3*D));

	% for each components, compute the fft of the covariance function
	% and its derivatives wrt its signal variance and lengthscale.
	fftCovD = zeros(T,D);
	dfftCovD = zeros(T,2,D);
	fftCovYMat = zeros(T,D);
	dfftCovY = zeros(T,3*D);
	for d = 1:D
		% get the spectrum for each component
		sig2d = sig2s(d);
		lengd = lengs(d);
		[fd,dfd] = getGPSEisoSpec(sig2d,lengd,T);
		fftCovD(:,d) = fd;
		dfftCovD(:,:,d) = dfd;

		% truncate the spectrum
		Md = M(d);
		if (Md>T)
        	error('Number of inducing points M must be smaller or equal T')
	    end
	    halfMd = floor(Md/2);
	    fdM = zeros(T,1);
	    fdM(1:halfMd) = fd(1:halfMd);
	    fdM(end-halfMd+1:end) = fd(end-halfMd+1:end);
	    dfdM = zeros(T,2);
	    dfdM(1:halfMd,:) = dfd(1:halfMd,:);
	    dfdM(end-halfMd+1:end,:) = dfd(end-halfMd+1:end,:);
	    

	    % compute the contribution of the dth component to the first row
	    % of P_yy
	    freqd = freqs(d);
	    wdt = 2*pi*freqd*(0:T-1)';
	    coswdt = cos(wdt);
	   	fftCovYMat(:,d) = coswdt.*ifft(fdM);

	   	dfftCovY(:,3*d-2) = fft(coswdt.*ifft(dfdM(:,1)));
	   	dfftCovY(:,3*d-1) = fft(coswdt.*ifft(dfdM(:,2)));
	   	dfftCovY(:,3*d) = fft(-wdt.*sin(wdt).*ifft(fdM));

	end
	fftCovYSum = sum(fftCovYMat,2);
	fftCovY = fft(fftCovYSum); 
end