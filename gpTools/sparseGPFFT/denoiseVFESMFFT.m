function [m, v] = denoiseVFESMFFT(y,params,M,uOrF,d)
	hypers = params(1:end-1);
	vary = exp(params(end));
	T = length(y);
	D = length(hypers)/3;
	if numel(M) == 1 && D>1
		M = repmat(M,D,1);
	end
	[fftCovD,dfftCovD,fftCovY,dfftCovY] = getVFESMSpecHelper(hypers,T,M);
	GammaY = fftCovY + vary;
	SigmaYinvY = ifft(fft(y)./(fftCovY+vary));
	%t = [0:ceil(T/2) -floor(T/2)+1:1:-1]';
	t = [0:T-1]';
	if strcmpi(uOrF,'u')
		freqd = exp(hypers(3*d));
		m1 = fft(cos(2*pi*freqd*t).*SigmaYinvY);
		m2 = fft(sin(2*pi*freqd*t).*SigmaYinvY);
		halfMd = floor(M(d)/2);
		indM  = [1:halfMd T-halfMd+1:T];
		mu1 = ifft(m1(indM).*fftCovD(indM,d));
		mu2 = ifft(m2(indM).*fftCovD(indM,d));
		m = real([mu1 mu2]);
		keyboard
	elseif strcmpi(uOrF,'f')
		freqd = exp(hypers(3*d));
		m1 = fft(cos(2*pi*freqd*t).*SigmaYinvY);
		m2 = fft(sin(2*pi*freqd*t).*SigmaYinvY);
		halfMd = floor(M(d)/2);
		indM  = [1:halfMd T-halfMd+1:T];
		indNotM = [(halfMd+1):(T-halfMd)];
		fd = fftCovD(:,d);
		fd(indNotM) = 0;
		mu1 = ifft(m1.*fd);
		mu2 = ifft(m2.*fd);
		m = real([mu1 mu2]);
	elseif strcmpi(uOrF,'uall') % TODO
	elseif strcmpi(uOrF,'fall')	% TODO
    else
    	m = zeros(T,1);
    	for d = 1:D
	    	freqd = exp(hypers(3*d));
			m1 = fft(cos(2*pi*freqd*t).*SigmaYinvY);
			m2 = fft(sin(2*pi*freqd*t).*SigmaYinvY);
			halfMd = floor(M(d)/2);
			indM  = [1:halfMd T-halfMd+1:T];
			indNotM = [(halfMd+1):(T-halfMd)];
			fd = fftCovD(:,d);
			fd(indNotM) = 0;
			mu1 = ifft(m1.*fd);
			mu2 = ifft(m2.*fd);
			m = m + real(mu1.*cos(2*pi*freqd*t) + mu2.*sin(2*pi*freqd*t));
		end
	end
    v = zeros(size(m)); %TODO
end