function [f,varargout] = getObjSparseGP(params,specy,M)
	T = length(specy);
	hypers = params(1:end-1);
	vary = exp(params(end));
	if numel(M) == 1 && length(hypers)/3>1
		M = repmat(M,length(hypers)/3,1);
	end
	[fftCovD,dfftCovD,fftCovY,dfftCovY] = getVFESMSpecHelper(hypers,T,M);
	
	D = length(M); % number of components in the mixture
	% objective function in frequency domain
	f1 = 1/2*sum(log(fftCovY+vary)) + 1/(2*T)*sum(specy./(fftCovY+vary));
	f2 = 0;
	for d = 1:D
		halfMd = floor(M(d)/2);
		indM  = [1:halfMd T-halfMd+1:T];
		indNotM = [(halfMd+1):(T-halfMd)];
		f2 = f2 + 1/(2*vary)*sum(fftCovD(indNotM,d));
	end
	f = f1+f2;

	% compute the derivatives
	if nargout>1
		df1 = 1./(2*(fftCovY+vary)) - 1/(2*T)*specy./((fftCovY+vary).^2);
		df = zeros(length(params),1);
		for d=1:D
			halfM = floor(M(d)/2);
			indM  = [1:halfM T-halfM+1:T];
			indNotM = [(halfM+1):(T-halfM)];
			df(d*3-2) = sum(df1.*dfftCovY(:,d*3-2)) + 1/(2*vary)*sum(dfftCovD(indNotM,1,d));
			df(d*3-1) = sum(df1.*dfftCovY(:,d*3-1)) + 1/(2*vary)*sum(dfftCovD(indNotM,2,d));
			df(d*3) = sum(df1.*dfftCovY(:,d*3));
		end

		df(end) = sum(df1)*vary - f2;
		df = real(df);
		varargout{1} = df/T;
	end
	f = real(f)/T;
end 