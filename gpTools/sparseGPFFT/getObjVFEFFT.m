function [f,varargout] = getObjVFEFFT(params,specy,kernel,M)
	T = length(specy);
	hypers = params(1:end-1);
	vary = exp(params(end));
	[fftCov,dfftCov] = getVFESpecHelper(kernel,hypers,T);

	halfM = floor(M/2);
	indM  = [1:halfM T-halfM+1:T];
	indNotM = [(halfM+1):(T-halfM)];

	% objective function in frequency domain
	f1 = 1/2*sum(log(fftCov(indM)+vary)) + 1/(2*T)*sum(specy(indM)./(fftCov(indM)+vary));
	f2 = 1/2*(T-M)*log(vary);
	f3 = 1/(2*T)*sum(specy(indNotM)./vary);
	f4 = 1/(2*vary)*sum(fftCov(indNotM));
	f = f1+f2+f3+f4;

	if nargout>1
		df1 = 1./(2*(fftCov(indM)+vary)) - 1/(2*T)*specy(indM)./((fftCov(indM)+vary).^2);
		D = length(params);
		df = zeros(D,1);
		for i=1:D-1
			df(i) = sum(df1.*dfftCov(indM,i)) + 1/(2*vary)*sum(dfftCov(indNotM,i));
		end
		df(D) = sum(df1)*vary + 1/2*(T-M) - f3 - f4;
		varargout{1} = df;
	end
end
