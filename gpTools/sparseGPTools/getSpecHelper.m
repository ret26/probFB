function [fftCov,dfftCov] = getSpecHelper(kernel,params,T)
	if (strcmpi(kernel.name,'SE'))
		% TODO
	elseif (strcmpi(kernel.name,'SM'))
		K = length(params)/3;
		% params of the covariance function
		sig2s = exp(params(1:3:3*K));
		lengs = exp(params(2:3:3*K));
		freqs = exp(params(3:3:3*K));
		[fftCov,dfftCov] = getGPSMSpec(sig2s,lengs,freqs,T);
	else
		error('kernel not implemented yet')
	end
end