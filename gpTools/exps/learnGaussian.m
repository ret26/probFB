function learnGaussian
	a = input('Learn Gaussian filters: 1 for synthetic, 2 for speech: ');
	switch a
		case 1
			learnGaussian_synthetic
		case 2
			learnGaussian_speech
		otherwise
			disp('unknown option')
	end
end

function learnGaussian_synthetic
	varx = [1 1.5 0.4 0.8 1.2]';
	freqx = [0.05 0.1 0.15 0.2 0.3]';
	lenx = [100 50 100 80 50]';
	T = 10000;
	ytrue = sampleGPSM(varx,lenx,freqx,T);
	vary = 0.3;
	ynoisy = ytrue + sqrt(vary)*randn(T,1);

	fftCov = getGPSMSpec(varx,lenx,freqx,T);
	fftYclean = abs(fft(ytrue)).^2;
	fftYnoisy = abs(fft(ynoisy)).^2;


	t = (1:T)';
	figure(1), 
	subplot(2,1,1), 
	plot(t,ynoisy,'color',[0.5 0.5 0.5]), hold on,
	plot(t,ytrue,'b','LineWidth',2);
	hold off;
	subplot(2,1,2),
	plot(t,fftYnoisy/T,'color',[0.5 0.5 0.5]), hold on,
	plot(t,fftYclean/T,'b');
	plot(t,fftCov,'g');
	%set(gca,'Yscale','log')
	hold off;
	
	setup.numIts = 500;
	setup.progress_chunk = setup.numIts;
	K = 5;
	params = initSMParams(ynoisy,K,100,0.5);
    var1 = exp(params(1:3:3*K));
    var1 = 0.5*ones(K,1);
    %len1 = exp(params(2:3:3*K));
    len1 = 100*ones(K,1);
    freq1 = exp(params(3:3:3*K));
    vary1 = 0.2;
	[varxEst,lenxEst,freqxEst,varyEst,info] = trainSMGP_freq(ynoisy,K,setup,varx,lenx,freqx,vary)
	[mf,vf] = denoiseSMGP_freq(varxEst,lenxEst,freqxEst,varyEst,ynoisy);
	fftCovEst = getGPSMSpec(varxEst,lenxEst,freqxEst,T);

	figure(2),
	subplot(2,1,1), 
	plot(t,ynoisy,'color',[0.5 0.5 0.5]), hold on,
	plot(t,ytrue,'b');
	plot(t,mf,'r-',t,mf+2*sqrt(vf),'r--',t,mf-2*sqrt(vf),'r--')
	hold off;

	subplot(2,1,2),
	plot(t,fftYnoisy/T,'color',[0.5 0.5 0.5]), hold on,
	plot(t,fftYclean/T,'b');
	plot(t,fftCov,'g');
	plot(t,fftCovEst,'r')
	%set(gca,'Yscale','log')
	hold off;
	keyboard
end

function learnGaussian_speech
	[ytrue,fs] = audioread('../../tsgp/datasets/audio/74 - Sentences.wav');
	ytrue = ytrue(8000:35000);
	T = length(ytrue);
	ytrue = ytrue/std(ytrue);
	ytrue = ytrue - mean(ytrue);
	ynoisy = ytrue + 0.4*randn(T,1);
	
	fftYclean = abs(fft(ytrue)).^2;
	fftYnoisy = abs(fft(ynoisy)).^2;
	
	t = (1:T)';
	figure(1), 
	subplot(2,1,1), 
	plot(t,ynoisy,'color',[0.5 0.5 0.5]), hold on,
	plot(t,ytrue,'b','LineWidth',2);
	hold off;
	subplot(2,1,2),
	plot(t,fftYnoisy/T,'color',[0.5 0.5 0.5]), hold on,
	plot(t,fftYclean/T,'b');
	%set(gca,'Yscale','log')
	hold off;
	
	setup.numIts = 500;
	setup.progress_chunk = setup.numIts;
	[varxEst,lenxEst,freqxEst,varyEst,info] = trainSMGP_freq(ynoisy,50,setup)
	[mf,vf] = denoiseSMGP_freq(varxEst,lenxEst,freqxEst,varyEst,ynoisy);
	fftCovEst = getGPSMSpec(varxEst,lenxEst,freqxEst,T);

	figure(2),
	subplot(2,1,1), 
	plot(t,ynoisy,'color',[0.5 0.5 0.5]), hold on,
	plot(t,ytrue,'b');
	plot(t,mf,'r-',t,mf+2*sqrt(vf),'r--',t,mf-2*sqrt(vf),'r--')
	hold off;

	subplot(2,1,2),
	plot(t,fftYnoisy/T,'color',[0.5 0.5 0.5]), hold on,
	plot(t,fftYclean/T,'b');
	plot(t,fftCovEst,'r')
	set(gca,'Yscale','log')
	hold off;

	p1 = audioplayer(ytrue,fs);
	p2 = audioplayer(ynoisy,fs);
	p3 = audioplayer(mf,fs);

	keyboard
end

