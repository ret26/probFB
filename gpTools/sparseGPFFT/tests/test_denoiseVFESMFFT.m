function test_denoiseVFESMFFT

	test_notraining_y_1
    test_notraining_y_2
	test_notraining_f
	test_notraining_u
end

function test_notraining_y_1
	T = 1000;
	varx = 1;
    lenx = 20;
    vary = 0.1;
    freqx = 0.1;
    ytrue = sampleGPSE(varx,lenx,T).*cos(2*pi*freqx*(0:T-1)') ...
    	+ sampleGPSE(varx,lenx,T).*sin(2*pi*freqx*(0:T-1)');
   	
    y = ytrue + sqrt(vary)*randn(T,1);
    M = 200;
    params = log([varx;lenx;freqx;vary]);
    [mf,vf] = denoiseVFESMFFT(y,params,M,'y',1);
    figure(1)
    plot(1:T,y,'-b',1:T,ytrue,'-r',1:T,mf,'-g',...
         1:T,mf+2*sqrt(vf),'--g',1:T,mf-2*sqrt(vf),'--g');
    legend('noisy','true','predicted');
end



function test_notraining_y_2
	T = 1000;
	varx = 0.1;
    lenx = 400;
    vary = 0.1;
    freqx = 0.1;
    ytrue = sampleGPSM(varx,lenx,freqx,T);
    y = ytrue + sqrt(vary)*randn(T,1);
    M = 200;
    params = log([0.101117;258;0.01;0.08]);
    [mf,vf] = denoiseVFESMFFT(y,params,M,'y',1);
    figure(2)
    plot(1:T,y,'-b',1:T,ytrue,'-r',1:T,mf,'-g',...
         1:T,mf+2*sqrt(vf),'--g',1:T,mf-2*sqrt(vf),'--g');
    legend('noisy','true','predicted');
end

function test_notraining_f
	T = 1000;
	varx = 1;
    lenx = 20;
    vary = 0.1;
    freqx = 0.1;
    a1 = sampleGPSE(varx,lenx,T);
    b1 = sampleGPSE(varx,lenx,T);
    ytrue = a1.*cos(2*pi*freqx*(0:T-1)') ...
    	+ b1.*sin(2*pi*freqx*(0:T-1)');
    y = ytrue + sqrt(vary)*randn(T,1);
    M = 200;
    params = log([varx;lenx;freqx;vary]);
    [mf,vf] = denoiseVFESMFFT(y,params,M,'f',1);
    figure(3)
    plot(1:T,y,'-b',1:T,ytrue,'-r',1:T,a1,'-m',1:T,mf(:,1),'-g',...
         1:T,mf(:,1)+2*sqrt(vf(:,1)),'--g',1:T,mf(:,1)-2*sqrt(vf(:,1)),'--g');
    legend('noisy','true','a true','a predicted');
    figure(4)
    plot(1:T,y,'-b',1:T,ytrue,'-r',1:T,b1,'-m',1:T,mf(:,2),'-g',...
         1:T,mf(:,2)+2*sqrt(vf(:,2)),'--g',1:T,mf(:,2)-2*sqrt(vf(:,2)),'--g');
    legend('noisy','true','b true','b predicted');
end

function test_notraining_u
	T = 1000;
	varx = 1;
    lenx = 20;
    vary = 0.1;
    freqx = 0.1;
    a1 = sampleGPSE(varx,lenx,T);
    b1 = sampleGPSE(varx,lenx,T);
    ytrue = a1.*cos(2*pi*freqx*(0:T-1)') ...
    	+ b1.*sin(2*pi*freqx*(0:T-1)');
    y = ytrue + sqrt(vary)*randn(T,1);
    M = 200;
    params = log([varx;lenx;freqx;vary]);
    [mf,vf] = denoiseVFESMFFT(y,params,M,'f',1);
    [mu,vu] = denoiseVFESMFFT(y,params,M,'u',1);
    xu = linspace(1,T,M);
    figure(5)
    plot(1:T,a1,'-m',1:T,mf(:,1),'-g',xu,mu(:,1),'-r',...
         1:T,mf(:,1)+2*sqrt(vf(:,1)),'--g',1:T,mf(:,1)-2*sqrt(vf(:,1)),'--g');
         %xu,mu(:,1)+2*sqrt(vu(:,1)),'--g',xu,mu(:,1)-2*sqrt(vu(:,1)),'--g');
    legend('a true','a predicted','a u');
    figure(6)
    plot(1:T,b1,'-m',1:T,mf(:,2),'-g',xu,mu(:,2),'-r',...
         1:T,mf(:,2)+2*sqrt(vf(:,2)),'--g',1:T,mf(:,2)-2*sqrt(vf(:,2)),'--g');
         %xu,mu(:,2)+2*sqrt(vu(:,2)),'--g',xu,mu(:,2)-2*sqrt(vu(:,2)),'--g');
    legend('b true','b predicted','b u');
end
