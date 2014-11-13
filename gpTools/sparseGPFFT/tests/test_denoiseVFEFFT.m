function test_denoiseVFEFFT
   % initTestSuite;
    test_notraining;
end

function test_notraining
    T = 1000;
    varx = 1;
    lenx = 20;
    vary = 0.1;
    ytrue = sampleGPSE(varx,lenx,T);
    y = ytrue + sqrt(vary)*randn(T,1);
    M = 100;
    oriSpec = getGPSESpec(lenx,T);
    [mf,vf] = denoiseVFEFFT(y,oriSpec,vary,M);
    figure
    plot(1:T,y,'-b',1:T,ytrue,'-r',1:T,mf,'-g',...
         1:T,mf+2*sqrt(vf),'--g',1:T,mf-2*sqrt(vf),'--g');
    legend('noisy','true','predicted');
end

function test_withtraining
end

