function test_suite = test_trainVFESMFFT
% written by Richard Turner
% modified by Thang Bui
% initTestSuite
test_train1
end

function test_train1

close all;

dispFigs=1;
randn('state',1)

notrials = 1;
%L = 10;
%lenxs = logspace(log10(10),log10(500),L);
L = 1;
lenxs = logspace(log10(50),log10(50),L);
%freqxs = [0.01 0.05 0.1 0.2 0.25 0.3 0.4];
freqxsUsed = zeros(notrials,L);
varxsUsed = zeros(notrials,L);
lenxsEst = zeros(notrials,L);
varxsEst = zeros(notrials,L);
freqxsEst = zeros(notrials,L);
varysEst = zeros(notrials,L);
vary = 0.2;

T = 1000;
freqxs = 0.1;
setup.numIts = 500;
setup.progress_chunk = setup.numIts;

M = 100;

for n = 1:notrials
    for l = 1:L
        fprintf('Trial %d/%d, l %d/%d\n',n,notrials,l,L);
        lenx = lenxs(l);
        %varx = 0.2+exp(randn);
        varx = 1;
        freqx = freqxs(randi(length(freqxs),1));
        freqxsUsed(n,l) = freqx;
        varxsUsed(n,l) = varx;
        
        %ytrue = sampleGPSM(varx,lenx,freqx,T);
        ytrue = sampleGPSE(varx,lenx,T).*cos(2*pi*freqx*(0:T-1)') + sampleGPSE(varx,lenx,T).*sin(2*pi*freqx*(0:T-1)');

        y = ytrue + sqrt(vary)*randn(T,1);
        
        params = [log(varx),log(lenx),log(freqx),log(vary)]';
        [params,info] = trainVFESMFFT(y,M,setup,params);

        lenxsEst(n,l) = exp(params(2));
        varxsEst(n,l) = exp(params(1));
        freqxsEst(n,l) = exp(params(3));
        varysEst(n,l) = exp(params(end));
    end
end

figure, hold on;
plot(lenxs,lenxs,'-r')
for l = 1:L
    plot(lenxs(l),lenxsEst(:,l),'.k');
end
plot(lenxs,mean(lenxsEst,1),'-b');
if notrials > 1
    plot(lenxs,mean(lenxsEst,1)+sqrt(var(lenxsEst,1)),'-b');
    plot(lenxs,mean(lenxsEst,1)-sqrt(var(lenxsEst,1)),'-b');
end
xlabel('lenx used'), ylabel('lenx est')


figure, hold on;
scatter(freqxsUsed(:),freqxsEst(:));
xlabel('freqx used'), ylabel('freqx est')

figure, hold on;
scatter(varxsUsed(:),varxsEst(:));
xlabel('varx used'), ylabel('varx est')


%tol = 1e-1;
%assertVectorsAlmostEqual(autoCor,autoCorEmp','absolute',tol,0)

dl = abs(ones(notrials,1)*lenxs-lenxsEst)./(ones(notrials,1)*lenxs);
fprintf('lenx - average percentage error %f %% \n',100*mean(dl(:)));
dl = abs(varxsUsed-varxsEst)./varxsUsed;
fprintf('varx - average percentage error %f %% \n',100*mean(dl(:)));
dl = abs(freqxsUsed-freqxsEst)./freqxsUsed;
fprintf('freqx - average percentage error %f %% \n',100*mean(dl(:)));


[mf,vf] = denoiseVFESMFFT(y,params,M,'y',1);
figure, hold on,
plot(1:T,y,'-','Color',[0.85 0.85 0.85]), 
plot(1:T,ytrue,'-b',1:T,mf,'-r',1:T,mf-2*sqrt(vf),'-r',1:T,mf+2*sqrt(vf),'--r')
legend('noisy','true','predicted');

fprintf('lenx %f, lenxEst %f\n', lenx, exp(params(2)))
fprintf('varx %f, varxEst %f\n', varx, exp(params(1)))
fprintf('freqx %f, freqxEst %f\n', freqx, exp(params(3)))
fprintf('vary %f, varyEst %f\n', vary, exp(params(end)))
keyboard
end