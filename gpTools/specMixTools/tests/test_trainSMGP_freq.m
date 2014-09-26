function test_suite = test_trainSMGP_freq
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
L = 30;
lenxs = logspace(log10(10),log10(1000),L);
%freqxs = [0.01 0.05 0.1 0.2 0.25 0.3 0.4];
freqxs = [0];
freqxsUsed = zeros(notrials,L);
varxsUsed = zeros(notrials,L);
lenxsEst = zeros(notrials,L);
varxsEst = zeros(notrials,L);
freqxsEst = zeros(notrials,L);
varysEst = zeros(notrials,L);

T = 1024*8;
setup.numIts = 500;
setup.progress_chunk = setup.numIts;


for n = 1:notrials
    for l = 1:L
        fprintf('Trial %d/%d, l %d/%d\n',n,notrials,l,L);
        lenx = lenxs(l);
        %varx = 0.2+exp(randn);
        varx = 0.2;
        freqx = freqxs(randi(length(freqxs),1));
        freqxsUsed(n,l) = freqx;
        varxsUsed(n,l) = varx;
        
        ytrue = sampleGPSM(varx,lenx,freqx,T);
        y = ytrue + 0.3*randn(T,1);
        
        %[varxEst,lenxEst,freqxEst,varyEst,info] = trainSMGP_freq(y,1,setup,log(varx),log(lenx),log(freqx),log(0.1));
        [varxEst,lenxEst,freqxEst,varyEst,info] = trainSMGP_freq(y,1,setup);
        lenxsEst(n,l) = lenxEst;
        varxsEst(n,l) = varxEst;
        freqxsEst(n,l) = freqxEst;
        varysEst(n,l) = varyEst;
    end
end

figure, hold on;
plot(lenxs,lenxs,'-r')
for l = 1:L
    plot(lenxs(l),lenxsEst(:,l),'.k');
end
plot(lenxs,mean(lenxsEst,1),'-b');
plot(lenxs,mean(lenxsEst,1)+sqrt(var(lenxsEst,1)),'-b');
plot(lenxs,mean(lenxsEst,1)-sqrt(var(lenxsEst,1)),'-b');
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

mf = denoiseSMGP_freq(varxEst,lenxEst,freqxEst,varyEst,y);
figure, plot(1:T,ytrue,'-b',1:T,mf,'-r')
legend('true','predicted');

fprintf('lenx %f, lenxEst %f\n', lenx, lenxEst)
fprintf('varx %f, varxEst %f\n', varx, varxEst)
fprintf('freqx %f, freqxEst %f\n', freqx, freqxEst)
keyboard
end