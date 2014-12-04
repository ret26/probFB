close all;
clear all;
clc
%resPath = './results/working/16_7/';
%mkdir(resPath);
exnFile = fopen('exceptions.txt','a'); % exception file

%% load data, setting
data = load('audio_subband.mat');
y = data.y;
T = length(y);
y = real(hilbert(y).*exp(-1i*2*pi*data.mu*[1:T]'));
yOri = y;
y = y+randn(T,1)/50; % add some noise to the original data

% only use first few samples for training and testing for now
y = y(8e3+(1:2e5));
y = resample(y,1,4);
yOri = yOri(8e3+(1:2e5));
yOri = resample(yOri,1,4);
%figure, plot(y)
T = length(y);

%% setting missing blocks
s = 80;
mInd = [400 1500 2400 3700 4800 9000 9400 ...
    1.16e4 1.3e4 1.7e4 1.8e4 2.02e4 2.14e4 ...
    2.85e4 2.9e4 3e4 3.55e4 3.8e4 3.96e4 ...
    4.52e4 4.6e4 4.72e4 4.78e4 4.83e4 4.87e4];
missingInd = zeros(T,1);
for m = 1:length(mInd);
    j = (mInd(m))+(1:randi(1)*s);
    j(j>T) = [];
    missingInd(j) = 1;
end
missingInd = missingInd==1;

%%
Xtrain = 1:T;
YtestOri = y; YtestOri(missingInd) = 0;
YtrainOri = y;
missingTrain = missingInd;
YtrainOri(missingTrain) = 0;

%% training and prediction

% parameters
tau2vec = [2 5 10 20];
ratevec = [1/4 1/5 1/10 1/20 1/25];
MM = [16 32 64 128 256 512 1024 1500];
MM = [1024];
approxDeg = 1:10;

l1 = length(tau2vec);
l2 = length(ratevec);
lm = length(MM);
ld = length(approxDeg);
noTrials = 1;
noEvals = 50;

methods = {'chain','local','fitc','var','ssgp','sdeKalman'};
methods = {'var'};
for k = 1:noTrials
    for m = 1:length(methods)
        method = methods{m};
        res = struct;
        if strcmpi(method,'chain')
            smse = zeros(l1,l2); % reconstruction error
            msll = zeros(l1,l2);
            trainTime = zeros(l1,l2); % training time
            testTime = zeros(l1,l2); % test time
            hypers = cell(l1,l2); % hypers
            for i = 1:l1
                for j = 1:l2
                    disp(['running trial ' num2str(k) '/' num2str(noTrials) ', '...
                        'method ' method ', ' ...
                        'i = ' num2str(i) '/' num2str(l1) ', '...
                        'j = ' num2str(j) '/' num2str(l2)]);
                    tau2 = tau2vec(i);
                    rate = ratevec(j);
                    tau1 = tau2/rate;
                    
                    K = floor(T/tau1);
                    Xtrain = 1:tau1*K;
                    Ytrain = YtrainOri(Xtrain);
                    Ytrue = yOri(Xtrain);
                    Ytest = YtestOri(Xtrain);
                    missing = missingTrain(Xtrain);
                    missingStack = reshape(missing,[tau1,K])';
                    
                    covfunc = {@covSEiso};
                    params = zeros(3,1);
                    params(1) = 50; % lengthscale
                    params(2) = sqrt(var(Ytrain)); % signal variance
                    params(3) = 1/2*sqrt(var(Ytrain)); % noise variance
                    theta_init = log(params);
                    try
                        tic
                        [theta_end,nlml] = trainSE(theta_init,covfunc,...
                            Xtrain',Ytrain,tau1,tau2,missingStack,noEvals);
                        trainTime(i,j) = toc;
                        tic
                        [fest,vest] = predictSE(theta_end,covfunc,...
                            Xtrain',Ytest,tau1,tau2,missingStack);
                        testTime(i,j) = toc;
                        
                        % data reconstruction loss
                        ytrue = yOri(missing);
                        yreco = fest(missing);
                        vreco = vest(missing);
                        smse(i,j) = smsError(ytrue,yreco);
                        meanTrain =  mean(Ytrain);
                        varTrain = var(Ytrain);
                        msll(i,j) = mslLoss(ytrue,yreco,vreco,meanTrain,varTrain);
                        hypers{i,j} = theta_end;
                        
                    catch exception
                        disp(exception);
                        trainTime(i,j) = NaN;
                        testTime(i,j) = NaN;
                        smse(i,j) = NaN;
                        msll(i,j) = NaN;
                        msg = [datestr(now) getReport(exception,'extended') '\n'];
                        fprintf(exnFile,msg);
                    end
                    
                end
            end
            res.msll = msll;
            res.smse = smse;
            res.tau2vec = tau2vec;
            res.ratevec = ratevec;
            res.trainTime = trainTime;
            res.testTime = testTime;
            res.hypers = hypers;
            res.Ytrain = YtrainOri;
        elseif strcmpi(method,'local')
            smse = zeros(l1,l2); % reconstruction error
            msll = zeros(l1,l2);
            trainTime = zeros(l1,l2); % training time
            testTime = zeros(l1,l2); % test time
            hypers = cell(l1,l2); % hypers
            for i = 1:l1
                for j = 1:l2
                    disp(['running trial ' num2str(k) '/' num2str(noTrials) ', '...
                        'method ' method ', ' ...
                        'i = ' num2str(i) '/' num2str(l1) ', '...
                        'j = ' num2str(j) '/' num2str(l2)]);
                    tau2 = tau2vec(i);
                    rate = ratevec(j);
                    tau1 = tau2/rate;
                    
                    K = floor(T/tau1);
                    Xtrain = 1:tau1*K;
                    Ytrain = YtrainOri(Xtrain);
                    Ytrue = yOri(Xtrain);
                    Ytest = YtestOri(Xtrain);
                    missing = missingTrain(Xtrain);
                    missingStack = reshape(missing,[tau1,K])';
                    
                    covfunc = {@covSEiso};
                    params = zeros(3,1);
                    params(1) = 50; % lengthscale
                    params(2) = sqrt(var(Ytrain)); % signal variance
                    params(3) = 1/2*sqrt(var(Ytrain)); % noise variance
                    theta_init = log(params);
                    try
                        tic
                        [theta_end,nlml] = trainSELocal(theta_init,covfunc,...
                            Xtrain',Ytrain,tau1,tau2,missingStack,noEvals);
                        trainTime(i,j) = toc;
                        tic
                        [fest,vest] = predictSELocal(theta_end,covfunc,...
                            Xtrain',Ytest,tau1,tau2,missingStack);
                        testTime(i,j) = toc;
                        
                        % data reconstruction loss
                        ytrue = yOri(missing);
                        yreco = fest(missing);
                        vreco = vest(missing);
                        smse(i,j) = smsError(ytrue,yreco);
                        meanTrain =  mean(Ytrain);
                        varTrain = var(Ytrain);
                        msll(i,j) = mslLoss(ytrue,yreco,vreco,meanTrain,varTrain);
                        hypers{i,j} = theta_end;
                        
                    catch exception
                        disp(exception);
                        trainTime(i,j) = NaN;
                        testTime(i,j) = NaN;
                        smse(i,j) = NaN;
                        msll(i,j) = NaN;
                        msg = [datestr(now) getReport(exception,'extended') '\n'];
                        fprintf(exnFile,msg);
                    end
                    
                end
            end
            res.msll = msll;
            res.smse = smse;
            res.tau2vec = tau2vec;
            res.ratevec = ratevec;
            res.trainTime = trainTime;
            res.testTime = testTime;
            res.hypers = hypers;
            res.Ytrain = YtrainOri;
        elseif strcmpi(method,'fitc')
            smse = zeros(lm,1); % reconstruction error
            msll = zeros(lm,1);
            trainTime = zeros(lm,1); % training time
            testTime = zeros(lm,1); % test time
            hypers = cell(lm,1); % hypers
            XtrainFITC = Xtrain;
            XtrainFITC(missingTrain) = [];
            YtrainFITC = YtrainOri;
            YtrainFITC(missingTrain) = [];
            XtestFITC = Xtrain(missingTrain);
            Ytrue = yOri(missingTrain);
            for i = 1:lm
                disp(['running trial ' num2str(k) '/' num2str(noTrials) ', '...
                    'method ' method ', ' ...
                    'i = ' num2str(i) '/' num2str(lm)]);
                covfunc = {@covSEiso};
                params = zeros(3,1);
                params(1) = 50; % lengthscale
                params(2) = sqrt(var(YtrainFITC)); % signal variance
                params(3) = 1/2*sqrt(var(YtrainFITC)); % noise variance
                
                theta_init = log(params);
                try
                    tic
                    [theta_end,nlml] = trainSEFITC(theta_init,covfunc,...
                        XtrainFITC',YtrainFITC,MM(i),noEvals);
                    trainTime(i) = toc;
                    tic
                    [fest,vest] = predictSECosFITC(theta_end,covfunc,XtrainFITC',...
                        YtrainFITC,XtestFITC',MM(i));
                    testTime(i) = toc;
                    % data reconstruction loss
                    smse(i) = smsError(Ytrue,fest);
                    meanTrain =  mean(YtrainFITC);
                    varTrain = var(YtrainFITC);
                    msll(i) = mslLoss(Ytrue,fest,vest,meanTrain,varTrain);
                    hypers{i} = theta_end;
                catch exception
                    disp(exception);
                    trainTime(i) = NaN;
                    testTime(i) = NaN;
                    smse(i) = NaN;
                    msll(i) = NaN;
                    msg = [datestr(now) getReport(exception,'extended') '\n'];
                    fprintf(exnFile,msg);
                end
                
            end
            res.msll = msll;
            res.smse = smse;
            res.M = MM;
            res.trainTime = trainTime;
            res.testTime = testTime;
            res.hypers = hypers;
            res.Ytrain = YtrainOri;
        elseif strcmpi(method,'var')
            smse = zeros(lm,1); % reconstruction error
            msll = zeros(lm,1);
            trainTime = zeros(lm,1); % training time
            testTime = zeros(lm,1); % test time
            hypers = cell(lm,1); % hypers
            XtrainFITC = Xtrain;
            XtrainFITC(missingTrain) = [];
            YtrainFITC = YtrainOri;
            YtrainFITC(missingTrain) = [];
            XtestFITC = Xtrain(missingTrain);
            Ytrue = yOri(missingTrain);
            for i = 1:lm
                disp(['running trial ' num2str(k) '/' num2str(noTrials) ', '...
                    'method ' method ', ' ...
                    'i = ' num2str(i) '/' num2str(lm)]);
                covfunc = {@covSEiso};
                params = zeros(3,1);
                params(1) = 50; % lengthscale
                params(2) = sqrt(var(YtrainFITC)); % signal variance
                params(3) = 1/2*sqrt(var(YtrainFITC)); % noise variance
                theta_init = log(params);
                try
                    [trainingTime,testTimei,hyp,nlml,fest,vest] = ...
                        evalVarSE(XtrainFITC',YtrainFITC,XtestFITC',MM(i),theta_init,covfunc,noEvals);
                    
                    trainTime(i) = trainingTime;
                    
                    testTime(i) = testTimei;
                    % data reconstruction loss
                    smse(i) = smsError(Ytrue,fest);
                    meanTrain =  mean(YtrainFITC);
                    varTrain = var(YtrainFITC);
                    msll(i) = mslLoss(Ytrue,fest,vest,meanTrain,varTrain);
                    hypers{i} = hyp;
                catch exception
                    disp(exception);
                    trainTime(i) = NaN;
                    testTime(i) = NaN;
                    smse(i) = NaN;
                    msll(i) = NaN;
                    msg = [datestr(now) getReport(exception,'extended') '\n'];
                    fprintf(exnFile,msg);
                end
                
            end
            res.msll = msll;
            res.smse = smse;
            res.M = MM;
            res.trainTime = trainTime;
            res.testTime = testTime;
            res.hypers = hypers;
            res.Ytrain = YtrainOri;
        elseif strcmpi(method,'ssgp')
            smse = zeros(lm,1); % reconstruction error
            msll = zeros(lm,1);
            trainTime = zeros(lm,1); % training time
            testTime = zeros(lm,1); % test time
            hypers = cell(lm,1); % hypers
            XtrainFITC = Xtrain;
            XtrainFITC(missingTrain) = [];
            YtrainFITC = YtrainOri;
            YtrainFITC(missingTrain) = [];
            XtestFITC = Xtrain(missingTrain);
            Ytrue = yOri(missingTrain);
            for i = 1:lm
                disp(['running trial ' num2str(k) '/' num2str(noTrials) ', '...
                    'method ' method ', ' ...
                    'i = ' num2str(i) '/' num2str(lm)]);
                covfunc = {@covSEiso};
                params = zeros(3,1);
                params(1) = 50; % lengthscale
                params(2) = sqrt(var(YtrainFITC)); % signal variance
                params(3) = 1/2*sqrt(var(YtrainFITC)); % noise variance
                theta_init = log(params);
                try
                    [trainingTime,testTimei,hyp,nlml,fest,vest] = ...
                        evalSSGPSE(XtrainFITC',YtrainFITC,XtestFITC',MM(i),theta_init,noEvals);
                    
                    trainTime(i) = trainingTime;
                    
                    testTime(i) = testTimei;
                    % data reconstruction loss
                    smse(i) = smsError(Ytrue,fest);
                    meanTrain =  mean(YtrainFITC);
                    varTrain = var(YtrainFITC);
                    msll(i) = mslLoss(Ytrue,fest,vest,meanTrain,varTrain);
                    
                    % envelop reconstruction loss
                    env1 = abs(hilbert(Ytrue));
                    env2 = abs(hilbert(fest));
                    hypers{i} = hyp;
                catch exception
                    disp(exception);
                    trainTime(i) = NaN;
                    testTime(i) = NaN;
                    smse(i) = NaN;
                    msll(i) = NaN;
                    msg = [datestr(now) getReport(exception,'extended') '\n'];
                    fprintf(exnFile,msg);
                end
                
            end
            res.msll = msll;
            res.smse = smse;
            res.M = MM;
            res.trainTime = trainTime;
            res.testTime = testTime;
            res.hypers = hypers;
            res.Ytrain = YtrainOri;
            
        elseif strcmpi(method,'sdeKalman')
            smse = zeros(ld,1); % reconstruction error
            msll = zeros(ld,1);
            trainTime = zeros(ld,1); % training time
            testTime = zeros(ld,1); % test time
            hypers = cell(ld,1); % hypers
            Nd = 2e4;
            Xtrain1 = Xtrain(1:Nd); % only consider 10k samples
            missingTrain1 = missingTrain(1:Nd);
            YtrainOri1 = YtrainOri(1:Nd);
            yOri1 = yOri(1:Nd);
            
            XtrainSDE = Xtrain1;
            XtrainSDE(missingTrain1) = [];
            YtrainSDE = YtrainOri1;
            YtrainSDE(missingTrain1) = [];
            XtestSDE = Xtrain1(missingTrain1);
            Ytrue = yOri1(missingTrain1);
            for i = 1:ld
                disp(['running trial ' num2str(k) '/' num2str(noTrials) ', '...
                    'method ' method ', ' ...
                    'i = ' num2str(i) '/' num2str(ld)]);
                covfunc = {@covSEiso};
                params = zeros(3,1);
                params(1) = 50; % lengthscale
                params(2) = sqrt(var(YtrainSDE)); % signal variance
                params(3) = 1/2*sqrt(var(YtrainSDE)); % noise variance
                theta_init = log(params);
                try
                    [trainingTime,testTimei,hyp,nlml,fest,vest] = ...
                        evalSSMSE(XtrainSDE',YtrainSDE,XtestSDE',approxDeg(i),theta_init,noEvals);
                    
                    trainTime(i) = trainingTime;
                    testTime(i) = testTimei;
                    % data reconstruction loss
                    smse(i) = smsError(Ytrue,fest);
                    meanTrain =  mean(YtrainSDE);
                    varTrain = var(YtrainSDE);
                    msll(i) = mslLoss(Ytrue,fest,vest,meanTrain,varTrain);
                    hypers{i} = hyp;
                catch exception
                    disp(exception);
                    trainTime(i) = NaN;
                    testTime(i) = NaN;
                    smse(i) = NaN;
                    msll(i) = NaN;
                    msg = [datestr(now) getReport(exception,'extended') '\n'];
                    fprintf(exnFile,msg);
                end
                
            end
            res.msll = msll;
            res.smse = smse;
            res.approxDeg = approxDeg;
            res.trainTime = trainTime;
            res.testTime = testTime;
            res.hypers = hypers;
            res.Ytrain = YtrainOri;
        end
        keyboard
        %save([resPath 'result_subbandSE_' method '_trial_' num2str(k) '.mat'],'res');
    end
end
