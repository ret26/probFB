close all;
clear all;
clc;

noiters = 200;
noInducing = 10;
%% load data
data = load('../../tsgp/datasets/audio/audio_subband.mat');
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


%% GP-subset of data
t = (1:T)';
ttrain = t(randperm(T,800));
ytrain = y(ttrain);
covfunc = {@covSEiso};
hypsubset.cov = log([50,0.5]);
hypsubset.lik = log(0.1);
hypsubset = minimize(hypsubset,'gp',noiters,'infExact',[],covfunc,@likGauss,ttrain,ytrain);
[mysubset,vysubset,mffull,vffull] = gp(hypsubset,'infExact',[],covfunc,[],ttrain,ytrain,t);
figure(1), hold on;
h1 = boundedline(t,mysubset,2*sqrt(vysubset),'b','alpha');
plot(t,yOri,'k-')
plot(t,y,'color',[0.8 0.8 0.8]);
title('GP - subset of data')
hold off;

%% GP-FFT % TODO
%{
options.numIts = 200;
options.progress_chunk = 100;
[lenxEst,varxEst,varyEst] = trainSEGPwithNoise(y,options,50,1,0.1);
[myfft,vyfft] = predictSEGP(y,lenxEst,varxEst,varyEst);
figure(1), hold on;
plot(t,yOri,'k-')
plot(t,y,'color',[0.8 0.8 0.8]);
h1 = boundedline(t,mysubset,2*sqrt(vysubset),'b','alpha');
title('GP - subset of data')
hold off;
%}

%% VFE-SE
t = (1:T)';
D = 1;
nu = noInducing;
hypvfe.cov = log([50;1]);
hypvfe.lik = log(0.1);
hypvfe.Xu = t(randperm(T,nu),:);
covfunc = {@covSEiso};
theta = vfeTrain([hypvfe.cov;hypvfe.lik;hypvfe.Xu(:)],covfunc,t,y,noiters);
hypvfe.cov = theta(1:D+1);
hypvfe.lik = theta(D+2);
hypvfe.Xu = reshape(theta(D+3:end),nu,D);
disp('predicting ...');
[myvfe,vyvfe] = vfePredict(covfunc,hypvfe,t,y,t);
vyvfe = vyvfe + exp(2*hypvfe.lik);
figure(2), hold on;
h3 = boundedline(t,myvfe,2*sqrt(vyvfe),'b','alpha');
plot(t,yOri,'k-')
plot(t,y,'color',[0.8 0.8 0.8]);
title('VFE - optimised u')
hold off;

%% VFE-FFT
setup.numIts = noiters;
setup.progress_chunk = noiters;
kernel.name = 'SE';
params = log([1;20;0.1]);
M = 2000;
[paramsOpt,info] = trainVFEFFT(y,kernel,M,setup,params);
oriSpec = getVFESpecHelper(kernel,paramsOpt(1:end-1),T);
[mfvfefft,vfvfefft] = denoiseVFEFFT(y,oriSpec,exp(paramsOpt(end)),M);
vyvfefft = vfvfefft + exp(paramsOpt(end));
figure(3), hold on;
h3 = boundedline(t,mfvfefft,2*sqrt(vyvfefft),'b','alpha');
plot(t,yOri,'k-')
plot(t,y,'color',[0.8 0.8 0.8]);
title('VFE FFT - regularly sampled u')
hold off;
%% AR(2) % TODO

%% Chain-SE
tau1 = 50;
tau2 = 5;

K = floor(T/tau1);
Xtrain = (1:tau1*K)';
Ytrain = y(Xtrain);
Ytrue = yOri(Xtrain);
missingStack = zeros([K,tau1])==1;

covfunc = {@covSEiso};
params = zeros(3,1);
params(1) = 50; % lengthscale
params(2) = sqrt(var(Ytrain)); % signal variance
params(3) = 1/2*sqrt(var(Ytrain)); % noise variance
theta_init = log(params);

[theta_end,nlml] = trainSE(theta_init,covfunc,...
    Xtrain,Ytrain,tau1,tau2,missingStack,noiters);

[mfchain,vfchain] = predictSE(theta_end,covfunc,...
    Xtrain,Ytrain,tau1,tau2,missingStack);
vychain = vfchain + exp(theta_end(end));

figure(4), hold on;
plot(t,yOri,'k-')
plot(t,y,'color',[0.8 0.8 0.8]);
h3 = boundedline(Xtrain,mfchain,2*sqrt(vychain),'b','alpha');
title('Chain - tau1=50, tau2=10')
hold off;

%% Chain-SE Local
tau1 = 50;
tau2 = 5;

K = floor(T/tau1);
Xtrain = (1:tau1*K)';
Ytrain = y(Xtrain);
Ytrue = yOri(Xtrain);
missingStack = zeros([K,tau1])==1;

covfunc = {@covSEiso};
params = zeros(3,1);
params(1) = 50; % lengthscale
params(2) = sqrt(var(Ytrain)); % signal variance
params(3) = 1/2*sqrt(var(Ytrain)); % noise variance
theta_init = log(params);

[theta_end,nlml] = trainSELocal(theta_init,covfunc,...
    Xtrain,Ytrain,tau1,tau2,missingStack,noiters);

[mfchain,vfchain] = predictSELocal(theta_end,covfunc,...
    Xtrain,Ytrain,tau1,tau2,missingStack);
vychain = vfchain + exp(theta_end(end));

figure(5), hold on;
h3 = boundedline(Xtrain,mfchain,2*sqrt(vychain),'b','alpha');
title('Local - tau1=50, tau2=10')
plot(t,yOri,'k-')
plot(t,y,'color',[0.8 0.8 0.8]);
hold off;

%% FITC-SE
params = log([50;1;0.1]);
covfunc = {@covSEiso};
M = 100;
[params_end] = trainSEFITC(params,covfunc,t,y,M,noiters);
[mffitc,vffitc] = predictSEFITC(params_end,covfunc,t,y,t,M);
vyfitc = vffitc + exp(2*params_end.lik);
figure(6), hold on;
h3 = boundedline(t,mffitc,2*sqrt(vyfitc),'b','alpha');
plot(t,yOri,'k-')
plot(t,y,'color',[0.8 0.8 0.8]);
title('FITC - linearly spaced u')
hold off;

%% SDE-Kalman -- be careful, can be very slow
Nsubset = 2000;
tsubset = (1:Nsubset)';
ysubset = y(tsubset);
params = log([50;1;0.1]);
approxDeg = 2;
[params_end] = trainSSM_SE(params,tsubset,ysubset,approxDeg,noiters);
[mysde,vfsde] = predictSSM_SE(params_end,tsubset,ysubset,tsubset,approxDeg);
vysde = vfsde + exp(2*params_end(end));

figure(7), hold on;
h3 = boundedline(tsubset,mysde,2*sqrt(vysde),'b','alpha');
plot(tsubset,yOri(tsubset),'k-')
plot(tsubset,ysubset,'color',[0.8 0.8 0.8]);
title('SDE - tau1=50, tau2=10')
hold off;
%% Sparse spectrum
params = log([50;1;0.1]);
M = 100;
spectralPoints = randn(M,1);
params = [params;spectralPoints];
params_end = minimize(params,'ssgpr',-noiters,t,y);
[myssgp,vyssgp] = ssgpr(params_end,t,y,t);
figure(8), hold on;
h3 = boundedline(t,myssgp,2*sqrt(vyssgp),'b','alpha');
plot(t,yOri,'k-')
plot(t,y,'color',[0.8 0.8 0.8]);
title('Sparse spectrum')
hold off;