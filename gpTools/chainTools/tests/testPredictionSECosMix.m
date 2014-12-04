% test prediction cosse
%resPath = '../results/working/17_6/';
resPath = '../probLDS/results/working/17_6/';
data1 = load([resPath 'result_SECos_two_bands_chain_trial_1.mat']);
res1 = data1.res;
ind1 = 3; ind2 = 3;
%%
data = load('five_bands.mat');
mu = data.mu;
y = data.y1+data.y2;
y = y/std(y);

T = length(y);
Yori = y;
y = Yori+randn(T,1)/100; % add some noise to the original data
% only use first few samples for training and testing for now
%y = y(1:3.5e4);
y = resample(y,1,4);
%figure, plot(y)
T = length(y);



%% randomly sample the envelop and set the missing blocks - not really!
s = 150;
mInd = [9200 1.12e4 1.4e4 1.55e4 1.74e4 2.35e4 2.65e4 2.8e4];
missingInd = zeros(T,1);
for m = 1:length(mInd);
    j = (mInd(m))+(1:randi(1)*s);
    j(j>T) = [];
    missingInd(j) = 1;
end
missingInd = missingInd==1;

%%
% Xtrain = 1:T;
YtestOri = y; YtestOri(missingInd) = 0;
YtrainOri = y;
missingTrain = missingInd;
%missingTrain = zeros(size(missingInd));
YtrainOri(missingTrain) = 0;

%%
tau2vec = res1.tau2vec;
ratevec = res1.ratevec;
tau2 = tau2vec(ind1);
rate = ratevec(ind2);
tau1 = tau2/rate;
hyp1 = res1.hypers{ind1,ind2};

covfunc = {@covSEiso};
K = floor(T/tau1);
Xtrain = 1:tau1*K;
Ytrain = YtrainOri(Xtrain);
Ytrue = y(Xtrain);
Ytest = YtestOri(Xtrain);
missing = missingTrain(Xtrain);
missingStack = reshape(missing,[tau1,K])';
noComps = 2;
[fest1,vest1] = predictSECosMix(hyp1,noComps,covfunc,...
    Xtrain',Ytest,tau1,tau2,missingStack);
figure(1), 
plot(Xtrain,Ytest,'g',Xtrain,Ytrue,'b',Xtrain,fest1,'k')


legend('Missing','True','Predicted');

ytrue = y(missing);
yreco1 = fest1(missing);
smse1 = smsError(ytrue,yreco1);
[smse1]

% envelop reconstruction loss
envOri = abs(hilbert(y));
env1 = abs(hilbert(fest1));
etrue = envOri(missing);
ereco1 = env1(missing);
smse3 = smsError(etrue,ereco1);
[smse3]


XtrainFITC = Xtrain;
XtrainFITC(missingTrain) = [];
YtrainFITC = YtrainOri;
YtrainFITC(missingTrain) = [];
XtestFITC = Xtrain(missingTrain);
Ytrue = y(missingTrain);
covfunc = {@covProd,{@covSEisoU,@covCos}};


figure, plot(Xtrain,envOri(1:length(Xtrain)),'r',Xtrain,env1,'g')
legend('True','Chain');

