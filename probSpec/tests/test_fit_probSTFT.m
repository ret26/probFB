function test_suite = test_fit_probSTFT
  initTestSuite;

% Tests:
%
% function [Lam,Var,Info] = fitAR2FB(y,D,varargin)
% 
  
function test_pSTFT

% Tests the function runs on data generated from a probabilistic spectrum

dispFigs = 1;
T = 10000;

% D = 3;
% CF = [1/32;1/20;1/10];
% DF = [1/128;1/50;1/100];
% varMa = [1;4;2];

D = 4;
fmaxTr = [1/32;1/20;1/10;1/15];
dfTr = [1/128;1/100;1/100;1/500];
varMaTr = [1;4;2;2];

vary = 0;

[omTr,lamxTr,varxTr] = freq2probSpec(fmaxTr,dfTr,varMaTr);

[y,X] = samplePSTFT(lamxTr,varxTr,omTr,vary,T);
[yHO,X] = samplePSTFT(lamxTr,varxTr,omTr,vary,T);

opts.verbose = 1;
opts.yHO = yHO;

[varxEst,lamxEst,omEst,Info] = fit_probSTFT(y,D,opts);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
specSig = abs(fft(y)).^2;
freqSig = linspace(0,1/2,floor(T/2));

specTr = zeros(D,T);
specEst = zeros(D,T);

freqs = linspace(0,1/2,T);

specTr = get_pSTFT_spec(freqs,lamxTr,varxTr,omTr);
specEst = get_pSTFT_spec(freqs,lamxEst,varxEst,omEst);
  
cumsumSig = 2*cumsum(specSig(1:floor(T/2))/T);
cumsumTr = cumsum(sum(specTr,1));
cumsumEst = cumsum(sum(specEst,1));

if dispFigs==1

  figure
  subplot(2,1,1)
  hold on

  ah1=plot(freqSig,specSig(1:floor(T/2))/T,'-k');
  for d=1:D
    ah2=plot(freqs,specTr(d,:),'-g');
    ah3=plot(freqs,specEst(d,:),'-m');
  end
  legend([ah1,ah2,ah3],'data','true','estimated')

  subplot(2,1,2)
  hold on
  ah1=plot(freqSig,cumsumSig,'-k');
  ah2=plot(freqs,cumsumTr,'-g');
  ah3=plot(freqs,cumsumEst,'-m');
  legend([ah1,ah2,ah3],'data','true','estimated')
end

tol = 5e-1;
assertElementsAlmostEqual(100+cumsumEst(1:2:end)',100+cumsumSig,'relative',tol,0)


function testSpeech

% Tests the function runs on speech - run using:
% $ runtests test_fitAR2FB:testSpeech
dispFigs=1;

File = '74 - Sentences.wav';
RngLim = [1.2,4.2];

File = '87 - Wind and Leaves.wav';
RngLim = [1.2,4.2];

File = '80 - Fire.wav';
RngLim = [0.2,1.5];

File = '05 - Bubbling stream.wav';
RngLim = [1.2,4.2]+3;

File = '07 - Bullfrogs.wav';
RngLim = [1,8.5];

File = '06 - Budgerigars.wav';
RngLim = [13,17];

File = '24 - Fire Burning.wav';
RngLim = [10,15];

File = '81 - Wind.wav';
RngLim = [.1,3];

File = '82 - Rain.wav';
RngLim = [10,15];

soundPath = '~/Music/Various Artists/Best Ever Sound Effects - Vol.3 - Sounds Of Nature Sound Effects/';%soundPath = '~/probFB/demos/signals/';
[yOrig,FSamp] = wavread([soundPath,File]);
%plot([1:length(yOrig(:,1))]/FSamp,yOrig(:,1))

RngLim = round(RngLim*FSamp);

y = yOrig(RngLim(1):RngLim(2),1);
y = y - mean(y);
y = y/sqrt(var(y));
yHO = yOrig(RngLim(2)+(RngLim(1):RngLim(2)));
yHO = yHO/sqrt(var(yHO));

D = 15;

opts.verbose = 1;
opts.yHO = yHO;

[varxEst,lamxEst,omEst,Info] = fit_probSTFT(y,D,opts);

Info.likeHO(end)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
T = length(y);
specSig = abs(fft(y)).^2;
freqSig = linspace(0,1/2,floor(T/2));

freqs = linspace(0,1/2,T);

specEst = get_pSTFT_spec(freqs,lamxEst,varxEst,omEst);


cumsumSig = 2*cumsum(specSig(1:floor(T/2))/T);
cumsumEst = cumsum(sum(specEst,1));

if dispFigs==1

  figure
  subplot(2,1,1)
  hold on

  ah1=plot(freqSig,specSig(1:floor(T/2))/T,'-k');
  for d=1:D
    ah3=plot(freqs,specEst(d,:),'-m');
  end
  legend([ah1,ah3],'data','estimated')
  set(gca,'yscale','log')
  
  subplot(2,1,2)
  hold on
  ah1=plot(freqSig,cumsumSig,'-k');
  ah3=plot(freqs,cumsumEst,'-m');
  legend([ah1,ah3],'data','estimated')
end

tol = 5e-1;
assertElementsAlmostEqual(100+cumsumEst(1:2:end)',100+cumsumSig,'relative',tol,0)

function testHOLikes

% compares the HO likelihoods on a bunch of signals - not a test
% but used to tune algorithmic parameters

dispFigs=1;

RngLim = zeros(9,2);

File{1} = '74 - Sentences.wav';
RngLim(1,:) = [1.2,4.2];

File{2} = '87 - Wind and Leaves.wav';
RngLim(2,:) = [1.2,4.2];

File{3} = '80 - Fire.wav';
RngLim(3,:) = [0.2,1.5];

File{4} = '05 - Bubbling stream.wav';
RngLim(4,:) = [1.2,4.2]+3;

File{5} = '07 - Bullfrogs.wav';
RngLim(5,:) = [1,8.5];

File{6} = '06 - Budgerigars.wav';
RngLim(6,:) = [13,17];

File{7} = '24 - Fire Burning.wav';
RngLim(7,:) = [10,15];

File{8} = '81 - Wind.wav';
RngLim(8,:) = [.1,3];

File{9} = '82 - Rain.wav';
RngLim(9,:) = [10,15];

K = size(RngLim,1);

soundPath = '~/Music/Various Artists/Best Ever Sound Effects - Vol.3 - Sounds Of Nature Sound Effects/';%soundPath = '~/probFB/demos/signals/';

% Preallocation of train and test likelihoods
Likes1 = zeros(K,1); 
HOLikes1 = zeros(K,1);
Likes2 = zeros(K,1); 
HOLikes2 = zeros(K,1);
tim1 = zeros(K,1);
tim2 = zeros(K,1);


for k=1:K

  [yOrig,FSamp] = wavread([soundPath,File{k}]);
  %plot([1:length(yOrig(:,1))]/FSamp,yOrig(:,1))

  RngLimCur = round(RngLim(k,:)*FSamp);

  y = yOrig(RngLimCur(1):RngLimCur(2),1);
  y = y - mean(y);
  y = y/sqrt(var(y));
  yHO = yOrig(RngLimCur(2)+(RngLimCur(1):RngLimCur(2)));
  yHO = yHO/sqrt(var(yHO));

  D = 15;

  opts.verbose = 1;
  opts.yHO = yHO;

  opts.maxT = 1000;
  opts.numIts = 10;
  opts.numLevels = 60;
  opts.vary_an = logspace(log10(1e-6),log10(1e-10),opts.numLevels);
  opts.bet = 300;
  opts.reassign = 1;
  opts.reassignThresh = 8;
  
  tic;
  [varxEst1,lamxEst1,omEst1,Info1] = fit_probSTFT(y,D,opts);
  tim1(k)=toc;

  % [fmax1,df1,varMa1] = probSpec2freq(omEst1,lamxEst1, varxEst1);
  % [val,ind] = sort(fmax1);
  % [fmax1(ind),df1(ind)]*1000
  
  opts.maxT = 1000;
  opts.numIts = 10;
  opts.numLevels = 40;
  opts.vary_an = logspace(log10(1e-6),log10(1e-10),opts.numLevels);
  opts.bet = 100;
  opts.reassign = 1;
  opts.reassignThresh = 8;
  
  tic;
  [varxEst2,lamxEst2,omEst2,Info2] = fit_probSTFT(y,D,opts);
  tim2(k)=toc;

  % [fmax2,df2,varMa2] = probSpec2freq(omEst2,lamxEst2, varxEst2);
  % [val,ind] = sort(fmax2);
  % [fmax2(ind),df2(ind)]*1000

  
  HOLikes1(k) = Info1.likeHO(end);
  Likes1(k) = Info1.likeUnReg(end);
  HOLikes2(k) = Info2.likeHO(end);
  Likes2(k) = Info2.likeUnReg(end);
  
  [HOLikes1,HOLikes2]  
  [Likes1,Likes2]
 % lamxEst1
 % lamxEst2
end

tot_tim1 = round(sum(tim1));
tot_tim2 = round(sum(tim2));
disp(['time taken overall: ',num2str(tot_tim1+tot_tim2),'s'])
disp(['time taken for first condition: ',num2str(tot_tim1),'s'])
disp(['time taken for second condition: ',num2str(tot_tim2),'s'])


figure
subplot(2,1,1)
title('held out likelihoods')
hold on
% plot(HOLikes1,HOLikes2,'or','markerfacecolor',[1,0,0])
% minLik = min([HOLikes1;HOLikes2]);
% maxLik = max([HOLikes1;HOLikes2]);
% plot([minLik,maxLik],[minLik,maxLik],'-k')
% set(gca,'xlim',[minLik,maxLik],'ylim',[minLik,maxLik])

dHOLik = HOLikes1-HOLikes2;
plot([1:K],HOLikes1-HOLikes2,'or','markerfacecolor',[1,0,0])
%minLik = min([HOLikes1;HOLikes2]);
%maxLik = max([HOLikes1;HOLikes2]);
%plot([minLik,maxLik],[minLi,maxLik],'-k')
set(gca,'xlim',[0.5,K+1/2],'ylim',[min(dHOLik),max(dHOLik)])
xlabel('data set')
ylabel('improvement in HO likelihood/nats per data point')
title('greater than zero is better')
plot([1/2,K+1/2],[0,0],'-k')

subplot(2,1,2)
title('overfitting measure')
hold on

plot([1:K],HOLikes1-Likes1,'ob','markerfacecolor',[0,0,1],'markersize',7)
plot([1:K],HOLikes2-Likes2,'or','markerfacecolor',[1,0,0])
legend('condition 1','condition 2')
xlabel('data set')
ylabel('amount of over-fitting (training lik - test lik)')
set(gca,'xlim',[0.5,K+1/2])
plot([1/2,K+1/2],[0,0],'-k')


mn1=mean(HOLikes1);
sig1=sqrt(var(HOLikes1));
mn2=mean(HOLikes2);
sig2=sqrt(var(HOLikes2));

disp(['mean like cond 1 = ',num2str(mn1),' +/- ',num2str(sig1)])
disp(['mean like cond 2 = ',num2str(mn2),' +/- ',num2str(sig2)])
nBetter = sum(HOLikes2<HOLikes1);
disp(['Second condition was better on ',num2str(nBetter),'/',num2str(K)])

keyboard
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% T = length(y);
% specSig = abs(fft(y)).^2;
% freqSig = linspace(0,1/2,floor(T/2));

% freqs = linspace(0,1/2,T);

% specEst = get_pSTFT_spec(freqs,lamxEst,varxEst,omEst);


% cumsumSig = 2*cumsum(specSig(1:floor(T/2))/T);
% cumsumEst = cumsum(sum(specEst,1));

% if dispFigs==1

%   figure
%   subplot(2,1,1)
%   hold on

%   ah1=plot(freqSig,specSig(1:floor(T/2))/T,'-k');
%   for d=1:D
%     ah3=plot(freqs,specEst(d,:),'-m');
%   end
%   legend([ah1,ah3],'data','estimated')
%   set(gca,'yscale','log')
  
%   subplot(2,1,2)
%   hold on
%   ah1=plot(freqSig,cumsumSig,'-k');
%   ah3=plot(freqs,cumsumEst,'-m');
%   legend([ah1,ah3],'data','estimated')
% end

% tol = 5e-1;
% assertElementsAlmostEqual(100+cumsumEst(1:2:end)',100+cumsumSig,'relative',tol,0)
