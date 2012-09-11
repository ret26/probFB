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
% % %File = '87 - Wind and Leaves.wav';
% %File = '80 - Fire.wav';
% % %File = '05 - Bubbling stream.wav';
% %File = '07 - Bullfrogs.wav';
% %File = '06 - Budgerigars.wav';
% %File = '24 - Fire Burning.wav';
% % %File = '81 - Wind.wav';
% % File = '82 - Rain.wav';

soundPath = '~/Music/Various Artists/Best Ever Sound Effects - Vol.3 - Sounds Of Nature Sound Effects/';%soundPath = '~/probFB/demos/signals/';
[yOrig,FSamp] = wavread([soundPath,File]);
RngLim = [FSamp*0.2+1,4.2*FSamp];
y = yOrig(RngLim(1):RngLim(2));
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
