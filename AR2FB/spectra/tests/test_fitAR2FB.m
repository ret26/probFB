function test_suite = test_fitAR2FB
  initTestSuite;

% Tests:
%
% function [Lam,Var,Info] = fitAR2FB(y,D,varargin)
% 
  
function testAR2

% Tests the function runs on data generated from an AR2FB

dispFigs = 1;
T = 10000;

% D = 3;
% CF = [1/32;1/20;1/10];
% DF = [1/128;1/50;1/100];
% varMa = [1;4;2];

D = 4;
CF = [1/32;1/20;1/10;1/15];
DF = [1/128;1/100;1/100;1/500];
varMa = [1;4;2;2];

vary = 0;

[LamTr,VarTr] = freq2AR2(CF,DF,varMa);

[y,X] = sampleAR2FB(LamTr,VarTr,vary,T);
[yHO,X] = sampleAR2FB(LamTr,VarTr,vary,T);

opts.verbose = 1;
opts.yHO = yHO;

[LamEst,VarEst,Info] = fitAR2FB(y,D,opts);

[CFEst,dFEst,mVar] = AR22freq(LamEst,VarEst);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
specSig = abs(fft(y)).^2;
freqSig = linspace(0,1/2,floor(T/2));


specTr = zeros(D,T);
specEst = zeros(D,T);

for d=1:D
  [freq,specTr(d,:),fMAX,SpecMAX,dF1,dF2] = getSpecAR2(LamTr(d,:),VarTr(d),T,[0,1/2]);


  [freq,specEst(d,:),fMAX,SpecMAX,dF1,dF2] = getSpecAR2(LamEst(d,:),VarEst(d),T,[0,1/2]);

end

cumsumSig = 2*cumsum(specSig(1:floor(T/2))/T);
cumsumTr = cumsum(sum(specTr,1));
cumsumEst = cumsum(sum(specEst,1));

if dispFigs==1

  figure
  subplot(2,1,1)
  hold on

  ah1=plot(freqSig,specSig(1:floor(T/2))/T,'-k');
  for d=1:D
    ah2=plot(freq,specTr(d,:),'-g');
    ah3=plot(freq,specEst(d,:),'-m');
  end
  legend([ah1,ah2,ah3],'data','true','estimated')

  subplot(2,1,2)
  hold on
  ah1=plot(freqSig,cumsumSig,'-k');
  ah2=plot(freq,cumsumTr,'-g');
  ah3=plot(freq,cumsumEst,'-m');
  legend([ah1,ah2,ah3],'data','true','estimated')
end

tol = 5e-1;
assertElementsAlmostEqual(100+cumsumEst(1:2:end)',100+cumsumSig,'relative',tol,0)


function testSpeech

% Tests the function runs on speech - run using:
% $ runtests test_fitAR2FB:testSpeech
dispFigs=1;
File = '74 - Sentences.wav';
soundPath = '~/Synchronised/probFB/demos/signals/';
[yOrig,FSamp] = wavread([soundPath,File]);
RngLim = [FSamp*0.2+1,4.2*FSamp];
y = yOrig(RngLim(1):RngLim(2));
y = y/sqrt(var(y));
yHO = yOrig(RngLim(2)+(RngLim(1):RngLim(2)));
yHO = yHO/sqrt(var(yHO));

% %File = '87 - Wind and Leaves.wav';
%File = '80 - Fire.wav';
% %File = '05 - Bubbling stream.wav';
%File = '07 - Bullfrogs.wav';
%File = '06 - Budgerigars.wav';
%File = '24 - Fire Burning.wav';
% %File = '81 - Wind.wav';
% File = '82 - Rain.wav';

% soundPath = '~/Music/Various Artists/Best Ever Sound Effects - Vol.3 - Sounds Of Nature Sound Effects/';
% [yOrig,FSamp] = wavread([soundPath,File]);
% yOrig(:,1) = yOrig(:,1) - mean(yOrig(:,1));
% y = yOrig(1:FSamp);
% y = y/sqrt(var(y));
% yHO = yOrig(1+FSamp:FSamp*2);
% yHO = yHO/sqrt(var(yHO));

D = 15;

opts.verbose = 1;
opts.yHO = yHO;

[LamEst,VarEst,Info] = fitAR2FB(y,D,opts);

Info.likeHO(end)

[CFEst,dFEst,mVar] = AR22freq(LamEst,VarEst);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
T = length(y);
specSig = abs(fft(y)).^2;
freqSig = linspace(0,1/2,floor(T/2));

specEst = zeros(D,T);

for d=1:D
  [freq,specEst(d,:),fMAX,SpecMAX,dF1,dF2] = getSpecAR2(LamEst(d,:),VarEst(d),T,[0,1/2]);

end

cumsumSig = 2*cumsum(specSig(1:floor(T/2))/T);
cumsumEst = cumsum(sum(specEst,1));

if dispFigs==1

  figure
  subplot(2,1,1)
  hold on

  ah1=plot(freqSig,specSig(1:floor(T/2))/T,'-k');
  for d=1:D
    ah3=plot(freq,specEst(d,:),'-m');
  end
  legend([ah1,ah3],'data','estimated')
  set(gca,'yscale','log')
  
  subplot(2,1,2)
  hold on
  ah1=plot(freqSig,cumsumSig,'-k');
  ah3=plot(freq,cumsumEst,'-m');
  legend([ah1,ah3],'data','estimated')
end

tol = 5e-1;
assertElementsAlmostEqual(100+cumsumEst(1:2:end)',100+cumsumSig,'relative',tol,0)
