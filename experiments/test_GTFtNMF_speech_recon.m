clear;

filename = 'trainDenoise11_D40';
load(['/home/ret26/data/probFB/nmf/',filename,'.mat'])
savename = [filename,'_recon_results_sentence_long_lenx2_750_mux2_75_varx2_100pc_its_75.mat'];

lenx2 = 750*ones(D,1);%200*lenx;
varx2 = 1*varx;
mux2 = 7.5+mux;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TESTING

 RngLimTest = round([fs*13.7+1,17.1*fs]);  
 gapPos = [1298,1985,2725,5850,6408,8000,9643,13130,14460,15710, ...
 	  19300,20600,22550, 32950,34420,35110,43800,44800,47810,49070,50000]; % Picks two sentences


%RngLimTest = round([fs*13.7+1,15.4*fs]);  
%% %gapPos = [1298,2100,2725,5400,6350,6900,8000,10000,13130,14660,16000,19000,20800,22550]; 
%gapPos = [1298,1985,2725,5850,6408,8000,9643,13130,14460,15710,19300,20600,22550]; % Picks a sentence

%RngLimTest = round([fs*13.73+1,14.28*fs]); gapPos = [1641,2725,5850,6408,8000];
%RngLimTest = round([fs*13.73+1,14*fs]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load signal and pre-process
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[y,fs] = wavread([soundPath,File,'.wav']); % reads in the file
yTest = y(RngLimTest(1):RngLimTest(2),1); 
yTest = resample(yTest, fs/DS, fs); % downsample the input
fs = fs/DS;
yTest = yTest/sqrt(var(yTest)); % rescale the input to unit variance
T = length(yTest);

% clean spectrogram computed here
ZTest = probFB(yTest,Lam1,Var1,om,0); 
ATest = abs(ZTest').^2;

[fmax,df,varMa] = probSpec2freq(om,Lam1,Var1);
tau = 1./df;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Reconstruction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%L = 10;
%gapLim = [5,300];
%NGaps = 10;


L = 10;
gapLim = [10,300];
NGaps = length(gapPos);

%gaps = ceil(logspace(log10(gapLim(1)),log10(gapLim(2)),L));
gaps = ceil(linspace(gapLim(1),gapLim(2),L));
%gapPos = gapLim(2)+ceil(rand(NGaps,1)*(T-2*gapLim(2)));

HInit = exp(randn(T,K))/10000;

numMethods = 7; % includes the un-denoised signal

snr_a = zeros(L,D,numMethods);
snr_loga = zeros(L,D,numMethods);
snr_y = zeros(L,numMethods);
pesq_y = zeros(L,numMethods);

optstNMF.numIts = 1000;
optstNMF.progress_chunk = 500;

%optsGTFtNMF.numIts = 75;
%optsGTFtNMF.progress_chunk = 25;

optsGTFtNMF.numIts = 75;
optsGTFtNMF.progress_chunk = 25;

optsNMF.numIts = 1000;
optsNMF.progress_chunk = 500;


% % Uncomment for quick running to test code
% optstNMF.numIts = 10;
% optstNMF.progress_chunk = 10;

% optsGTFtNMF.numIts = 10;
% optsGTFtNMF.progress_chunk = 10;

% optsNMF.numIts = 10;
% optsNMF.progress_chunk = 10;


% parameters for reconstruction of the signal for NMF/tNMF
spec = get_probFB_spec(Lam1,Var1,om,0,T);
lamRec = 0; varyRec = (1-lamRec^2);
yInit = randn(T,1)/1000;
Ys = zeros(T,numMethods,L);
cnt = 0;

load(savename)

for l=[10,1,5,2,6,9,3,8,4,7]
  
  % display progress to user
  cnt= cnt+1;
  disp(['Progress ',num2str(cnt),'/',num2str(L)]);

  % set the state 
  randn('state',1);
  
  % signal with missing data
  yGap = yTest;
  ind = [];
  for ng=1:NGaps
    ind = [ind,gapPos(ng)+[-ceil(gaps(l)/2):+ceil(gaps(l)/2)]];
  end

  yGap(ind) = 0;

  % the NMF variants need to ignore regions of the spectrogram
  % which overlap with the gap (due to the analysis filter)
  
  varyCur = zeros(T,D);
  fr = 0.05; % determined empirically to be optimal bandwidth
             % fraction to ignore when doing tNMF
  for d=1:D
    for ng=1:NGaps
      minind = max([ gapPos(ng)-ceil(gaps(l)/2+tau(d)*fr),1]);
      maxind = min([ gapPos(ng)+ceil(gaps(l)/2+tau(d)*fr),T]);
      indCur = [minind:maxind];
      varyCur(indCur,d) = 1e5;
    end
  end
  
  ZGap = probFB(yGap,Lam1,Var1,om,0); 
  AGap = abs(ZGap').^2;
     
  % denoise using NMF
  [HTest_fp,infoFP] = nmf_inf_fp(AGap,WEst1,HInit,varyCur);
  [HTestNMF,info] =  nmf_inf(AGap,WEst1,HTest_fp,[],[],[],varyCur,optsNMF);
  AReconNMF = HTestNMF*WEst1;
  [yReconNMF,aRec1,info3] = recon_FB_mag(yInit,AReconNMF.^0.5,spec,[],lamRec,varyRec);
  
  % denoise using temporal NMF
  [HTesttNMF,info] = tnmf_inf(AGap,WEst2,HInit,lenx,mux,varx,varyCur,optstNMF);
  ARecontNMF = HTesttNMF*WEst2;
  [yRecontNMF,aRec2,info4] = recon_FB_mag(yInit,ARecontNMF.^0.5,spec,[],lamRec,varyRec);

  varys = 1e-4*ones(T,1);
  varys(ind) = 1e5;
  
  % denoise using untrained GTF
  ZRecon = probFB(yGap,Lam1,Var1,om,varys); 
  yReconGTF_UT = sum(real(ZRecon'),2);
  ZTemp = probFB(yReconGTF_UT,Lam1,Var1,om,0); 
  AReconGTF_UT = abs(ZTemp').^2;

  % denoise using GTF
  ZRecon = probFB(yGap,Lam1Fit,Var1Fit,omFit,varys); 
  yReconGTF = sum(real(ZRecon'),2);
  ZTemp = probFB(yReconGTF,Lam1,Var1,om,0); 
  AReconGTF = abs(ZTemp').^2;

  % denoise using GTFtNMF
 % HInitGTFtNMF = HTesttNMF;
 % HInitGTFtNMF = exp(-8*ones(T,K)); %*mux2'
  HInitGTFtNMF = exp(ones(T,1)*mux2');
  [HEst_GTFtNMF,mnV,covV,info5] = GTFtNMF_inf(yGap,WEst3, ...
					      HInitGTFtNMF,Lam1,Var1, ...
					      om,varys,lenx2,mux2,varx2,optsGTFtNMF);
  
  yReconGTFtNMF = sum(real(mnV'),2);
  ZTemp = probFB(yReconGTFtNMF,Lam1,Var1,om,0); 
  AReconGTFtNMF = abs(ZTemp').^2;

  % denoise using GTFtNMF with trained front end filters
  %HInitGTFtNMF = HTesttNMF;
  HInitGTFtNMF = HEst_GTFtNMF;
  [HEst_GTFtNMF2,mnV,covV,info6] = GTFtNMF_inf(yGap,WEst4, ...
					      HInitGTFtNMF,Lam14,Var14, ...
					      om4,varys,lenx2,mux2,varx2,optsGTFtNMF);
  
  yReconGTFtNMF2 = sum(real(mnV'),2);
  ZTemp = probFB(yReconGTFtNMF2,Lam1,Var1,om,0); 
  AReconGTFtNMF2 = abs(ZTemp').^2;

  
  %AReconGTFtNMF = abs(mnV').^2;
  
  Ys(:,:,l) = single([yGap,yReconNMF,yRecontNMF,yReconGTF_UT,yReconGTF,yReconGTFtNMF,yReconGTFtNMF2]);
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % figure out the amount NMF has denoised the spectrogram by
  
  snr_a(l,:,1) = snr(ATest(ind,:),AGap(ind,:));
  snr_a(l,:,2) = snr(ATest(ind,:),AReconNMF(ind,:));
  snr_a(l,:,3) = snr(ATest(ind,:),ARecontNMF(ind,:));
  snr_a(l,:,4) = snr(ATest(ind,:),AReconGTF_UT(ind,:));
  snr_a(l,:,5) = snr(ATest(ind,:),AReconGTF(ind,:));    
  snr_a(l,:,6) = snr(ATest(ind,:),AReconGTFtNMF(ind,:));
  snr_a(l,:,7) = snr(ATest(ind,:),AReconGTFtNMF2(ind,:));

  deltaSNR = 1e-5; % apply loudness floor

  % code to check that the floor is only getting rid of a
  % relatively small proportion of the spectrogram
  %
  % figure
  % subplot(2,1,1)
  % imagesc(log(ATest)')
  
  % subplot(2,1,2)
  % [ATest_floor,frac] = loud_floor(ATest,deltaSNR);
  % frac
  % imagesc(log(ATest_floor)')

  
  snr_loga(l,:,1) = snr(log10(loud_floor(ATest(ind,:),deltaSNR)), ...
			log10(loud_floor(AGap(ind,:),deltaSNR)));
  snr_loga(l,:,2) = snr(log10(loud_floor(ATest(ind,:),deltaSNR)), ...
			log10(loud_floor(AReconNMF(ind,:),deltaSNR)));
  snr_loga(l,:,3) = snr(log10(loud_floor(ATest(ind,:),deltaSNR)), ...
			log10(loud_floor(ARecontNMF(ind,:),deltaSNR)));
  snr_loga(l,:,4) = snr(log10(loud_floor(ATest(ind,:),deltaSNR)), ...
			log10(loud_floor(AReconGTF_UT(ind,:),deltaSNR)));
  snr_loga(l,:,5) = snr(log10(loud_floor(ATest(ind,:),deltaSNR)), ...
			log10(loud_floor(AReconGTF(ind,:),deltaSNR)));
  snr_loga(l,:,6) = snr(log10(loud_floor(ATest(ind,:),deltaSNR)), ...
			log10(loud_floor(AReconGTFtNMF(ind,:),deltaSNR)));
  snr_loga(l,:,7) = snr(log10(loud_floor(ATest(ind,:),deltaSNR)), ...
			log10(loud_floor(AReconGTFtNMF2(ind,:),deltaSNR)));

  snr_y(l,1) = snr(yTest(ind,:),yGap(ind,:));
  snr_y(l,2) = snr(yTest(ind,:),yReconNMF(ind,:));
  snr_y(l,3) = snr(yTest(ind,:),yRecontNMF(ind,:));
  snr_y(l,4) = snr(yTest(ind,:),yReconGTF_UT(ind,:));
  snr_y(l,5) = snr(yTest(ind,:),yReconGTF(ind,:));
  snr_y(l,6) = snr(yTest(ind,:),yReconGTFtNMF(ind,:));
  snr_y(l,7) = snr(yTest(ind,:),yReconGTFtNMF2(ind,:));

  pesq_y(l,1) = pesq(yTest, yGap, fs);
  pesq_y(l,2) = pesq(yTest, yReconNMF, fs);
  pesq_y(l,3) = pesq(yTest, yRecontNMF, fs);
  pesq_y(l,4) = pesq(yTest, yReconGTF_UT, fs);
  pesq_y(l,5) = pesq(yTest, yReconGTF, fs);
  pesq_y(l,6) = pesq(yTest, yReconGTFtNMF, fs);
  pesq_y(l,7) = pesq(yTest, yReconGTFtNMF2, fs);


  disp(savename)    
  snr_a_mn = squeeze(mean(snr_a,2))
  snr_loga_mn = squeeze(mean(snr_loga,2))
  snr_y
  pesq_y
  
  % figure
  % subplot(5,1,1)
  % plot(yTest)

  % subplot(5,1,2)
  % plot(yReconGTFtNMF)
  
  % subplot(5,1,3)
  % plot(yRecontNMF)

  % subplot(5,1,4)
  % plot(yReconGTF)

  % subplot(5,1,5)
  % plot(yGap)
  
  % figure
  % subplot(5,1,1)
  % imagesc(log(ATest)')

  % subplot(5,1,2)
  % imagesc(log(AReconGTFtNMF)')

  % subplot(5,1,3)
  % imagesc(log(ARecontNMF)')

  % subplot(5,1,4)
  % imagesc(log(AReconGTF)')

  % subplot(5,1,5)
  % imagesc(log(AGap)')

  % drawnow
  % pause(0.1)
  

  save(['/home/ret26/data/probFB/nmf/',savename],'yTest','fs','WEst1','lenx','mux','varx','optstNMF','Lam1Fit','Var1Fit','omFit','WEst2','WEst3','Lam1','Var1','om','WEst4','Lam14','Var14','om4','varys','lenx2','mux2','varx2','optsGTFtNMF','Ys','snr_a_mn','snr_loga_mn','snr_a','snr_loga','snr_y','pesq_y','gaps','gapPos')
  %  keyboard
end

%Ys = single(Ys);

%clear ATest AReconGTFtNMF ARecontNMF  AReconGTF AGap mnV ...
%    covV ZRecon ZTemp HTestNMF HTesttNMF HInitGTFtNMF
%save(['/home/ret26/data/probFB/nmf/',savename])

figure
subplot(2,1,1)
hold on
title('envelope')
plot(gaps*1000/fs,mean(snr_a(:,:,2),2),'--k')
plot(gaps*1000/fs,mean(snr_a(:,:,3),2),'-k')
plot(gaps*1000/fs,mean(snr_a(:,:,4),2),'--b')
plot(gaps*1000/fs,mean(snr_a(:,:,5),2),'-b')
plot(gaps*1000/fs,mean(snr_a(:,:,6),2),'--r')
plot(gaps*1000/fs,mean(snr_a(:,:,7),2),'-r')
legend('NMF','tNMF','GTF untrained','GTF','GTFtNMF','GTFtNMF train FE')
xlabel('gap /ms')
ylabel('snr a')

subplot(2,1,2)
hold on
title('log-envelope')
plot(gaps*1000/fs,mean(snr_loga(:,:,2),2),'--k')
plot(gaps*1000/fs,mean(snr_loga(:,:,3),2),'-k')
plot(gaps*1000/fs,mean(snr_loga(:,:,4),2),'--b')
plot(gaps*1000/fs,mean(snr_loga(:,:,5),2),'-b')
plot(gaps*1000/fs,mean(snr_loga(:,:,6),2),'--r')
plot(gaps*1000/fs,mean(snr_loga(:,:,7),2),'-r')
legend('NMF','tNMF','GTF untrained','GTF','GTFtNMF','GTFtNMF train FE')
xlabel('gap /ms')
ylabel('snr log(a)')


figure
subplot(2,1,1)
hold on
title('waveform')
plot(gaps*1000/fs,pesq_y(:,2)-pesq_y(:,1),'--k')
plot(gaps*1000/fs,pesq_y(:,3)-pesq_y(:,1),'-k')
plot(gaps*1000/fs,pesq_y(:,4)-pesq_y(:,1),'--b')
plot(gaps*1000/fs,pesq_y(:,5)-pesq_y(:,1),'-b')
plot(gaps*1000/fs,pesq_y(:,6)-pesq_y(:,1),'--r')
plot(gaps*1000/fs,pesq_y(:,7)-pesq_y(:,1),'-r')
legend('NMF','tNMF','GTF untrained','GTF','GTFtNMF','GTFtNMF train FE')
xlabel('gap /ms')
ylabel('PSEQ improvement')

subplot(2,1,2)
hold on
plot(gaps*1000/fs,snr_y(:,2)-snr_y(:,1),'--k')
plot(gaps*1000/fs,snr_y(:,3)-snr_y(:,1),'-k')
plot(gaps*1000/fs,snr_y(:,4)-snr_y(:,1),'--b')
plot(gaps*1000/fs,snr_y(:,5)-snr_y(:,1),'-b')
plot(gaps*1000/fs,snr_y(:,6)-snr_y(:,1),'--r')
plot(gaps*1000/fs,snr_y(:,7)-snr_y(:,1),'-r')
legend('NMF','tNMF','GTF untrained','GTF','GTFtNMF','GTFtNMF train FE')
xlabel('gap /ms')
ylabel('SNR improvement')
