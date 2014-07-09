clear;

%filename = 'trainDenoise7_D40_K20';
%filename = 'trainDenoise5_D40';
%filename = 'trainDenoise8_D40_long'; %
%filename = 'trainDenoise11_D40'; %
filename = 'trainDenoise10_D40'; %
load(['/home/ret26/data/probFB/nmf/',filename,'.mat'])

lenx2 = 750*ones(K,1);%lenx2 = 200*lenx;
varx2 = 1.0*varx;
mux2 = 7.5+mux;
savename = [filename,'_denoise_results_sentence_long_lenx2_750_mux2_75_varx_100pc_its_40_15.mat'];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TESTING

RngLimTest = round([fs*13.7+1,17.1*fs]);  % two sentences
%RngLimTest = round([fs*13.7+1,15.4*fs]);  % Picks a sentence
%RngLimTest = round([fs*13.73+1,14.28*fs]);
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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Denoising
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

L = 10;
varys = [logspace(log10(1e-2),log10(10),L)];

%L=1;
%varys = 1;

HInit = exp(randn(T,K))/10000;

numMethods = 7; % includes the un-denoised signal
snr_a = zeros(L,D,numMethods);
snr_loga = zeros(L,D,numMethods);
snr_y = zeros(L,numMethods);
pesq_y = zeros(L,numMethods);

optstNMF.numIts = 1000;
optstNMF.progress_chunk = 500;

optsGTFtNMF.numIts = 40;
optsGTFtNMF.progress_chunk = 15;

optsGTFtNMF2.numIts = 15;
optsGTFtNMF2.progress_chunk = 15;

optsNMF.numIts = 1000;
optsNMF.progress_chunk = 500;


% Uncomment for quick running to test code
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
%load('/home/ret26/data/probFB/nmf/trainDenoise5_D40_denoise_results_sentence_lenx2_750_mux2_75_varx_80pc.mat');
%load('/home/ret26/data/probFB/nmf/trainDenoise11_D40_denoise_results_sentence_medium_new_init_lenx2_750_mux2_75_varx_100pc_its_75.mat');

%load(savename)
for l=[10,9,8,5,1,2,6,3,4,7]
  
  % display progress to user
  cnt= cnt+1;
  disp(['Progress ',num2str(cnt),'/',num2str(L)]);

  % set the state 
  randn('state',1);
  
  % noisy signal
  yNoisy = yTest + randn(T,1)*sqrt(varys(l));
  ZNoisy = probFB(yNoisy,Lam1,Var1,om,0); 
  ANoisy = abs(ZNoisy').^2;

  % figure out the noise levels -- this is done empirically, but we
  % could do this analytically if needed
  filt_noise = probFB(randn(T,1)*sqrt(varys(l)),Lam1,Var1,om,0);
  varyCur = ones(T,1)*var(filt_noise');
  
  % denoise using NMF
  [HTest_fp,infoFP] = nmf_inf_fp(ANoisy,WEst1,HInit,varyCur);
  [HTestNMF,info] =  nmf_inf(ANoisy,WEst1,HTest_fp,[],[],[],varyCur,optsNMF);
  ADenoiseNMF = HTestNMF*WEst1;
  [yDenoiseNMF,aRec1,info3] = recon_FB_mag(yInit,ADenoiseNMF.^0.5,spec,[],lamRec,varyRec);
  
  % denoise using temporal NMF
  [HTesttNMF,info] = tnmf_inf(ANoisy,WEst2,HInit,lenx,mux,varx,varyCur,optstNMF);
  ADenoisetNMF = HTesttNMF*WEst2;
  [yDenoisetNMF,aRec2,info4] = recon_FB_mag(yInit,ADenoisetNMF.^0.5,spec,[],lamRec,varyRec);


  % denoise using untrained GTF
  ZDenoise = probFB(yNoisy,Lam1,Var1,om,varys(l)); 
  yDenoiseGTF_UT = sum(real(ZDenoise'),2);
  ZTemp = probFB(yDenoiseGTF_UT,Lam1,Var1,om,0); 
  ADenoiseGTF_UT = abs(ZTemp').^2;

  % denoise using GTF
  ZDenoise = probFB(yNoisy,Lam1Fit,Var1Fit,omFit,varys(l)); 
  yDenoiseGTF = sum(real(ZDenoise'),2);
  ZTemp = probFB(yDenoiseGTF,Lam1,Var1,om,0); 
  ADenoiseGTF = abs(ZTemp').^2;


  % denoise using GTFtNMF
 % HInitGTFtNMF = HTesttNMF;
  HInitGTFtNMF = exp(ones(T,1)*mux2');
%HInitGTFtNMF = exp(-12*ones(T,K)); %*mux2'
[HEst_GTFtNMF,mnV,covV,info5] = GTFtNMF_inf(yNoisy,WEst3, ...
					      HInitGTFtNMF,Lam1,Var1, ...
					      om,varys(l),lenx2,mux2,varx2,optsGTFtNMF);

  
  yDenoiseGTFtNMF = sum(real(mnV'),2);
  ZTemp = probFB(yDenoiseGTFtNMF,Lam1,Var1,om,0); 
  ADenoiseGTFtNMF = abs(ZTemp').^2;

    % denoise using GTFtNMF w/ trained front end filters
  HInitGTFtNMF = HEst_GTFtNMF;
%  HInitGTFtNMF = exp(ones(T,1)*mux2'-2);
%HInitGTFtNMF = exp(-12*ones(T,K)); %*mux2'
[HEst_GTFtNMF2,mnV,covV,info5] = GTFtNMF_inf(yNoisy,WEst4, ...
					      HInitGTFtNMF,Lam14,Var14, ...
					      om4,varys(l),lenx2,mux2,varx2,optsGTFtNMF2);

  
  yDenoiseGTFtNMF2 = sum(real(mnV'),2);
  ZTemp = probFB(yDenoiseGTFtNMF2,Lam1,Var1,om,0); 
  ADenoiseGTFtNMF2 = abs(ZTemp').^2;

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  Ys(:,:,l) = single([yNoisy,yDenoiseNMF,yDenoisetNMF,yDenoiseGTF_UT,yDenoiseGTF,yDenoiseGTFtNMF,yDenoiseGTFtNMF2]);
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % figure out the amount NMF has denoised the spectrogram by
  
  snr_a(l,:,1) = snr(ATest,ANoisy);
  snr_a(l,:,2) = snr(ATest,ADenoiseNMF);
  snr_a(l,:,3) = snr(ATest,ADenoisetNMF);
  snr_a(l,:,4) = snr(ATest,ADenoiseGTF_UT);
  snr_a(l,:,5) = snr(ATest,ADenoiseGTF);
  snr_a(l,:,6) = snr(ATest,ADenoiseGTFtNMF);
  snr_a(l,:,7) = snr(ATest,ADenoiseGTFtNMF2);

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

  
  snr_loga(l,:,1) = snr(log10(loud_floor(ATest,deltaSNR)), ...
			log10(loud_floor(ANoisy,deltaSNR)));
  snr_loga(l,:,2) = snr(log10(loud_floor(ATest,deltaSNR)), ...
			log10(loud_floor(ADenoiseNMF,deltaSNR)));
  snr_loga(l,:,3) = snr(log10(loud_floor(ATest,deltaSNR)), ...
			log10(loud_floor(ADenoisetNMF,deltaSNR)));
  snr_loga(l,:,4) = snr(log10(loud_floor(ATest,deltaSNR)), ...
			log10(loud_floor(ADenoiseGTF_UT,deltaSNR)));
  snr_loga(l,:,5) = snr(log10(loud_floor(ATest,deltaSNR)), ...
			log10(loud_floor(ADenoiseGTF,deltaSNR)));
  snr_loga(l,:,6) = snr(log10(loud_floor(ATest,deltaSNR)), ...
			log10(loud_floor(ADenoiseGTFtNMF,deltaSNR)));
  snr_loga(l,:,7) = snr(log10(loud_floor(ATest,deltaSNR)), ...
			log10(loud_floor(ADenoiseGTFtNMF2,deltaSNR)));


  snr_y(l,1) = snr(yTest,yNoisy);
  snr_y(l,2) = snr(yTest,yDenoiseNMF);
  snr_y(l,3) = snr(yTest,yDenoisetNMF);
  snr_y(l,4) = snr(yTest,yDenoiseGTF_UT);
  snr_y(l,5) = snr(yTest,yDenoiseGTF);
  snr_y(l,6) = snr(yTest,yDenoiseGTFtNMF);
  snr_y(l,7) = snr(yTest,yDenoiseGTFtNMF2);

  pesq_y(l,1) = pesq(yTest, yNoisy, fs);
  pesq_y(l,2) = pesq(yTest, yDenoiseNMF, fs);
  pesq_y(l,3) = pesq(yTest, yDenoisetNMF, fs);
  pesq_y(l,4) = pesq(yTest, yDenoiseGTF_UT, fs);
  pesq_y(l,5) = pesq(yTest, yDenoiseGTF, fs);
  pesq_y(l,6) = pesq(yTest, yDenoiseGTFtNMF, fs);
  pesq_y(l,7) = pesq(yTest, yDenoiseGTFtNMF2, fs);

  disp(savename)  
  snr_a_mn = squeeze(mean(snr_a,2))
  snr_loga_mn = squeeze(mean(snr_loga,2))
  snr_y
  pesq_y
  
  % figure
  % subplot(5,1,1)
  % plot(yTest)

  % subplot(5,1,2)
  % plot(yDenoiseGTFtNMF)
  
  % subplot(5,1,3)
  % plot(yDenoisetNMF)

  % subplot(5,1,4)
  % plot(yDenoiseGTF)

  % subplot(5,1,5)
  % plot(yNoisy)
  
  % figure
  % subplot(5,1,1)
  % imagesc(log(ATest)')

  % subplot(5,1,2)
  % imagesc(log(ADenoiseGTFtNMF)')

  % subplot(5,1,3)
  % imagesc(log(ADenoisetNMF)')

  % subplot(5,1,4)
  % imagesc(log(ADenoiseGTF)')

  % subplot(5,1,5)
  % imagesc(log(ANoisy)')

  % drawnow
  % pause(0.1)

  save(['/home/ret26/data/probFB/nmf/',savename],'yTest','fs','WEst1','lenx','mux','varx','optstNMF','Lam1Fit','Var1Fit','omFit','WEst2','WEst3','Lam1','Var1','om','WEst4','Lam14','Var14','om4','varys','lenx2','mux2','varx2','optsGTFtNMF','Ys','snr_a_mn','snr_loga_mn','snr_a','snr_loga','snr_y','pesq_y')

  
  %  keyboard
end

figure
subplot(2,1,1)
hold on
title('envelope')
plot(mean(snr_a(:,:,1),2),mean(snr_a(:,:,2)-snr_a(:,:,1),2),'--k')
plot(mean(snr_a(:,:,1),2),mean(snr_a(:,:,3)-snr_a(:,:,1),2),'-k')
plot(mean(snr_a(:,:,1),2),mean(snr_a(:,:,4)-snr_a(:,:,1),2),'--b')
plot(mean(snr_a(:,:,1),2),mean(snr_a(:,:,5)-snr_a(:,:,1),2),'-b')
plot(mean(snr_a(:,:,1),2),mean(snr_a(:,:,6)-snr_a(:,:,1),2),'--r')
plot(mean(snr_a(:,:,1),2),mean(snr_a(:,:,7)-snr_a(:,:,1),2),'-r')
legend('NMF','tNMF','GTF untrained','GTF','GTFtNMF','GTFtNMF train FE')
xlabel('snr a before')
ylabel('snr a improvement')

subplot(2,1,2)
hold on
title('log-envelope')
plot(mean(snr_loga(:,:,1),2),mean(snr_loga(:,:,2)-snr_loga(:,:,1),2),'--k')
plot(mean(snr_loga(:,:,1),2),mean(snr_loga(:,:,3)-snr_loga(:,:,1),2),'-k')
plot(mean(snr_loga(:,:,1),2),mean(snr_loga(:,:,4)-snr_loga(:,:,1),2),'--b')
plot(mean(snr_loga(:,:,1),2),mean(snr_loga(:,:,5)-snr_loga(:,:,1),2),'-b')
plot(mean(snr_loga(:,:,1),2),mean(snr_loga(:,:,6)-snr_loga(:,:,1),2),'--r')
plot(mean(snr_loga(:,:,1),2),mean(snr_loga(:,:,7)-snr_loga(:,:,1),2),'-r')

legend('NMF','tNMF','GTF untrained','GTF','GTFtNMF','GTFtNMF train FE')
xlabel('snr log(a) before')
ylabel('snr log(a) improvement')


figure
subplot(2,1,1)
hold on
title('waveform')
plot(pesq_y(:,1),pesq_y(:,2)-pesq_y(:,1),'--k')
plot(pesq_y(:,1),pesq_y(:,3)-pesq_y(:,1),'-k')
plot(pesq_y(:,1),pesq_y(:,4)-pesq_y(:,1),'--b')
plot(pesq_y(:,1),pesq_y(:,5)-pesq_y(:,1),'-b')
plot(pesq_y(:,1),pesq_y(:,6)-pesq_y(:,1),'--r')
plot(pesq_y(:,1),pesq_y(:,7)-pesq_y(:,1),'-r')

legend('NMF','tNMF','GTF untrained','GTF','GTFtNMF','GTFtNMF train FE')
xlabel('PSEQ before')
ylabel('PSEQ improvement')

subplot(2,1,2)
hold on
plot(snr_y(:,1),snr_y(:,2)-snr_y(:,1),'--k')
plot(snr_y(:,1),snr_y(:,3)-snr_y(:,1),'-k')
plot(snr_y(:,1),snr_y(:,4)-snr_y(:,1),'--b')
plot(snr_y(:,1),snr_y(:,5)-snr_y(:,1),'-b')
plot(snr_y(:,1),snr_y(:,6)-snr_y(:,1),'--r')
plot(snr_y(:,1),snr_y(:,7)-snr_y(:,1),'-r')

legend('NMF','tNMF','GTF untrained','GTF','GTFtNMF','GTFtNMF train FE')
xlabel('SNR before')
ylabel('SNR improvement')


%[mean(snr_a-snr_orig_a,2)';
%mean(snr_loga-snr_orig_loga,2)']
%[pesq_before,pesq_after-pesq_before]
%[snr_orig_y,snr_y-snr_orig_y]







