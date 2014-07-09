% makes the denoising figure for the transactions on signal
% processing paper

clear;

printOutput = 0;

loadDir = '/home/rich/Data/probFB/nmf/';
simName = 'trainDenoise5_D40_denoise_results_sentence_long_lenx2_750_mux2_75_varx_100pc_its_100.mat'

load([loadDir,simName])

saveDir = '/home/rich/Synchronised/Writings/mpad/tsp/sounds/';
saveName = 'denoising'

Ys= double(Ys);

[T,M,L] = size(Ys);

for l=1:L

  yCur1  = Ys(:,1,l)/(1.1*max(abs(Ys(:,1,l))));
  yCur2  = Ys(:,4,l)/(1.1*max(abs(Ys(:,4,l))));
  yCur3  = Ys(:,5,l)/(1.1*max(abs(Ys(:,5,l))));
  yCur0  = yTest/(1.1*max(abs(yTest)));

%  scale = 1.01*max(abs([yCur1;yCur2;yCur3;yCur0]));
  scale = 1;
  wavwrite(yCur0/scale,fs,[saveDir,saveName,'_true_',num2str(l),'.wav']); 
  wavwrite(yCur1/scale,fs,[saveDir,saveName,'_noisy_',num2str(l),'.wav']); 
  wavwrite(yCur2/scale,fs,[saveDir,saveName,'_GTF_',num2str(l),'.wav']); 
  wavwrite(yCur3/scale,fs,[saveDir,saveName,'_GTFtNMF_',num2str(l),'.wav']); 

end
