% makes the denoising figure for the transactions on signal
% processing paper

clear;

printOutput = 1;

loadDir = '/home/rich/Data/probFB/nmf/';
%simName = 'trainDenoise5_D40_denoise_results_sentence_lenx2_750.mat';
%simName = 'trainDenoise5_D40_denoise_results_sentence_lenx2_750_mux2_75_varx_100pc_its_75.mat';


%simName = 'trainDenoise5_D40_denoise_results_sentence_lenx2_750_mux2_75_varx_100pc_its_100.mat';
%simName = 'trainDenoise5_D40_denoise_results_sentence_long_lenx2_750_mux2_75_varx_100pc_its_100.mat'
%simName = 'trainDenoise8_D40_long_denoise_results_sentence_long_lenx2_750_mux2_75_varx_100pc_its_100.mat'
%simName = 'trainDenoise11_D40_denoise_results_sentence_medium_new_init_more_its_lenx2_750_mux2_75_varx_100pc_its_75.mat';
%simName = 'trainDenoise11_D40_denoise_results_sentence_medium_new_init_lenx2_750_mux2_75_varx_100pc_its_75.mat';
%simName = 'trainDenoise11_D40_denoise_results_sentence_med_lenx2_750_mux2_75_varx_100pc_its_15.mat';

%simName = 'trainDenoise10_D40_denoise_results_sentence_long_lenx2_750_mux2_75_varx_100pc_its_15.mat';

%simName = 'trainDenoise10_D40_denoise_results_sentence_long_lenx2_750_mux2_75_varx_100pc_its_25.mat';


% SECOND BEST
%simName = 'trainDenoise11_D40_denoise_results_sentence_med_lenx2_750_mux2_75_varx_100pc_its_25.mat';

% BEST SO FAR
simName = 'trainDenoise10_D40_denoise_results_sentence_long_lenx2_750_mux2_75_varx_100pc_its_40_15.mat'

load([loadDir,simName])

saveDir = '/home/rich/Synchronised/Writings/mpad/tsp/figs/';
saveName = 'denoising_fig_new.eps'

Ys= double(Ys);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fonts
FontName = 'Times';
FSsm = 6;
FSmed = 7;
FSlg = 9;
%FSvlg = 13;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Line styles
colNoisy = [1,1,1]*0.7;
colGTFtNMF = [1,0,0];
colGTF = [0,0,1];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% preliminary figures for looking at things

% figure 
% subplot(2,1,1)
% hold on
% title('waveform')  
% plot(pesq_y(:,1),pesq_y(:,2)-pesq_y(:,1),'--k')
% plot(pesq_y(:,1),pesq_y(:,3)-pesq_y(:,1),'-k')
% plot(pesq_y(:,1),pesq_y(:,4)-pesq_y(:,1),'-b')
% plot(pesq_y(:,1),pesq_y(:,5)-pesq_y(:,1),'-r')
% legend('NMF','tNMF','GTF','GTFtNMF')
% xlabel('PSEQ before')
% ylabel('PSEQ improvement')

% subplot(2,1,2)
% hold on
% plot(snr_y(:,1),snr_y(:,2)-snr_y(:,1),'--k')
% plot(snr_y(:,1),snr_y(:,3)-snr_y(:,1),'-k')
% plot(snr_y(:,1),snr_y(:,4)-snr_y(:,1),'-b')
% plot(snr_y(:,1),snr_y(:,5)-snr_y(:,1),'-r')
% legend('NMF','tNMF','GTF','GTFtNMF')
% xlabel('SNR before')
% ylabel('SNR improvement')                          

% figure
% hold on
% plot(Ys(:,1,7),'-','color',colNoisy)
% plot(Ys(:,4,7),'-','color',colGTF)
% plot(Ys(:,5,7),'-','color',colGTFtNMF)
% plot(yTest,'-k','linewidth',2)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Snippets and noise-levels to plot

nsip_lim = [15260,15600;  2037,2329; 6011,6234];

noise_level = [6,8];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Axis positions

left = 0.075;
right = 0.005;
top = 0.0075;
bottom = 0.075;

vspace = 0.1;
vspaceSN = 0.035;
hspace = 0.1;

width = (1-left-right-2*hspace)/3;
heightTop = 0.25;
height = (1-top-bottom-2*vspaceSN-vspace-heightTop)/3;

across = [width+hspace,0,0,0];

pos11 = [left,1-top-heightTop,width,heightTop];
pos12 = pos11+across;
pos13 = pos12+across;

widthSnip = (1-left-right-hspace)/2;
acrossSnip = [widthSnip+hspace,0,0,0];
down = -[0,height+vspaceSN,0,0];

pos21 = [left,1-top-heightTop-vspace-height,widthSnip,height];
pos22 = pos21+acrossSnip;
%pos23 = pos22+across;

pos31 = pos21+down;
pos32 = pos31+acrossSnip;
%pos33 = pos32+across;

pos41 = pos31+down;
pos42 = pos41+acrossSnip;
%pos43 = pos42+across;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Setup the figure settings

figure1=figure;
PP = [0,0,17.00,12.40]; %*** paper position in centimeters
PS = PP(end-1:end); % paper size in centimeters

set(figure1,'paperpositionmode','manual','paperposition', ...
        PP,'papersize',PS, 'paperunits','centimeters');

% So the figure is the same size on the screen as when it is printed:
pu = get(gcf,'PaperUnits');
pp = get(gcf,'PaperPosition');
set(gcf,'Units',pu,'Position',pp)

PR = PS(1)/PS(2);

colormap gray

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% axis 11 -- snr y

ax11 = axes('position',pos11); % produce axis
hold on

delta1= snr_y(:,2)-snr_y(:,1);
delta2= snr_y(:,3)-snr_y(:,1);
delta3 = snr_y(:,4)-snr_y(:,1);
delta4 = snr_y(:,5)-snr_y(:,1);
delta5 = snr_y(:,6)-snr_y(:,1);
delta6 = snr_y(:,7)-snr_y(:,1);

disp(['waveform improvement (lGTF/lGTFtNMF)   ',num2str(mean(delta6-delta4)),'dB'])
disp(['waveform improvement (GTF/lGTF)   ',num2str(mean(delta4-delta3)),'dB'])
disp(['waveform improvement (GTFtNMF/lGTFtNMF)   ',num2str(mean(delta6-delta5)),'dB'])

h01 = plot([1,0],[1,1],'--','color',[1,1,1]*0.6);
h11 = plot([1,0],[1,1],'-','color',[1,1,1]*0.6);

h1 =plot(snr_y(:,1),snr_y(:,2)-snr_y(:,1),'-.','color',[0,0.6,0]);
h2 =plot(snr_y(:,1),snr_y(:,3)-snr_y(:,1),'-','color',[0,0.6,0]);
h3 =plot(snr_y(:,1),snr_y(:,4)-snr_y(:,1),'--b');
h4 =plot(snr_y(:,1),snr_y(:,5)-snr_y(:,1),'-b');
h5 =plot(snr_y(:,1),snr_y(:,6)-snr_y(:,1),'--r');
h6 =plot(snr_y(:,1),snr_y(:,7)-snr_y(:,1),'-r');

%if printOutput == 1
%  legend([h1,h2,h4,h6,h01,h11]',{'NMF','tNMF','GTF','GTFtNMF','unadapted filters','adapted filters'},'location','southwest')
%end

delete(h01); delete(h11);

for k=1:length(noise_level)
  plot(snr_y(noise_level(k),1),delta4(noise_level(k)),'ob')
  plot(snr_y(noise_level(k),1),delta6(noise_level(k)),'or')
end

%legend('NMF','tNMF','GTF','GTFtNMF','location','southwest')
xlabel('SNR before /dB','fontname',FontName,'FontSize',FSlg)
ylabel('SNR improvement /dB','fontname',FontName,'FontSize',FSmed)                          
set(gca,'fontname',FontName,'FontSize',FSsm)
set(gca,'TickDir','out')
set(gca,'xlim',[min(snr_y(:,1))-1,max(snr_y(:,1))])
%set(gca,'ylim',[min([delta1;delta2;delta3;delta4;delta5;delta6]),max([delta1;delta2;delta3;delta4;delta5;delta6])])
set(gca,'ylim',[max([0,min([delta1;delta2;delta3;delta4;delta5;delta6])]),max([delta1;delta2;delta3;delta4;delta5;delta6])])


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% axis 12 -- pesq y

ax12 = axes('position',pos12); % produce axis
hold on

delta1= pesq_y(:,2)-pesq_y(:,1);
delta2= pesq_y(:,3)-pesq_y(:,1);
delta3 = pesq_y(:,4)-pesq_y(:,1);
delta4 = pesq_y(:,5)-pesq_y(:,1);
delta5 = pesq_y(:,6)-pesq_y(:,1);
delta6 = pesq_y(:,7)-pesq_y(:,1);

disp(['pesq improvement over GTF   ',num2str(mean(delta6-delta4)),'units'])
disp(['pesq improvement (GTF/lGTF)   ',num2str(mean(delta4-delta3)),'dB'])
disp(['pesq improvement (GTFtNMF/lGTFtNMF)   ',num2str(mean(delta6-delta5)),'dB'])

for k=1:length(noise_level)
  plot(pesq_y(noise_level(k),1),delta4(noise_level(k)),'ob')
  plot(pesq_y(noise_level(k),1),delta6(noise_level(k)),'or')
end

plot(pesq_y(:,1),pesq_y(:,2)-pesq_y(:,1),'-.','color',[0,0.6,0])
plot(pesq_y(:,1),pesq_y(:,3)-pesq_y(:,1),'-','color',[0,0.6,0])
plot(pesq_y(:,1),pesq_y(:,4)-pesq_y(:,1),'--b')
plot(pesq_y(:,1),pesq_y(:,5)-pesq_y(:,1),'-b')
plot(pesq_y(:,1),pesq_y(:,6)-pesq_y(:,1),'--r')
plot(pesq_y(:,1),pesq_y(:,7)-pesq_y(:,1),'-r')

xlabel('PESQ before','fontname',FontName,'FontSize',FSlg)
ylabel('PESQ improvement','fontname',FontName,'FontSize',FSmed)               
set(gca,'fontname',FontName,'FontSize',FSsm)
set(gca,'TickDir','out')

set(gca,'xlim',[min(pesq_y(:,1))-0.1,max(pesq_y(:,1))])
set(gca,'ylim',[min([delta1;delta2;delta3;delta4;delta5;delta6]),max([delta1;delta2;delta3;delta4;delta5;delta6])])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% axis 13 -- log amn
delta1= snr_loga_mn(:,2)-snr_loga_mn(:,1);
delta2= snr_loga_mn(:,3)-snr_loga_mn(:,1);
delta3 = snr_loga_mn(:,4)-snr_loga_mn(:,1);
delta4 = snr_loga_mn(:,5)-snr_loga_mn(:,1);
delta5 = snr_loga_mn(:,6)-snr_loga_mn(:,1);
delta6 = snr_loga_mn(:,7)-snr_loga_mn(:,1);


ax13 = axes('position',pos13); % produce axis
hold on
plot(snr_loga_mn(:,1),snr_loga_mn(:,2)-snr_loga_mn(:,1),'-.','color',[0,0.6,0])
plot(snr_loga_mn(:,1),snr_loga_mn(:,3)-snr_loga_mn(:,1),'-','color',[0,0.6,0])
plot(snr_loga_mn(:,1),snr_loga_mn(:,4)-snr_loga_mn(:,1),'--b')
plot(snr_loga_mn(:,1),snr_loga_mn(:,5)-snr_loga_mn(:,1),'-b')
plot(snr_loga_mn(:,1),snr_loga_mn(:,6)-snr_loga_mn(:,1),'--r')
plot(snr_loga_mn(:,1),snr_loga_mn(:,7)-snr_loga_mn(:,1),'-r')

for k=1:length(noise_level)
  plot(snr_loga_mn(noise_level(k),1),delta4(noise_level(k)),'ob')
  plot(snr_loga_mn(noise_level(k),1),delta6(noise_level(k)),'or')
end

xlabel('SNR log-spec before /dB','fontname',FontName,'FontSize',FSlg)
ylabel('SNR log-spec improvement /dB','fontname',FontName,'FontSize',FSmed) ...
    
    
disp(['loga SNR improvement over tNMF   ',num2str(mean(delta6-delta2)),'dB'])
disp(['loga SNR improvement (GTF/lGTF)   ',num2str(mean(delta4-delta3)),'dB'])
disp(['loga SNR improvement (GTFtNMF/lGTFtNMF)   ',num2str(mean(delta6-delta5)),'dB'])

set(gca,'fontname',FontName,'FontSize',FSsm)
set(gca,'TickDir','out')

set(gca,'xlim',[min(snr_loga_mn(:,1))-1,max(snr_loga_mn(:,1))])
set(gca,'ylim',[min([delta1;delta2;delta3;delta4;delta5;delta6]),max([delta1;delta2;delta3;delta4;delta5;delta6])])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SNIPPETS

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Snippet 1

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% axis 21 -- snippet 1, noise level 1
ax21 = axes('position',pos21); % produce axis
hold on

hold on
tim1 = 1000*(1:1+nsip_lim(1,2)-nsip_lim(1,1))/fs;
plot(tim1,Ys(nsip_lim(1,1):nsip_lim(1,2),1,noise_level(1)),'-','color',colNoisy)
plot(tim1,yTest(nsip_lim(1,1):nsip_lim(1,2)),'-k','linewidth',1)
plot(tim1,Ys(nsip_lim(1,1):nsip_lim(1,2),5,noise_level(1)),'-','color',colGTF)
plot(tim1,Ys(nsip_lim(1,1):nsip_lim(1,2),7,noise_level(1)),'-','color',colGTFtNMF)

%xlabel('time /ms','fontname',FontName,'FontSize',FSlg)
ylabel('y_t','fontname',FontName,'FontSize',FSlg)               
set(gca,'fontname',FontName,'FontSize',FSsm)
set(gca,'TickDir','out')
set(gca,'xlim',[tim1(1),tim1(end)])

ylim = [-1,1]*1.01*max(abs(Ys(nsip_lim(1,1):nsip_lim(1,2),1,noise_level(1))));

set(gca,'ylim',ylim)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% axis 22 -- snippet 1, noise level 2
ax22 = axes('position',pos22); % produce axis
hold on

hold on
tim1 = 1000*(1:1+nsip_lim(1,2)-nsip_lim(1,1))/fs;
plot(tim1,Ys(nsip_lim(1,1):nsip_lim(1,2),1,noise_level(2)),'-','color',colNoisy)
plot(tim1,yTest(nsip_lim(1,1):nsip_lim(1,2)),'-k','linewidth',1)
plot(tim1,Ys(nsip_lim(1,1):nsip_lim(1,2),5,noise_level(2)),'-','color',colGTF)
plot(tim1,Ys(nsip_lim(1,1):nsip_lim(1,2),7,noise_level(2)),'-','color',colGTFtNMF)

%xlabel('time /ms','fontname',FontName,'FontSize',FSlg)
ylabel('y_t','fontname',FontName,'FontSize',FSlg)               
set(gca,'fontname',FontName,'FontSize',FSsm)
set(gca,'TickDir','out')
set(gca,'xlim',[tim1(1),tim1(end)])

ylim = [-1,1]*1.01*max(abs(Ys(nsip_lim(1,1):nsip_lim(1,2),1,noise_level(2))));

set(gca,'ylim',ylim)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Snippet 2

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% axis 31 -- snippet 2, noise level 1
ax31 = axes('position',pos31); % produce axis
hold on

hold on
tim1 = 1000*(1:1+nsip_lim(2,2)-nsip_lim(2,1))/fs;
plot(tim1,Ys(nsip_lim(2,1):nsip_lim(2,2),1,noise_level(1)),'-','color',colNoisy)
plot(tim1,yTest(nsip_lim(2,1):nsip_lim(2,2)),'-k','linewidth',1)
plot(tim1,Ys(nsip_lim(2,1):nsip_lim(2,2),5,noise_level(1)),'-','color',colGTF)
plot(tim1,Ys(nsip_lim(2,1):nsip_lim(2,2),7,noise_level(1)),'-','color',colGTFtNMF)

%xlabel('time /ms','fontname',FontName,'FontSize',FSlg)
ylabel('y_t','fontname',FontName,'FontSize',FSlg)               
set(gca,'fontname',FontName,'FontSize',FSsm)
set(gca,'TickDir','out')
set(gca,'xlim',[tim1(1),tim1(end)])

ylim = [-1,1]*1.01*max(abs(Ys(nsip_lim(2,1):nsip_lim(2,2),1,noise_level(1))));

set(gca,'ylim',ylim)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% axis 32 -- snippet 2, noise level 2
ax32 = axes('position',pos32); % produce axis
hold on

hold on
tim1 = 1000*(1:1+nsip_lim(2,2)-nsip_lim(2,1))/fs;
plot(tim1,Ys(nsip_lim(2,1):nsip_lim(2,2),1,noise_level(2)),'-','color',colNoisy)
plot(tim1,yTest(nsip_lim(2,1):nsip_lim(2,2)),'-k','linewidth',1)
plot(tim1,Ys(nsip_lim(2,1):nsip_lim(2,2),5,noise_level(2)),'-','color',colGTF)
plot(tim1,Ys(nsip_lim(2,1):nsip_lim(2,2),7,noise_level(2)),'-','color',colGTFtNMF)

%xlabel('time /ms','fontname',FontName,'FontSize',FSlg)
ylabel('y_t','fontname',FontName,'FontSize',FSlg)               
set(gca,'fontname',FontName,'FontSize',FSsm)
set(gca,'TickDir','out')
set(gca,'xlim',[tim1(1),tim1(end)])

ylim = [-1,1]*1.01*max(abs(Ys(nsip_lim(2,1):nsip_lim(2,2),1,noise_level(2))));

set(gca,'ylim',ylim)



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Snippet 3

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% axis 41 -- snippet 3, noise level 1
ax41 = axes('position',pos41); % produce axis
hold on

hold on
tim1 = 1000*(1:1+nsip_lim(3,2)-nsip_lim(3,1))/fs;
plot(tim1,Ys(nsip_lim(3,1):nsip_lim(3,2),1,noise_level(1)),'-','color',colNoisy)
plot(tim1,yTest(nsip_lim(3,1):nsip_lim(3,2)),'-k','linewidth',1)
plot(tim1,Ys(nsip_lim(3,1):nsip_lim(3,2),5,noise_level(1)),'-','color',colGTF)
plot(tim1,Ys(nsip_lim(3,1):nsip_lim(3,2),7,noise_level(1)),'-','color',colGTFtNMF)

xlabel('time /ms','fontname',FontName,'FontSize',FSlg)
ylabel('y_t','fontname',FontName,'FontSize',FSlg)               
set(gca,'fontname',FontName,'FontSize',FSsm)
set(gca,'TickDir','out')
set(gca,'xlim',[tim1(1),tim1(end)])

ylim = [-1,1]*1.01*max(abs(Ys(nsip_lim(3,1):nsip_lim(3,2),1,noise_level(1))));

set(gca,'ylim',ylim)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% axis 42 -- snippet 1, noise level 2
ax42 = axes('position',pos42); % produce axis
hold on

hold on
tim1 = 1000*(1:1+nsip_lim(3,2)-nsip_lim(3,1))/fs;
plot(tim1,Ys(nsip_lim(3,1):nsip_lim(3,2),1,noise_level(2)),'-','color',colNoisy)
plot(tim1,yTest(nsip_lim(3,1):nsip_lim(3,2)),'-k','linewidth',1)
plot(tim1,Ys(nsip_lim(3,1):nsip_lim(3,2),5,noise_level(2)),'-','color',colGTF)
plot(tim1,Ys(nsip_lim(3,1):nsip_lim(3,2),7,noise_level(2)),'-','color',colGTFtNMF)

xlabel('time /ms','fontname',FontName,'FontSize',FSlg)
ylabel('y_t','fontname',FontName,'FontSize',FSlg)               
set(gca,'fontname',FontName,'FontSize',FSsm)
set(gca,'TickDir','out')
set(gca,'xlim',[tim1(1),tim1(end)])

ylim = [-1,1]*1.01*max(abs(Ys(nsip_lim(3,1):nsip_lim(3,2),1,noise_level(2))));

set(gca,'ylim',ylim)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Printing figure

if printOutput==1
  print(figure1,'-depsc','-painters',[saveDir,saveName],'-loose');
  %plot2svg([filename,'ColBar','.svg'],figure2,'png')
  %str = ['! gv ',filename,'ColBar','.eps &'];
  %eval(str)

end

snr_a_mn 
snr_loga_mn
snr_y
pesq_y       



% How much difference does the floor make on the log amplitude SNR calculation:
% snr_loga_new = snr_loga;

% ZTest = probFB(yTest,Lam1,Var1,om,0);
% ATest = abs(ZTest').^2;

% for l=1:10
% %    Ys(:,:,l) = single([yNoisy,yDenoiseNMF,yDenoisetNMF,yDenoiseGTF_UT,yDenoiseGTF,yDenoiseGTFtNMF,yDenoiseGTFtNMF2]);

%     yDenoiseGTF = Ys(:,5,l);
%     yDenoiseGTFtNMF = Ys(:,6,l);
%     yDenoiseGTFtNMF2 = Ys(:,7,l);
    
%     ZTemp = probFB(yDenoiseGTF,Lam1,Var1,om,0);
%     ADenoiseGTF = abs(ZTemp').^2;
    
%     ZTemp = probFB(yDenoiseGTFtNMF,Lam1,Var1,om,0);
%     ADenoiseGTFtNMF = abs(ZTemp').^2;
 
%     ZTemp = probFB(yDenoiseGTFtNMF2,Lam1,Var1,om,0);
%     ADenoiseGTFtNMF2 = abs(ZTemp').^2;
    
%     deltaSNR = 1e-5;
    
%     snr_loga_new(l,:,5) = snr(log10(loud_floor(ATest,deltaSNR)), ...
% 			  log10(loud_floor(ADenoiseGTF,deltaSNR)));
%     snr_loga_new(l,:,6) = snr(log10(loud_floor(ATest,deltaSNR)), ...
% 			  log10(loud_floor(ADenoiseGTFtNMF,deltaSNR)));
%     snr_loga_new(l,:,7) = snr(log10(loud_floor(ATest,deltaSNR)), ...
% 			  log10(loud_floor(ADenoiseGTFtNMF2,deltaSNR)));
  
% end

% snr_loga_mn 
% snr_loga_new_mn = squeeze(mean(snr_loga_new,2))


% keyboard