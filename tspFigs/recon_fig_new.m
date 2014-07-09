% makes the denoising figure for the transactions on signal
% processing paper

clear;

printOutput = 1;

loadDir = '/home/rich/Data/probFB/nmf/';
%simName = 'trainDenoise5_D40_recon_results_sentence_lenx2_750.mat';
%simName = 'trainDenoise5_D40_recon_results_sentence_lenx2_750_varx2_half.mat';
%simName = ...
%    'trainDenoise5_D40_recon_results_sentence_D50_lenx2_750_mux2_75_varx2_100pc_its_100.mat';

%simName = 'trainDenoise5_D40_recon_results_sentence_more_gaps_lenx2_750_mux2_75_varx2_100pc_its_100.mat';

%simName = 'trainDenoise5_D40_recon_results_sentence_long_lenx2_750_mux2_75_varx2_100pc_its_75.mat';

%simName = 'trainDenoise11_D40_recon_results_sentence_medium_lenx2_750_mux2_75_varx2_100pc_its_75.mat';

simName = 'trainDenoise11_D40_recon_results_sentence_long_lenx2_750_mux2_75_varx2_100pc_its_75.mat';


load([loadDir,simName])

saveDir = '/home/rich/Synchronised/Writings/mpad/tsp/figs/';
saveName = 'reconstruction_fig_new_2.eps'

Ys= double(Ys);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fonts
FontName = 'Times';
FSsm = 6;
FSmed = 7;
FSlg = 9;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Line styles
colNoisy = [1,1,1]*0.9;
colGTFtNMF = [1,0,0];
colGTF = [0,0,1];
missCol = [.8,.8,0.8];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% preliminary figures for checking things look reasonable




plotGaps = 0;

if plotGaps==1

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
% plot(Ys(:,1,trial1),'-','color',colNoisy)
% plot(Ys(:,4,10),'-','color',colGTF)
% plot(Ys(:,5,10),'-','color',colGTFtNMF)
% plot(yTest,'-k','linewidth',2)

numGaps = length(gapPos);
trial1 = 8;
trial2 = 5;
snrGTF = zeros(numGaps,2);
snrGTFtNMF = zeros(numGaps,2);

for k=1:numGaps

  duration = max([gaps(trial1),gaps(trial2)]);
  indCur = gapPos(k)+[-ceil(1*duration):ceil(1*duration)];

  snrGTFtNMF1 =snr(yTest(indCur),Ys(indCur,7,trial1));
  snrGTF1 =snr(yTest(indCur),Ys(indCur,5,trial1));
  snrGTFtNMF2 =snr(yTest(indCur),Ys(indCur,7,trial2));
  snrGTF2 =snr(yTest(indCur),Ys(indCur,5,trial2));
  
  figure
  subplot(2,1,1)
  hold on
  title(['snr GTFtNMF: ',num2str(snrGTFtNMF1),'  snr GTF: ',num2str(snrGTF1)])
  plot(Ys(indCur,1,trial1),'-','color',colNoisy)
  plot(Ys(indCur,5,trial1),'-','color',colGTF)
  plot(Ys(indCur,7,trial1),'-','color',colGTFtNMF)
  plot(yTest(indCur),'-k','linewidth',2)

  subplot(2,1,2)
  hold on
  title(['snr GTFtNMF: ',num2str(snrGTFtNMF2),'  snr GTF: ',num2str(snrGTF2)])
  plot(Ys(indCur,1,trial2),'-','color',colNoisy)
  plot(Ys(indCur,5,trial2),'-','color',colGTF)
  plot(Ys(indCur,7,trial2),'-','color',colGTFtNMF)
  plot(yTest(indCur),'-k','linewidth',2)

  snrGTF(k,:) = [snrGTF1,snrGTF2];
  snrGTFtNMF(k,:) = [snrGTFtNMF1,snrGTFtNMF2];
end

figure
subplot(1,2,1)
hold on
plot(snrGTF(:,1),'.b')
plot(snrGTFtNMF(:,1),'.r')

subplot(1,2,2)
hold on
plot(snrGTF(:,2),'.b')
plot(snrGTFtNMF(:,2),'.r')

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Snippets and gap durations to plot


%[14,18,17,19||5,21,11]

nsip_lim = [gapPos(14)-200,gapPos(14)+200;
	    gapPos(19)-200,gapPos(19)+200;
	    gapPos(5)-200,gapPos(5)+200];

% nsip_lim = [gapPos(12)-200,gapPos(12)+200;
% 	    gapPos(10)-200,gapPos(10)+200;
% 	    gapPos(5)-200,gapPos(5)+200];

gaps_plot = [5,8];

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
% raw numbers

snr_a_mn 
snr_loga_mn
snr_y
pesq_y       

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

gapsms = gaps*1000/fs;
delta1= snr_y(:,2);
delta2= snr_y(:,3);
delta3 = snr_y(:,4);
delta4 = snr_y(:,5);
delta5 = snr_y(:,6);
delta6 = snr_y(:,7);

disp(['waveform improvement (lGTF/lGTFtNMF)   ',num2str(mean(delta6-delta4)),'dB'])
disp(['waveform improvement (GTF/lGTF)   ',num2str(mean(delta4-delta3)),'dB'])
disp(['waveform improvement (GTFtNMF/lGTFtNMF)   ',num2str(mean(delta6-delta5)),'dB'])

%disp(['waveform improvement over GTF   ',num2str(mean(delta6-delta4)),'dB'])


%plot(gapsms,snr_y(:,2),'--k')
h01 = plot([1,0],[1,1],'--','color',[1,1,1]*0.6);
h11 = plot([1,0],[1,1],'-','color',[1,1,1]*0.6);

h1=plot(gapsms,snr_y(:,3),'-','color',[0,0.6,0]);
h2=plot(gapsms,snr_y(:,4),'--b');
h3=plot(gapsms,snr_y(:,5),'-b');
h4=plot(gapsms,snr_y(:,6),'--r');
h5=plot(gapsms,snr_y(:,7),'-r');

for k=1:length(gaps_plot)
  plot(gapsms(gaps_plot(k)),snr_y(gaps_plot(k),5),'ob')
  plot(gapsms(gaps_plot(k)),snr_y(gaps_plot(k),7),'or')
end

%legend('NMF','tNMF','GTF','GTFtNMF','location','northeast')
legend([h1,h3,h5,h01,h11]',{'tNMF','GTF','GTFtNMF','unadapted filters','adapted filters'},'location','northeast')
delete(h01); delete(h11);

xlabel('missing region /ms','fontname',FontName,'FontSize',FSlg)
ylabel('SNR /dB','fontname',FontName,'FontSize',FSlg)                          
set(gca,'fontname',FontName,'FontSize',FSsm)
set(gca,'TickDir','out')
set(gca,'xlim',[min(gapsms)-1,max(gapsms)+1])
set(gca,'ylim',[min([delta2;delta3;delta4;delta5;delta6]),max([delta2;delta3;delta4;delta5;delta6])])


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% axis 12 -- pesq y

ax12 = axes('position',pos12); % produce axis
hold on

delta1= pesq_y(:,2);
delta2= pesq_y(:,3);
delta3 = pesq_y(:,4);
delta4 = pesq_y(:,5);
delta5 = pesq_y(:,6);
delta6 = pesq_y(:,7);


disp(['pesq improvement over GTF   ',num2str(mean(delta6-delta4)),'units'])
disp(['pesq improvement (GTF/lGTF)   ',num2str(mean(delta4-delta3)),'dB'])
disp(['pesq improvement (GTFtNMF/lGTFtNMF)   ',num2str(mean(delta6-delta5)),'dB'])

%disp(['pesq improvement over GTF   ',num2str(mean(delta6-delta4)),'units'])


for k=1:length(gaps_plot)
  plot(gapsms(gaps_plot(k)),delta4(gaps_plot(k)),'ob')
  plot(gapsms(gaps_plot(k)),delta6(gaps_plot(k)),'or')
end

%plot(gapsms,pesq_y(:,2),'--k')
plot(gapsms,pesq_y(:,3),'-','color',[0,0.6,0])
plot(gapsms,pesq_y(:,4),'--b')
plot(gapsms,pesq_y(:,5),'-b')
plot(gapsms,pesq_y(:,6),'--r')
plot(gapsms,pesq_y(:,7),'-r')

xlabel('missing region /ms','fontname',FontName,'FontSize',FSlg)
ylabel('PESQ','fontname',FontName,'FontSize',FSlg)               
set(gca,'fontname',FontName,'FontSize',FSsm)
set(gca,'TickDir','out')

set(gca,'xlim',[min(gapsms)-1,max(gapsms)+1])
set(gca,'ylim',[min([delta2;delta3;delta4;delta5;delta6]),max([delta2;delta3;delta4;delta5;delta6])])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% axis 13 --  amn
% delta1= snr_a_mn(:,2);
% delta2= snr_a_mn(:,3);
% delta3 = snr_a_mn(:,4);
% delta4 = snr_a_mn(:,5);


% ax13 = axes('position',pos13); % produce axis
% hold on
% plot(gapsms,snr_a_mn(:,2),'--k')
% plot(gapsms,snr_a_mn(:,3),'-k')
% plot(gapsms,snr_a_mn(:,4),'-b')
% plot(gapsms,snr_a_mn(:,5),'-r')


delta1= snr_loga_mn(:,2);
delta2= snr_loga_mn(:,3);
delta3 = snr_loga_mn(:,4);
delta4 = snr_loga_mn(:,5);
delta5 = snr_loga_mn(:,6);
delta6 = snr_loga_mn(:,7);

disp(['loga SNR improvement over tNMF   ',num2str(mean(delta6-delta2)),'dB'])
disp(['loga SNR improvement (GTF/lGTF)   ',num2str(mean(delta4-delta3)),'dB'])
disp(['loga SNR improvement (GTFtNMF/lGTFtNMF)   ',num2str(mean(delta6-delta5)),'dB'])

%disp(['loga SNR improvement over tNMF   ',num2str(mean(delta6-delta2)),'dB'])

ax13 = axes('position',pos13); % produce axis
hold on
%plot(gapsms,snr_loga_mn(:,2),'--k')
plot(gapsms,snr_loga_mn(:,3),'-','color',[0,0.6,0])
plot(gapsms,snr_loga_mn(:,4),'--b')
plot(gapsms,snr_loga_mn(:,5),'-b')
plot(gapsms,snr_loga_mn(:,6),'--r')
plot(gapsms,snr_loga_mn(:,7),'-r')

for k=1:length(gaps_plot)
  plot(gapsms(gaps_plot(k)),delta4(gaps_plot(k)),'ob')
  plot(gapsms(gaps_plot(k)),delta6(gaps_plot(k)),'or')
end

xlabel('missing region /ms','fontname',FontName,'FontSize',FSlg)
ylabel('SNR log-spec /dB','fontname',FontName,'FontSize',FSlg)               
set(gca,'fontname',FontName,'FontSize',FSsm)
set(gca,'TickDir','out')

set(gca,'xlim',[min(gapsms)-1,max(gapsms)+1])
set(gca,'ylim',[min([delta2;delta3;delta4;delta5;delta6]),max([delta2;delta3;delta4;delta5;delta6])])


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % axis 13 -- log amn
% delta1= snr_loga_mn(:,2);
% delta2= snr_loga_mn(:,3);
% delta3 = snr_loga_mn(:,4);
% delta4 = snr_loga_mn(:,5);


% ax13 = axes('position',pos13); % produce axis
% hold on
% plot(gapsms,snr_loga_mn(:,2),'--k')
% plot(gapsms,snr_loga_mn(:,3),'-k')
% plot(gapsms,snr_loga_mn(:,4),'-b')
% plot(gapsms,snr_loga_mn(:,5),'-r')

% for k=1:length(gaps_plot)
%   plot(gapsms(gaps_plot(k)),delta3(gaps_plot(k)),'ob')
%   plot(gapsms(gaps_plot(k)),delta4(gaps_plot(k)),'or')
% end

% xlabel('missing region /ms','fontname',FontName,'FontSize',FSlg)
% ylabel('SNR log-spec /dB','fontname',FontName,'FontSize',FSlg)               
% set(gca,'fontname',FontName,'FontSize',FSsm)
% set(gca,'TickDir','out')

% set(gca,'xlim',[min(gapsms)-1,max(gapsms)+1])
% set(gca,'ylim',[min([delta1;delta2;delta3;delta4]),max([delta1;delta2;delta3;delta4])])

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

for l=1:length(gapPos)
  xs = gapPos(l)-nsip_lim(1,1)+[-gaps(gaps_plot(1))/2,gaps(gaps_plot(1))/2];
  xs =  1000*[xs(1),xs(2),xs(2),xs(1),xs(1)]/fs;
  ys = [-100,-100,100,100,-100];
  pat=patch(xs,ys,missCol,'edgecolor','none');    
end

set(gca,'Layer','top') 
%plot(tim1,Ys(nsip_lim(1,1):nsip_lim(1,2),1,gaps_plot(1)),'-','color',colNoisy)
plot(tim1,yTest(nsip_lim(1,1):nsip_lim(1,2)),'-k','linewidth',1)
plot(tim1,Ys(nsip_lim(1,1):nsip_lim(1,2),5,gaps_plot(1)),'-','color',colGTF)
plot(tim1,Ys(nsip_lim(1,1):nsip_lim(1,2),7,gaps_plot(1)),'-','color',colGTFtNMF)

%xlabel('time /ms','fontname',FontName,'FontSize',FSlg)
ylabel('y_t','fontname',FontName,'FontSize',FSlg)               
set(gca,'fontname',FontName,'FontSize',FSsm)
set(gca,'TickDir','out')
set(gca,'xlim',[tim1(1),tim1(end)])

ylim = [-1,1]*1.01*max(abs(Ys(nsip_lim(1,1):nsip_lim(1,2),1,gaps_plot(1))));

set(gca,'ylim',ylim)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% axis 22 -- snippet 1, noise level 2
ax22 = axes('position',pos22); % produce axis
hold on

hold on
tim1 = 1000*(1:1+nsip_lim(1,2)-nsip_lim(1,1))/fs;
%plot(tim1,Ys(nsip_lim(1,1):nsip_lim(1,2),1,gaps_plot(2)),'-','color',colNoisy)

for l=1:length(gapPos)
  xs = gapPos(l)-nsip_lim(1,1)+[-gaps(gaps_plot(2))/2,gaps(gaps_plot(2))/2];
  xs =  1000*[xs(1),xs(2),xs(2),xs(1),xs(1)]/fs;
  ys = [-100,-100,100,100,-100];
  pat=patch(xs,ys,missCol,'edgecolor','none');    
end

set(gca,'Layer','top') 
plot(tim1,yTest(nsip_lim(1,1):nsip_lim(1,2)),'-k','linewidth',1)
plot(tim1,Ys(nsip_lim(1,1):nsip_lim(1,2),5,gaps_plot(2)),'-','color',colGTF)
plot(tim1,Ys(nsip_lim(1,1):nsip_lim(1,2),7,gaps_plot(2)),'-','color',colGTFtNMF)

%xlabel('time /ms','fontname',FontName,'FontSize',FSlg)
ylabel('y_t','fontname',FontName,'FontSize',FSlg)               
set(gca,'fontname',FontName,'FontSize',FSsm)
set(gca,'TickDir','out')
set(gca,'xlim',[tim1(1),tim1(end)])

ylim = [-1,1]*1.01*max(abs(Ys(nsip_lim(1,1):nsip_lim(1,2),1,gaps_plot(2))));

set(gca,'ylim',ylim)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Snippet 2

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% axis 31 -- snippet 2, noise level 1
ax31 = axes('position',pos31); % produce axis
hold on

hold on
tim1 = 1000*(1:1+nsip_lim(2,2)-nsip_lim(2,1))/fs;

for l=1:length(gapPos)
  xs = gapPos(l)-nsip_lim(2,1)+[-gaps(gaps_plot(1))/2,gaps(gaps_plot(1))/2];
  xs =  1000*[xs(1),xs(2),xs(2),xs(1),xs(1)]/fs;
  ys = [-100,-100,100,100,-100];
  pat=patch(xs,ys,missCol,'edgecolor','none');    
end

set(gca,'Layer','top') 
%plot(tim1,Ys(nsip_lim(2,1):nsip_lim(2,2),1,gaps_plot(1)),'-','color',colNoisy)
plot(tim1,yTest(nsip_lim(2,1):nsip_lim(2,2)),'-k','linewidth',1)
plot(tim1,Ys(nsip_lim(2,1):nsip_lim(2,2),5,gaps_plot(1)),'-','color',colGTF)
plot(tim1,Ys(nsip_lim(2,1):nsip_lim(2,2),7,gaps_plot(1)),'-','color',colGTFtNMF)

%xlabel('time /ms','fontname',FontName,'FontSize',FSlg)
ylabel('y_t','fontname',FontName,'FontSize',FSlg)               
set(gca,'fontname',FontName,'FontSize',FSsm)
set(gca,'TickDir','out')
set(gca,'xlim',[tim1(1),tim1(end)])

ylim = [-1,1]*1.01*max(abs(Ys(nsip_lim(2,1):nsip_lim(2,2),1,gaps_plot(1))));

set(gca,'ylim',ylim)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% axis 32 -- snippet 2, noise level 2
ax32 = axes('position',pos32); % produce axis
hold on

hold on
tim1 = 1000*(1:1+nsip_lim(2,2)-nsip_lim(2,1))/fs;

for l=1:length(gapPos)
  xs = gapPos(l)-nsip_lim(2,1)+[-gaps(gaps_plot(2))/2,gaps(gaps_plot(2))/2];
  xs =  1000*[xs(1),xs(2),xs(2),xs(1),xs(1)]/fs;
  ys = [-100,-100,100,100,-100];
  pat=patch(xs,ys,missCol,'edgecolor','none');    
end

set(gca,'Layer','top') 
%plot(tim1,Ys(nsip_lim(2,1):nsip_lim(2,2),1,gaps_plot(2)),'-','color',colNoisy)
plot(tim1,yTest(nsip_lim(2,1):nsip_lim(2,2)),'-k','linewidth',1)
plot(tim1,Ys(nsip_lim(2,1):nsip_lim(2,2),5,gaps_plot(2)),'-','color',colGTF)
plot(tim1,Ys(nsip_lim(2,1):nsip_lim(2,2),7,gaps_plot(2)),'-','color',colGTFtNMF)

%xlabel('time /ms','fontname',FontName,'FontSize',FSlg)
ylabel('y_t','fontname',FontName,'FontSize',FSlg)               
set(gca,'fontname',FontName,'FontSize',FSsm)
set(gca,'TickDir','out')
set(gca,'xlim',[tim1(1),tim1(end)])

ylim = [-1,1]*1.01*max(abs(Ys(nsip_lim(2,1):nsip_lim(2,2),1,gaps_plot(2))));

set(gca,'ylim',ylim)



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Snippet 3

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% axis 41 -- snippet 3, noise level 1
ax41 = axes('position',pos41); % produce axis
hold on

hold on
tim1 = 1000*(1:1+nsip_lim(3,2)-nsip_lim(3,1))/fs;

for l=1:length(gapPos)
  xs = gapPos(l)-nsip_lim(3,1)+[-gaps(gaps_plot(1))/2,gaps(gaps_plot(1))/2];
  xs =  1000*[xs(1),xs(2),xs(2),xs(1),xs(1)]/fs;
  ys = [-100,-100,100,100,-100];
  pat=patch(xs,ys,missCol,'edgecolor','none');    
end

set(gca,'Layer','top') 
%plot(tim1,Ys(nsip_lim(3,1):nsip_lim(3,2),1,gaps_plot(1)),'-','color',colNoisy)
plot(tim1,yTest(nsip_lim(3,1):nsip_lim(3,2)),'-k','linewidth',1)
plot(tim1,Ys(nsip_lim(3,1):nsip_lim(3,2),5,gaps_plot(1)),'-','color',colGTF)
plot(tim1,Ys(nsip_lim(3,1):nsip_lim(3,2),7,gaps_plot(1)),'-','color',colGTFtNMF)

xlabel('time /ms','fontname',FontName,'FontSize',FSlg)
ylabel('y_t','fontname',FontName,'FontSize',FSlg)               
set(gca,'fontname',FontName,'FontSize',FSsm)
set(gca,'TickDir','out')
set(gca,'xlim',[tim1(1),tim1(end)])

ylim = [-1,1]*1.01*max(abs(Ys(nsip_lim(3,1):nsip_lim(3,2),1,gaps_plot(1))));

set(gca,'ylim',ylim)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% axis 42 -- snippet 1, noise level 2
ax42 = axes('position',pos42); % produce axis
hold on

hold on
tim1 = 1000*(1:1+nsip_lim(3,2)-nsip_lim(3,1))/fs;

for l=1:length(gapPos)
  xs = gapPos(l)-nsip_lim(3,1)+[-gaps(gaps_plot(2))/2,gaps(gaps_plot(2))/2];
  xs =  1000*[xs(1),xs(2),xs(2),xs(1),xs(1)]/fs;
  ys = [-100,-100,100,100,-100];
  pat=patch(xs,ys,missCol,'edgecolor','none');    
end

set(gca,'Layer','top') 
%plot(tim1,Ys(nsip_lim(3,1):nsip_lim(3,2),1,gaps_plot(2)),'-','color',colNoisy)
plot(tim1,yTest(nsip_lim(3,1):nsip_lim(3,2)),'-k','linewidth',1)
plot(tim1,Ys(nsip_lim(3,1):nsip_lim(3,2),5,gaps_plot(2)),'-','color',colGTF)
plot(tim1,Ys(nsip_lim(3,1):nsip_lim(3,2),7,gaps_plot(2)),'-','color',colGTFtNMF)

xlabel('time /ms','fontname',FontName,'FontSize',FSlg)
ylabel('y_t','fontname',FontName,'FontSize',FSlg)               
set(gca,'fontname',FontName,'FontSize',FSsm)
set(gca,'TickDir','out')
set(gca,'xlim',[tim1(1),tim1(end)])

ylim = [-1,1]*1.01*max(abs(Ys(nsip_lim(3,1):nsip_lim(3,2),1,gaps_plot(2))));

set(gca,'ylim',ylim)

  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Printing figure

if printOutput==1
  print(figure1,'-depsc','-painters',[saveDir,saveName],'-loose');
  %plot2svg([filename,'ColBar','.svg'],figure2,'png')
  %str = ['! gv ',filename,'ColBar','.eps &'];
  %eval(str)
end
