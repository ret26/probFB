% plot parameters of the model trained on a speech sound for the
% TSP paper

clear;

printOutput = 1;

loadDir = '/home/rich/Data/probFB/nmf/';
%simName = 'trainDenoise9_wide_D40_med.mat';
%simName = 'trainDenoise10_D40.mat';
simName = 'trainDenoise11_D40.mat';
load([loadDir,simName])

saveDir = '/home/rich/Synchronised/Writings/mpad/tsp/figs/';
saveName = 'parameters_fig_train.eps'


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% rows are features

meanom = WEst4*om4;
[val,ind] = sort(meanom);
WEst4_sort = WEst4(ind,:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure1=figure;
PP = [0,0,7.70,8.5]; %*** paper position in centimeters
PS = PP(end-1:end); % paper size in centimeters

set(figure1,'paperpositionmode','manual','paperposition', ...
        PP,'papersize',PS, 'paperunits','centimeters');

% So the figure is the same size on the screen as when it is printed:
pu = get(gcf,'PaperUnits');
pp = get(gcf,'PaperPosition');
set(gcf,'Units',pu,'Position',pp)

PR = PS(1)/PS(2);

colormap gray
colormap(flipud(colormap))
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

freqLim = [100,fs/2];
D = length(Var14);
NumPoints = 4000;
freqs = logspace(log10(0.001),log10(0.5),NumPoints);

left = 0.11;
right = 0.03;
top = 0.07;
bottom = 0.10;
vspaceSM =  0.01;
hspace = 0.11;
vspace = 0.15;

w1 = 0.45;
h1 = 0.45;

h2 = 1-top-bottom-h1-vspace;
w2 = 1-left-right;

pos11 = [left,1-top-h1,w1,h1];
pos21 = [left,bottom,w2,h2];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
indPlot = [7,12,23,37];%[7,16,25,33];
M = length(indPlot);

hk = (1-bottom-h1-(M-1)*vspaceSM)/M;
width = 1-w1-left-right-hspace;

posk = [left+w1+hspace,1-top-h1,width,hk];
up = [0,hk+vspaceSM,0,0];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fonts
FontName = 'Times';
FSvsm = 5;
FSsm = 6;
FSmed = 7;
FSlg = 8;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ax11 = axes('position',pos11);
hold on
imagesc(log(WEst4_sort+1e-4))

box on;
set(gca,'xlim',[1/2,D+1/2],'ylim',[+1/2,D+1/2])
xlabel('filter index','fontname',FontName,'FontSize',FSlg)
ylabel('feature index','fontname',FontName,'FontSize',FSlg)
set(gca,'TickDir','out')
set(gca,'fontname',FontName,'FontSize',FSsm)
xtick = [1,10,20,30,40]';
set(gca,'ytick',xtick,'yticklabel',num2str(xtick))
set(gca,'xtick',xtick,'xticklabel',num2str(xtick))
colorbar('location','NorthOutside','fontname',FontName,'FontSize',FSvsm);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ax21 = axes('position',pos21);
hold on
spec = get_pSTFT_spec(freqs,Lam14,Var14,om4);
sum_spec = sum(spec);
sum_spec = max(spec(:))*sum_spec/max(sum_spec);
ah1= plot(freqs*fs, sum_spec,'-k','linewidth',1.5);
set(gca,'TickDir','out')
for d=1:D,
  ah2=plot(freqs*fs,spec(d,:));
end
set(gca,'fontname',FontName,'FontSize',FSsm)

hal=legend([ah2,ah1],{'log-filter responses','summed response'},'location','northeast','fontname',FontName,'FontSize',FSsm,'boxoff');

%posleg = get(hal,'position');
%poslegnew = posleg + [0.2,0,-0.2,0];
%set(hal,'position',poslegnew);
%keyboard

% if logSpec==1
%   set(gca,'yscale','log')
% end
set(gca,'xscale','log')
set(gca,'xlim',freqLim)
%set(gca,'ylim',[0,max(spec(:))*1.01])

set(gca,'yscale','log')
spec_ylim = [max(spec(:))/1000,max(spec(:))*1.1];
set(gca,'ylim',spec_ylim)

%set(gca,'visible','off')
xtick = [125,500,2000,8000];

set(gca,'xtick',xtick,'xticklabel',num2str(xtick'))
set(gca,'yticklabel','')
xlabel('frequency /Hz','fontname',FontName,'FontSize',FSlg)
ylabel('log-filter responses','fontname',FontName,'FontSize',FSlg)
set(gca,'TickDir','out')

for m=1:M
  
  vark_eff = WEst4_sort(indPlot(m),:)'.^2.*Var14;
  spec = get_pSTFT_spec(freqs,Lam14,vark_eff,om4);
    
  axk = axes('position',posk);
  hold on

  sum_spec = sum(spec);
  sum_spec = max(spec(:))*sum_spec/max(sum_spec);

  plot(freqs*fs,sum_spec,'-k','linewidth',1.5)

  for d=1:D,
    plot(freqs*fs,spec(d,:))
  end
  
set(gca,'fontname',FontName,'FontSize',FSsm)
  
  % if logSpec==1
  %   set(gca,'yscale','log')
  % end
 set(gca,'TickDir','out')
  set(gca,'xscale','log')
  set(gca,'xlim',freqLim)
%  set(gca,'visible','off')
  posk = posk+up;
%  set(gca,'xticklabel','','yticklabel','')

  ylabel(['feature ',num2str(indPlot(m))],'fontname',FontName,'FontSize',FSvsm)

  set(gca,'yscale','log')
spec_ylim_cur =  [max(sum_spec)/100,max(sum_spec)*1.1];
  set(gca,'ylim',spec_ylim_cur)

 if m==1
   xlabel('frequency /Hz','fontname',FontName,'FontSize',FSlg)
   set(gca,'yticklabel','')
   xtick = [125,500,2000,8000];
   set(gca,'xtick',xtick,'xticklabel',num2str(xtick'))
 else
   set(gca,'xticklabel','','yticklabel','')
 end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Printing figure

if printOutput==1
  print(figure1,'-depsc','-painters',[saveDir,saveName],'-loose');
  %plot2svg([filename,'ColBar','.svg'],figure2,'png')
  %str = ['! gv ',filename,'ColBar','.eps &'];
  %eval(str)

end
