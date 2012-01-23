function figH = plotYX2(y,X2,tit,FSamp)

% function figH = plotYX2(y,X2,tit,FSamp)
%
% Plots a waveform and the extracted X2 values
%
% INPUTS
% y = waveform size [T,1]
% X2 = latent modulation variables [T,K]
% tit = title of the figure (string)
% FSamp = sample rate
%
% OUTPUTS
% figH = figure handle

left = .1;
right = .03;
top = .1;
bottom = .1;
vspace = .05;

[T,K] = size(X2);

height = (1-top-bottom-vspace*K)/(K+1);
width = 1-left-right;

pos(1,:) = [left,1-top-height,width,height];
down = -[0,height+vspace,0,0];

time = (0:T-1)/FSamp;
tLim = [min(time),max(time)];

for k=1:K
  pos(k+1,:) = pos(k,:)+down;
end

figH = figure;

ax1 = axes('position',pos(1,:));
plot(time,y,'-k')
set(gca,'xlim',tLim,'ylim',[-1.01,1.01]*max(abs(y)))
ylabel('y')
set(gca,'xticklabel','')
title(tit)

for k=1:K
  axk = axes('position',pos(1+k,:));
  plot(time,X2(:,k),'-k')
  set(gca,'xlim',tLim,'ylim',[-1.01,1.01]*max(abs(X2(:,k))))
  ylabel(['X2_',num2str(k)])

  if k<K
    set(gca,'xticklabel','')
  else
    xlabel('time')
  end
end
