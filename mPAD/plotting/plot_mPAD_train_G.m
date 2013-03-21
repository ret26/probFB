  function plot_mPAD_train_G(y,X2,X2old,Params)  
  
%  function plot_mPAD_train_G(y,X2,X2old,Params)  
   
  [G,Lam1,Var1,Len2,Var2,Mu2,vary,om] = unpackParamsMPAD(Params);

  T = size(y,1);
  Tx = size(X2,1);
  
  A = log(1+exp(X2(1:T,:)*G'+ones(T,1)*Mu2'));

  [ObjA,Xfin,Pfin] = kalman_mPAD_FB(Params,y,A');
  
  [X1,covX1] = getFBLDSOutput(Xfin,Pfin);
  X1 = X1';

  [D,K] = size(G);

  Params.G
  
  clf;

  top = 0.01;
  left = 0.1;
  right = 0.01;
  bottom = 0.1;
  vspace = 0.01;
  hspace = 0.1;
  wspace = 0.2;
  
  
  height1 = (1-top-bottom-vspace*(K-1))/K;
  height2 = (1-top-bottom-vspace*(D-1))/D;
  width = (1-left-right-2*hspace-wspace)/2;
  posCur1 = [left+wspace+hspace,1-top-height1,width,height1];
  posCur2 = [left+wspace+width+2*hspace,1-top-height2,width,height2];
  down1 = [0,-height1-vspace,0,0];
  down2 = [0,-height2-vspace,0,0];
  posW = [left,1/2-wspace/2,wspace,wspace];
  
  axW = axes('position',posW);
  imagesc(G)
  
  for k = 1:K
    %kPlot = 1+2*(k-1);
    %subplot(K,2,kPlot)
    ax1k = axes('position',posCur1);
    posCur1 = posCur1+down1;
    hold on   
    
    if sum(G(:,k))<1
      sign = -1;
    else
      sign = 1;
    end

    plot(sign*X2old(1:T,k),'-k','linewidth',2)
    plot(sign*X2(1:T,k),'-r','linewidth',1)
    set(gca,'xlim',[1,T],'ylim',[min(sign*X2old(1:T,k)),max(sign*X2old(1:T,k))+1e-9])
    
    ylabel(['X2_',num2str(k)])

    if  k==K
      legend('old','new')
      xlabel('time')
    else
      set(gca,'xticklabel','')
    end
    
  end
  
  
  
    
  for d=1:D
    %dPlot = 2*d;
    %subplot(D,2,dPlot)

    ax2k = axes('position',posCur2);
    posCur2 = posCur2+down2;

    hold on
    Ychan = real(X1).*A;
    plot(Ychan(:,d),'-k')
    plot(A(:,d),'-r','linewidth',2)
    ylabel(['channel_',num2str(d)])
    set(gca,'xlim',[1,T],'ylim',[min(Ychan(:,d)),max(Ychan(:,d))+1e-7])
    
    if d==D
      xlabel('time')
    else
      set(gca,'xticklabel','')
    end
  end
  
  drawnow
