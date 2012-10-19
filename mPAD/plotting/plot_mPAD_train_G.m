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
  
  for k = 1:K
    kPlot = 1+2*(k-1);
    subplot(K,2,kPlot)
    hold on
    
    if sum(G(:,k))<1
      sign = -1;
    else
      sign = 1;
    end
    
    plot(sign*X2old(1:T,k),'-k','linewidth',2)
    plot(sign*X2(1:T,k),'-r','linewidth',1)
    legend('old','new')
    ylabel(['X2_',num2str(k)])
    xlabel('time')
  end
  
  for d=1:D
    dPlot = 2*d;
    subplot(D,2,dPlot)
    hold on
    Ychan = real(X1).*A;
    plot(Ychan(:,d),'-k')
    plot(A(:,d),'-r','linewidth',2)
    ylabel(['channel_',num2str(d)])
    xlabel('time')

  end
  
  drawnow
