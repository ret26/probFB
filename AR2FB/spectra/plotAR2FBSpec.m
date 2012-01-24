function figH = plotAR2FBSpec(Lam,Var,varargin)
  
  % function figH = plotAR2FBSpec(Lam,Var)
  % or 
  % function figH = plotAR2FBSpec(Lam,Var,FSamp,logSpec)

  % Plots the spectra of a filterbank of AR(2)
  % processes
  
  % INPUTS
  % Lam = AR dynamics [D,2]
  % Var = AR noise variances [D,1]
  % Optional Inputs:
  % FSamp = sample rate  
  % logSpec = 1 if log-spectra to be plotted
  % 
  % OUTPUTS
  % figH = figure handle
  
    
  if nargin>2
    FSamp = varargin{1};
  else
    FSamp=1;
  end

  if nargin>3
    logSpec = varargin{2};
  else
    logSpec=0;
  end

  D = length(Var);
  NumPoints = 4000;
  Spec = zeros(D,NumPoints);
  
  for d=1:D  
    [Freqs,Spec(d,:),fMAX(d),SpecMAX(d),dF1(d),dF2(d)] ... 
        = getSpecAR2(Lam(d,:),Var(d),NumPoints,[0,0.5]);
  end
  
  figH = figure;
  hold on;
  
  for d=1:D,
    plot(Freqs*FSamp,Spec(d,:))
    plot(FSamp*fMAX(d)*[1,1],[0,SpecMAX(d)],'-g')
    plot([dF1(d),dF2(d)]*FSamp,SpecMAX(d)/2*[1,1],'-g')

  end
  
  plot(Freqs*FSamp,max(SpecMAX)*sum(Spec)/max(sum(Spec)),'-m','linewidth',3)
 
  if logSpec==1
    set(gca,'yscale','log')
  end
 %set(gca,'xscale','log')