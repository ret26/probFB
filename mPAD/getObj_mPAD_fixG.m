function [Obj,dObj] = getObj_mPAD_fixG(X2G,y,fftCov,Params,GScale);

% function [Obj,dObj] = getObj_mPAD_fixG(X2G,y,fftCov,Params,GScale);
%
% Gets the mPAD objective function and derivatives wrt X2
% and G.
%
% INPUTS
% X2G = Transformed envelopes [Tx*K,1] (gets reshaped
%      into [Tx,K]) and G [D*K,1](gets reshaped to [D,K])
% y = Obsevations [T,1] 
% fftCov = FFTs of the covariance matrices, [Tx,K]
% Params = structure containing the following parameters
%   Lam1 = Carrier dynamics [D,2]
%   Var1 = Carrier noise [1,D]
%   Len2 = Modulator Length Scales [K,1] (not required)
%   Var2 = Modulator Noise [1,K]
%   vary = observation noise  
%
% OUTPUTS
% Obj = objective
% dObjdX2 = derivative of the objective by X2 and G, [2(T-1)*K,1] 
%
% This function is too long - in the future it would be better to
% split it up into smaller modules.

  D = length(Params.Var1);
  [Tx,K] = size(fftCov);  
  T = length(y);

  G = reshape(X2G(Tx*K+1:Tx*K+K*D),[D,K]);

  % Rescale G
  scale = diag(sqrt(diag(G'*G)));
  G = (G/scale)*diag(GScale);
  Params.G = G;  

  [G,Lam1,Var1,Len2,Var2,Mu2,vary] = unpackParamsMPAD(Params);

  X2 = reshape(X2G(1:Tx*K,1),[Tx,K]);    

  A = log(1+exp(X2(1:T,:)*G'+ones(T,1)*Mu2'));
  
  dAdX2 = 1./(1+exp(-X2(1:T,:)*G'-ones(T,1)*Mu2'));

  %%%%%%%%%%%%%%%%%%%%%%%%
  % Likelihood component

  % Calculate the sufficient statistics of the
  % posterior distribution over the carriers p(X1|Amp,Y)

  [ObjA,Xfin,Pfin] = kalman_mPAD_FB(Params,y,A');
  
  % extract the necessary statistics
  % it would be nice to speed up the following

  X1 = zeros(T,D);
  X1X1 = zeros(T,D,D);
    
  for d=1:D
    e = 2*d-1;
    X1(:,d) =  reshape(Xfin(1,e,:),[T,1]); 
  end
 
  for d=1:D
    for dd=1:D
      e = 2*d-1; 
      ee = 2*dd-1;
      X1X1(:,d,dd) = X1(:,d).*X1(:,dd)+reshape(Pfin(e,ee,:),[T,1,1]);
    end
  end
  
    % The following is a costly step - better to loop over D
  % X1X1A = tensor(X1X1,[1,-1,2],A,[1,-1]); % replaced by...

  X1X1A = zeros(T,D);

  for d=1:D
    X1X1A(:,d) = sum(squeeze(X1X1(:,:,d)).*A,2);
  end

%  dX = 1/vary*(X1X1A-tensor(y,1,X1,[1,2])).*dAdX2;
  dX = 1/vary*(X1X1A-(y*ones(1,D)).*X1).*dAdX2;
  dObjAdX2 = [dX*G;zeros(Tx-T,K)];
  
%  dObjAdGhat = tensor(dX,[-1,2],X2(1:T,:),[-1,1])';  
  dObjAdGhat = dX'*X2(1:T,:);
  
  dObjAdG = zeros(D,K);

  for k=1:K
    dObjAdG(:,k) = (GScale(k)*eye(D)-G(:,k)*G(:,k)'/GScale(k))/scale(k,k)*dObjAdGhat(:,k);
  end

  %%%%%%%%%%%%%%%%%%%%%%%%%
  % Prior component
  
  ObjB = 0;
  dObjBdX2 = zeros(Tx,K);
  
  % Could vectorise the following over k using tensor, but
  % I don't think this is the limiting step
  
tiny = 1e-6;

for k=1:K
  fftx2k = fft(X2(:,k));
  InvCovX = ifft(fftx2k./(fftCov(:,k)+tiny));
  ObjB = ObjB+1/2*InvCovX'*(X2(:,k))/Var2(k);
  dObjBdX2(:,k) = InvCovX/Var2(k);
end

%%%%%%%%%%%%%%%%%%%%

Obj = (-ObjA+ObjB)/T;

dObj = ([dObjAdX2(:)+dObjBdX2(:);dObjAdG(:)])/T;

%Obj = ObjB;
%dObj = [dObjBdX2(:);zeros(D*K,1)];

% Obj = -ObjA;
% dObj = [dObjAdX2(:);dObjAdG(:)];