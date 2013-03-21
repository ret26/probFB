function [Obj,dObj] = getObj_mPAD_fixG_basis_func(z2G,y,Params,GScale,tol);

% function [Obj,dObj] = getObj_mPAD_fixG_basis_func(z2G,y,Params,GScale,tol);
%
% Gets the mPAD objective function and derivatives wrt z and
% G. Uses the basis function representation of the transformed
% envelopes rather than the function-space view. 
%
% INPUTS
% z2G = Transformed envelopes basis function coefficients [T*K,1] (gets reshaped
%      into [T,K]) and G [D*K,1] (gets reshaped to [D,K])
% y = Obsevations [T,1] 
% Params = structure containing the following parameters
%   Lam1 = Carrier dynamics [D,2]
%   Var1 = Carrier noise [1,D]
%   Len2 = Modulator Length Scales [K,1] (not required)
%   Var2 = Modulator Noise [1,K]
%   vary = observation noise  
% GScale = fixed norm of the G weight matrix [K,1]
% tol = tolerance parameter which sets the extent of the Gaussian
%       basis functions
%
% OUTPUTS
% Obj = objective
% dObjdX2 = derivative of the objective by  Z2 and G, [T*K,1] 
%
% This function is too long - in the future it would be better to
% split it up into smaller modules.

  % read out the dimensions of the model
  D = length(Params.Var1);
  T = length(y);  
  K = length(Params.Len2);

  % readout the basis function coefficients
  Z2 = reshape(z2G(1:T*K),[T,K]);
  
  % Rescale G and read out other parameters
  G = reshape(z2G(T*K+1:T*K+K*D),[D,K]);
  scale = diag(sqrt(diag(G'*G)));
  G = (G/scale)*diag(GScale);

  % this is messy and dangerous - the version of G stored in the
  % parameter structure is old until this point
  Params.G = G;  
  [G,Lam1,Var1,Len2,Var2,Mu2,vary] = unpackParamsMPAD(Params);

  %%%%%%
  % form envelopes from basis functions

  for k=1:K
    tau = ceil(Len2(k)/sqrt(2)*tol);
    bh = sqrt(Var2(k)*sqrt(2)/(Len2(k)*sqrt(pi))); % basis function height
                                           % to get correct variance
    bas{k} = bh*exp(-1/(Len2(k)^2)*([-tau:1:tau]').^2);
  end
  
  X2 = zeros(T,K);
  for k=1:K
    X2(:,k) =  conv(Z2(:,k),bas{k},'same');
  end
 
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
  dX = ((X1X1A-(y*ones(1,D)).*X1).*dAdX2)/vary;
  dZ = dX*G;
  
  dObjAdZ2 = zeros(T,K);
  for k=1:K
    dObjAdZ2(:,k) = conv(dZ(:,k),bas{k},'same');
  end
  
%  dObjAdGhat = tensor(dX,[-1,2],X2(1:T,:),[-1,1])';  
  dObjAdGhat = dX'*X2(1:T,:);
  dObjAdG = zeros(D,K);

  for k=1:K
    dObjAdG(:,k) = (GScale(k)*eye(D)-G(:,k)*G(:,k)'/GScale(k))/scale(k,k)*dObjAdGhat(:,k);
  end

  %%%%%%%%%%%%%%%%%%%%%%%%%
  % Prior component

  ObjB = 1/2*sum(Z2(:).^2);
  dObjBdZ2 = Z2;


%%%%%%%%%%%%%%%%%%%%

Obj = (-ObjA+ObjB)/T;
dObj = ([dObjAdZ2(:)+dObjBdZ2(:);dObjAdG(:)])/T;

% Obj = ObjB;
% dObj = [dObjBdZ2(:);zeros(D*K,1)];

%Obj = -ObjA;
%dObj = [dObjAdZ2(:);dObjAdG(:)];