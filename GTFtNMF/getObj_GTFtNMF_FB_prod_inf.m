function [Obj,dObj] = getObj_GTFtNMF_FB_prod_inf(z,W,y,lamv,varv,omv, ...
						   vary,lenx,mux,varx,preZ,tol);

% function [Obj,dObj] =
% getObj_GTFtNMF_FB_prod_inf(z,W,y,lamv,varv,omv,vary,lenx,mux,varx,preZ,tol);
% 
% GTF-tNMF with product filter bank formulation for optimising
% over the envelopes. I.e. for performing inference only
%
% y_t = sum_d Re(a_{d,t }v_{d,t}) + sigy*\eta_t
% v_{d,t} = lamv_d*exp(i omv_d)*v_{d,t-1} +  \sigma_d \epsilon_{d,t}
% amp_{d,t} = sqrt(H*W/2)
% W_{k,d} = exp(logW_{k,d})/sum_d exp(logW_{k,d})
%
% Each of the temporal basis functions in NMF is given by: 
% H_{t,k} = exp(X_{t,k})  where X{1:T,k} ~ mux_k + GP(lenx_k,varx_k)
%
% The basis function representation of the GP used for the temporal
% basis functions is employed rather than the function-space view
% so that:
%
% bas_k \propto sqrt(varx_k)*exp(-1/(lenx_k^2) t^2)
% X_{1:T,k} = convolve(Z_{1:K,k},bas_k) 
% 
%
% INPUTS
% z = vectorised matrices of temporal basis functions [T*K,1]
% y = signal [T,1]
% W = spectral weights [K,D]
% lamv = dynamical AR parameters [D,1]
% varv = dynamical noise parameters [D,1]
% omv = mean frequencies of the sinusoids [D,1]
% vary = observation noise [T,1]
% lenx = squared exponential length-scale [1,K]
% mux = steady state mean of the log temporal priors [1,K]
% varx = steady state variance of the log temporal priors [1,K]
% preZ = precision of Z
% tol = tolerance parameter which sets the extent of the Gaussian
%       basis functions
%
% OUTPUTS
% Obj = objective
% dObjdX2 = derivative of the objective by  z, [T*K,1] 
%
% This function is too long - in the future it would be better to
% split it up into smaller modules.

  % read out the dimensions of the model  
  D = length(lamv);
  T = length(y);  
  K = length(lenx);
  
  if length(vary)==1
      vary = vary*ones(T,1);
  end
  
  % form the temporal basis functions
  Z = reshape(z(1:T*K),[T,K]);

  for k=1:K
    tau = ceil(lenx(k)/sqrt(2)*tol);
    bh = sqrt(varx(k)*sqrt(2)/(lenx(k)*sqrt(pi))); % basis function height
                                                   % to get correct variance
    bas{k} = bh*exp(-1/(lenx(k)^2)*([-tau:1:tau]').^2);
  end
  
  H = zeros(T,K);
  for k=1:K
     H(:,k) =  exp(conv(Z(:,k),bas{k},'same')+mux(k));
  end
 
  % form the envelopes
  A = (1/2*H*W).^(1/2);

  %%%%%%%%%%%%%%%%%%%%%%%%
  % Likelihood component

  % Calculate the sufficient statistics of the
  % posterior distribution over the carriers p(X1|Amp,Y)

  [ObjA,Xfin,Pfin] = kalman_GTFtNMF_FB_prod(y,A',lamv,varv,omv,vary);
  
  % extract the necessary statistics
  % it would be nice to speed up the following

  X1 = zeros(T,D);
  X1X1 = zeros(T,D,D);
    
  for d=1:D
    e = 2*d-1;
    X1(:,d) =  A(:,d).*reshape(Xfin(1,e,:),[T,1]); 
  end
 
  for d=1:D
    for dd=1:D
      e = 2*d-1; 
      ee = 2*dd-1;
      X1X1(:,d,dd) = X1(:,d).*X1(:,dd)+A(:,d).*A(:,dd).*reshape(Pfin(e,ee,:),[T,1,1]);
    end
  end
  
%  dAdX = 1./(1+exp(-X2(1:T,:)*G'-ones(T,1)*Mu2'));

  dY = 1./(4*vary*ones(1,D)).*(sum(X1X1,3)-(y*ones(1,D)).*X1)./A.^2; %*W).*H;
  dX = (dY*W').*H;
  
  dObjAdZ = zeros(T,K);
  for k=1:K
    dObjAdZ(:,k) = conv(dX(:,k),bas{k},'same');
  end
  
  %%%%%%%%%%%%%%%%%%%%%%%%%
  % Prior component

  ObjB = 1/2*sum(Z(:).^2)*preZ;
  dObjBdZ = Z*preZ;

%%%%%%%%%%%%%%%%%%%%

Obj = (-ObjA+ObjB)/T;
dObj = (dObjAdZ(:)+dObjBdZ(:))/T;

% Obj = ObjB;
% dObj = [dObjBdZ(:);zeros(D*K,1)];

%Obj = -ObjA;
%dObj = [dObjAdZ(:);dObj1dW(:)];