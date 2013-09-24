function [obj,dobj] = getObj_GTFtNMF_FB_hier_EM_inf(z,W,VV,RVV,lamv,varv,lenx,mux,varx,preZ,tol);

% function [obj,dobj] = getObj_GTFtNMF_FB_hier_EM_inf(z,W,VV,RVV,lamv,varv,lenx,mux,varx,tol);
% 
% EM inference of GTF-tNMF with hierarchical filter bank formulation
% for optimising both over the envelopes and the weights. I.e. for
% performing inference only.
%
% y_t = sum_d Re(v_{d,t}) + sigy*\eta_t
% v_{d,t} = amp_{d,t}/amp_{d,t-1}*lamv_d*exp(i omv_d)*v_{d,t-1} 
%                              + amp_{d,t} \sigma_d \epsilon_{d,t}
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
% z = vectorised matrices of temporal basis functions and W [T*K+D*K,1]
% W = spectral tNMF weights [K,D]
% VV = marginal sufficient statistics of the filter bank coefficients
%      <v_{d,t}^T*v_{d,t}> [T,D]
% RVV = lag-1 marginal sufficient statistics of the filter bank
%       coefficients <v_{d,t}^T*R*v_{d,t-1}> [T-1,D]
% lamv = dynamical AR parameters [D,1]
% varv = dynamical noise parameters [D,1]
% lenx = squared exponential length-scale [1,K]
% mux = steady state mean of the log temporal priors [1,K]
% varx = steady state variance of the log temporal priors [1,K]
% preZ = precision of Zs
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
  [T,D] = size(VV);
  K = length(lenx);
  
  % form the temporal basis functions
  Z = reshape(z,[T,K]);

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
  dyn2 = (1+lamv(:).^2)./(2*varv(:));
  var0 = varv(:)./(1-lamv(:).^2);
  dyn3 = lamv(:).^2./(2*varv(:))+1./(2*var0);
  dyn4 = 1./(2*varv(:));

  obj1 = 2*sum(log(A(:)))+sum(VV(2:T-1,:)./(A(2:T-1,:).^2),1)*dyn2+...
	 (VV(T,:)./A(T,:).^2)*dyn4+(VV(1,:)./A(1,:).^2)*dyn3 ...
	 -2*sum(RVV./(A(2:T,:).*A(1:T-1,:)),1)*dyn4;
  
  dobj1dA = 2./A;
  
  dobj1dA(1,:) = dobj1dA(1,:)-2*VV(1,:)./(A(1,:).^3)*diag(dyn3) ...
                 +2*RVV(1,:)./(A(1,:).^2.*A(2,:))*diag(dyn4);

  dobj1dA(2:T-1,:) =  dobj1dA(2:T-1,:) ...
              -2*VV(2:T-1,:)./(A(2:T-1,:).^3)*diag(dyn2) + ...
               2*((RVV(1:T-2,:)./A(1:T-2,:)+RVV(2:T-1,:)./A(3:T,:)) ...
		  ./A(2:T-1,:).^2)*diag(dyn4);

  dobj1dA(T,:) = dobj1dA(T,:)-2*VV(T,:)./(A(T,:).^3)*diag(dyn4) ...
         +2*RVV(T-1,:)./(A(T,:).^2.*A(T-1,:))*diag(dyn4);
    
  
  dY = dobj1dA./(4*A);
  dX = (dY*W').*H;
  
  dobj1dZ = zeros(T,K);
  for k=1:K
    dobj1dZ(:,k) = conv(dX(:,k),bas{k},'same');
  end
  
  %%%%%%%%%%%%%%%%%%%%%%%%%
  % Prior component

  obj2 = 1/2*sum(Z(:).^2)*preZ;
  dobj2dZ = Z*preZ;

  %%%%%%%%%%%%%%%%%%%%

  obj = (obj1+obj2)/T;
  dobj = ([dobj1dZ(:)+dobj2dZ(:)])/T;
  
