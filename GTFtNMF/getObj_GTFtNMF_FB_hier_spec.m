function [Obj,dObj] = getObj_GTFtNMF_FB_hier_spec(zlogwspec,y, ...
						   vary,lenx,mux,varx,preZ,tol);

% function [Obj,dObj] =
% getObj_GTFtNMF_FB_hier(zlogw,y,lamv,varv,omv,vary,lenx,mux,varx,preZ,tol);
% 
% GTF-tNMF with hierarchical filter bank formulation for optimising
% both over the envelopes and the weights. I.e. for performing
% inference and learning.
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
% zlogWspec = vectorised matrices of temporal basis functions, W and
%         spectral parameters [T*K+D*K+3*D,1]
% y = signal [T,1]
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
% dObjdX2 = derivative of the objective by  z and logW, [T*K+D*K,1] 
%
% This function is too long - in the future it would be better to
% split it up into smaller modules.



  % read out the dimensions of the model  
  T = length(y);  
  K = length(lenx);
  D = (length(zlogwspec)-T*K)/(3+K);

  % form the normalised spectral weights
  W = reshape(exp(zlogwspec(T*K+1:T*K+K*D)),[K,D]);
  W = diag(1./sum(W,2))*W; %normalised weights
  
  % extract spectral parameters which we note are:
  %  
  % lamv = dynamical AR parameters [D,1]
  % varv = dynamical noise parameters [D,1]
  % omv = mean frequencies of the sinusoids [D,1]
  
  q1 = zlogwspec(1+T*K+K*D:T*K+K*D+D);
  varv = exp(q1);
  
  q2 = zlogwspec(1+T*K+K*D+D:T*K+K*D+2*D);
  lamv = 1./(1+exp(-q2));

  q3 = zlogwspec(1+T*K+K*D+2*D:T*K+K*D+3*D);
  cosomv = (1-exp(-q3))./(1+exp(-q3));
  omv = angle(cosomv + i* sqrt(1-cosomv.^2));

 % varv'
 % lamv'
 % omv'

  
  % form the temporal basis functions
  Z = reshape(zlogwspec(1:T*K),[T,K]);

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

  [ObjA,Xfin,Pfin,VV,VVlag1,VVlag2] = kalman_GTFtNMF_FB_hier_spec(y,A',lamv,varv,omv,vary);
  
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
  
%  dAdX = 1./(1+exp(-X2(1:T,:)*G'-ones(T,1)*Mu2'));
  dY = 1/(4*vary)*(sum(X1X1,3)-(y*ones(1,D)).*X1)./A.^2; %*W).*H;
  dX = (dY*W').*H;
  
  dObjAdZ = zeros(T,K);
  for k=1:K
    dObjAdZ(:,k) = conv(dX(:,k),bas{k},'same');
  end
  
  dObjAdWhat = H'*dY;
  dObj1dW = dObjAdWhat.*W;
  dObj1dW = dObj1dW-diag(sum(dObj1dW,2))*W;

  %%%%%%%%%%%%%%%%%%%%%%%%%
  % Prior component

  ObjB = 1/2*sum(Z(:).^2)*preZ;
  dObjBdZ = Z*preZ;

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-(1-lamv.^2)./(2*varv).*VV(:,1)
%keyboard

%keyboard

sumZZ1 = sum(VV(:,2:T)./(A(2:T,:)'.^2),2);
sumZZ2 = sum(VV(:,1:T-1)./(A(1:T-1,:)'.^2),2);
sumZZlag1 = sum(VVlag1./(A(1:T-1,:)'.*A(2:T,:)'),2);
sumZZlag2 = sum(VVlag2./(A(1:T-1,:)'.*A(2:T,:)'),2);

alp1 = lamv.*cos(omv);
alp2 = lamv.*sin(omv);


dObjAdlogvarv = (T-1) - sumZZ1./(2*varv) ...
                 - lamv.^2.*sumZZ2./(2*varv) ...
                 + alp1.*sumZZlag1./varv ...
                 - alp2.*sumZZlag2./varv ...
                 + 1-(1-lamv.^2)./(2*varv).*VV(:,1)./A(1,:)'.^2; % last bit is
								 % the initial
								 % state contribution
 dObjAdalp1 =   alp1.*sumZZ2./varv - sumZZlag1./varv ...
                -alp1.*VV(:,1)./A(1,:)'.^2./varv + 2*alp1./(1-lamv.^2);

 dObjAdalp2 =   alp2.*sumZZ2./varv + sumZZlag2./varv ...
                -alp2.*VV(:,1)./A(1,:)'.^2./varv + 2*alp2./(1-lamv.^2);

 dObjAdq2 = (cos(omv).*dObjAdalp1+sin(omv).*dObjAdalp2)./(exp(q2/2)+exp(-q2/2)).^2;
 
 dObjAdq3 = lamv.*(dObjAdalp1*2 + dObjAdalp2.*(exp(-q3/2)-exp(q3/2))...
		   )./(exp(q3/2)+exp(-q3/2)).^2;
 
 %%%%%%%%%%%%%%%%%%%%

Obj = (-ObjA+ObjB)/T;
dObj = ([dObjAdZ(:)+dObjBdZ(:);dObj1dW(:);dObjAdlogvarv;dObjAdq2;dObjAdq3])/T;

% Obj = ObjB;
% dObj = [dObjBdZ(:);zeros(D*K,1)];

%Obj = -ObjA;
%dObj = [dObjAdZ(:);dObj1dW(:)];