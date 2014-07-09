function  [lik,Xfin,Pfin] = kalman_GTFtNMF_FB_hier(y,Amp,lamx,varx,omx,vary,varargin);

% function [lik,Xfin,Pfin] = kalman_mPAD_FB_hier(y,Amp,lamx,varx,omx,vary)
% 
% Implements the probabilistic filter bank (see Cemgil and Godsill,
% 2005) with time-varying weights or amplitudes using the hierarchical
% version of the model though, rather than the product version. (A
% variant which is required for GTF-tNMF.) Returns the likelihood and
% sufficient statistics of the Gaussian LDS. Based on Zoubin's
% code. Modified by Richard Turner.
%
% The model can be written as an LDS:
% x_{t}|x_{t-1} ~ Norm(A_t x_{t-1},Q_t)
% y_{t}|x_{t} ~ Norm(C x_{t},R) 
% x_1 ~ Norm(x0,P0)  
%
% where the dynamics are time-varying:
%
% R_dt = a_{1,t}/a_{1,t-1} \lambda_d [cos(\om_d),-sin(\om_d)
%                                     sin(\om_d), cos(\om_d)]
% A_t = [R_{1,t},   0   ,...,   0
%           0    R_{2,t},...,   0
%           0       0   ,...,R_{D,t}]
%
% Q_t = [];
%
% see test_kalman_GTFtNMF_FB_hier.m for unit tests.
% 
% INPUTS:
% y = signal [T,1]
% Amp = matrix of time-varying amplitudes, size [T,D]
% lamx = dynamical AR parameters [D,1]
% varx = dynamical noise parameters [D,1]
% omx = mean frequencies of the sinusoids [D,1]
% vary = observation noise (either scalar or vector [T,1])
%
% OPTIONAL INPUTS:
% verbose = binary scalar, if set to 1 displays progress
%           information
%
% OUTPUTS
% lik = likelihood
% X1 = posterior means of the filter coefficients, size [N,K,T]
% X1X1 = posterior covariance of the filter coefficients, size [K,K,T]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Massage parameters into Kalman form
K = 2*length(omx);
T = length(y);
A = zeros(K);

C = reshape([ones(1,K/2);zeros(1,K/2)],[1,K]);

T_R = length(vary);

x0 = zeros(K,1);

mVar = varx(:)./(1-lamx(:).^2);
temp3 = [mVar'.*Amp(:,1)'.^2;
	 mVar'.*Amp(:,1)'.^2];

P0 = diag(temp3(:));

Y = reshape(y,[1,1,T]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[N D T]=size(Y);
K=length(x0);
tiny=exp(-700);
I=eye(K);
const=(2*pi)^(-D/2);
problem=0;
lik=0;

Xpre=zeros(N,K);   % P(x_t | y_1 ... y_{t-1})
Xcur=zeros(N,K,T);   % P(x_t | y_1 ... y_t)
Xfin=zeros(N,K,T);   % P(x_t | y_1 ... y_T)    given all outputs

Ppre=zeros(K,K,T);
Pcur=zeros(K,K,T);
Pfin=zeros(K,K,T); 

Pt=zeros(K,K); 
Pcov=zeros(K,K); 
Kcur=zeros(K,D);
invP=zeros(D,D);
J=zeros(K,K,T);

if nargin<=6
  verbose = 0 ;
else
  verbose = varargin{1};
end

if nargin<=7
  KF = 0 ;
else
  KF = varargin{2};
end

if verbose==1
  if KF==1
    disp('Kalman Filtering')
  else
    disp('Kalman Smoothing')
  end
end

%%%%%%%%%%%%%%%
% FORWARD PASS

%R=R+(R==0)*tiny;

Xpre=ones(N,1)*x0';
Ppre(:,:,1)=P0;

CntInt=T/5; % T / 2*number of progress values displayed

for t=1:T,
  

  %%%%%%%%
  
  if T_R>1
    R = vary(t);
    invR = inv(R);
  else
    R = vary;
    invR = inv(R);
  end

  if verbose==1&mod(t-1,CntInt)==0
    fprintf(['Progress ',num2str(floor(50*t/T)),'%%','\r'])
  end

  % if (K<D)
  %   temp1= R\C;%  rdiv(C,R);
  %   temp2=temp1*Ppre(:,:,t); % inv(R)*C*Ppre
  %   temp3=C'*temp2;
  %  % temp4=inv(I+temp3)*temp1';
  %   temp4=(I+temp3)\temp1';

  %   invP=invR-temp2*temp4; 
  %   CP= temp1' - temp3*temp4;    % C'*invP
  % else
%    temp1=diag(R)+C*Ppre(:,:,t)*C';
    temp1=R+C*Ppre(:,:,t)*C';
    L = chol(temp1);
%%    invP=inv(temp1);
%%    CP=C'*invP;
%    CP=C'/temp1;
%  end;

%  Kcur=Ppre(:,:,t)*CP;
  Kcur=((Ppre(:,:,t)*C')/L)/L';
  KC=Kcur*C;
  Ydiff=Y(:,:,t)-Xpre*C';
  Xcur(:,:,t)=Xpre+Ydiff*Kcur'; 
  % a more numerically stable version of:   Pcur(:,:,t)=Ppre(:,:,t)-KC*Ppre(:,:,t)
  eye_KC = eye(K)-KC;
  Pcur(:,:,t)=eye_KC*Ppre(:,:,t)*eye_KC'+Kcur*R*Kcur';
%  keyboard

  if (t<T)

    % compute the current setting of the dynamics from the amplitudes
    for l=1:K/2
      ind = [1:2]+(l-1)*2;
      lamEff = lamx(l)*Amp(l,t+1)/Amp(l,t);
      A(ind,ind) = lamEff*[cos(omx(l)),-sin(omx(l));sin(omx(l)),cos(omx(l))];
    end
    
    ASq = Amp(:,t+1).^2;
    temp2 = [varx(:)'.*ASq';
	     varx(:)'.*ASq'];
    
    Q = diag(temp2(:));
        
    Xpre=Xcur(:,:,t)*A';
    Ppre(:,:,t+1)=A*Pcur(:,:,t)*A'+Q;
    
    [LPpre,p] = chol(Ppre(:,:,t+1));
    if p>0
      disp('filtering')
      keyboard  
    end

  end;

  % calculate likelihood
  
  % Old version of the code
  %  detiP=sqrt(det(invP));
  % if (isreal(detiP) & detiP>0)
  %   lik=lik+N*log(detiP)-0.5*sum(sum(Ydiff.*(Ydiff*invP)));
  % else
  %   problem=1;
  % end;

  logdetiP=-sum(log(diag(L)));
  lik=lik+N*logdetiP-0.5*sum(sum(Ydiff.*(Ydiff/L)/L'));

end;  

lik=lik+N*T*log(const);

%%%%%%%%%%%%%%%
% BACKWARD PASS
  
t=T; 
Xfin(:,:,t)=Xcur(:,:,t);
Pfin(:,:,t)=Pcur(:,:,t); 

for t=(T-1):-1:1
  if verbose==1&mod(t+1,CntInt)==0
    fprintf(['Progress ',num2str(50+floor(50*(T-t+1)/T)),'%%','\r'])
  end

  % compute the current setting of the dynamics from the amplitudes
  for l=1:K/2
    ind = [1:2]+(l-1)*2;
    lamEff = Amp(l,t+1)/Amp(l,t)*lamx(l);
    A(ind,ind) = lamEff*[cos(omx(l)),-sin(omx(l));sin(omx(l)),cos(omx(l))];
  end
  
  %%%%%%%%

  [L,p] = chol(Ppre(:,:,t+1));  

  if p>0
    disp('smoothing')
    keyboard
  end
%  J(:,:,t)=Pcur(:,:,t)*A'/Ppre(:,:,t+1);
  J(:,:,t)=((Pcur(:,:,t)*A')/L)/L';
  Xfin(:,:,t)=Xcur(:,:,t)+(Xfin(:,:,t+1)-Xcur(:,:,t)*A')*J(:,:,t)';
  Pfin(:,:,t)=Pcur(:,:,t)+J(:,:,t)*(Pfin(:,:,t+1)-Ppre(:,:,t+1))*J(:,:,t)';
  
end


if verbose==1
  fprintf('                                        \r')
end

