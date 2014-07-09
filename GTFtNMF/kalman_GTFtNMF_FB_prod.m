function  [lik,Xfin,Pfin] = kalman_GTFtNMF_FB_prod(y,Amp,lamx,varx,omx,vary,varargin);

% function [lik,Xfin,Pfin] = kalman_mPAD_FB_prod(y,Amp,lamx,varx,omx,vary)
% 
% Implements the probabilistic filter bank (see Cemgil and Godsill,
% 2005) with time-varying weights or amplitudes using the product
% version of the model though, rather than the hierarchical version. (A
% variant which is required for GTF-tNMF.) Returns the likelihood and
% sufficient statistics of the Gaussian LDS. Based on Zoubin's
% code. Modified by Richard Turner.
%
% The model can be written as an LDS:
% x_{t}|x_{t-1} ~ Norm(A x_{t-1},Q_t)
% y_{t}|x_{t} ~ Norm(C_t x_{t},R) 
% x_1 ~ Norm(x0,P0)  
%
% where the dynamics are fixed:
%
% A = \lambda_d [cos(\om_d),-sin(\om_d)
%                   sin(\om_d), cos(\om_d)]
% C_t = [Amp_{1,t},   0   ,...,   0
%           0    Amp_{2,t},...,   0
%           0       0   ,...,Amp_{D,t}]
%
% Q_t = [];
%
% see test_kalman_GTFtNMF_FB_prod.m for unit tests.
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

T_R = length(vary);

x0 = zeros(K,1);

mVar = varx(:)./(1-lamx(:).^2);
temp3 = [mVar';
	 mVar'];

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

%Ppre=zeros(K,K,T);
Lpre=zeros(K,K,T);
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
%Ppre(:,:,1)=P0;
Lpre(:,:,1)=chol(P0);

CntInt=T/5; % T / 2*number of progress values displayed

% compute the setting of the dynamics from the amplitudes
for l=1:K/2
  ind = [1:2]+(l-1)*2;
  A(ind,ind) = lamx(l)*[cos(omx(l)),-sin(omx(l));sin(omx(l)),cos(omx(l))];
end

temp2 = [varx(:)';
	 varx(:)'];

Q = diag(temp2(:));
    
for t=1:T,
  
  %%%%%%%%%
  % compute the current setting of the weights
  C = reshape([Amp(:,t)';zeros(1,K/2)],[1,K]);

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

%    temp1=R+C*Ppre(:,:,t)*C';
    CLpre = C*Lpre(:,:,t)';
    temp1=R+CLpre*CLpre';  %C*Lpre(:,:,t)'*C';
    [L,p1] = chol(temp1);

%    if p1>0
%      disp('error filtering 1')
%      keyboard
%    end
%%    invP=inv(temp1);
%%    CP=C'*invP;
%    CP=C'/temp1;
%  end;

%  Kcur=Ppre(:,:,t)*CP;
  Kcur=(Lpre(:,:,t)'*CLpre')/temp1;

  KC=Kcur*C;
  Ydiff=Y(:,:,t)-Xpre*C';
  Xcur(:,:,t)=Xpre+Ydiff*Kcur'; 
  % a more numerically stable version of:   Pcur(:,:,t)=Ppre(:,:,t)-KC*Ppre(:,:,t)
%  eye_KC = (eye(K)-KC);
%  Pcur(:,:,t)=eye_KC*Ppre(:,:,t)*eye_KC'+Kcur*R*Kcur';

  eye_KC_Lpre = (eye(K)-KC)*Lpre(:,:,t)';
  Pcur(:,:,t)=eye_KC_Lpre*eye_KC_Lpre'+Kcur*R*Kcur';
  [Lcur,p2] = chol(Pcur(:,:,t));

%    if p2>0
%      disp('error filtering 2')
%      keyboard
%    end

%  keyboard

  if (t<T)
        
    Xpre=Xcur(:,:,t)*A';
%    Ppre(:,:,t+1)=A*Pcur(:,:,t)*A'+Q;
%    Lpre(:,:,t+1) = chol(A*Pcur(:,:,t)*A'+Q);

    ALcur = A*Lcur';
    [Lpre(:,:,t+1),p3] = chol(ALcur*ALcur'+Q);

    J(:,:,t)=((Pcur(:,:,t)*A')/Lpre(:,:,t+1))/Lpre(:,:,t+1)';

 %   if p3>0
 %     disp('error filtering 3')
 %     keyboard
 %   end

 %   cc = cond(Lpre(:,:,t+1));

 %   if cc>1e15
 %     disp('condition number Lpre filter')
 %     keyboard
 %   end

 %   [LPpre,p] = chol(Ppre(:,:,t+1));

  end;

  % calculate likelihood
  
  % Old version of the code
  %  detiP=sqrt(det(invP));
  % if (isreal(detiP) & detiP>0)
  %   lik=lik+N*log(detiP)-0.5*sum(sum(Ydiff.*(Ydiff*invP)));
  % else
  %   problem=1;
  % end;

%  logdetiP=-sum(log(diag(L)));
%  lik=lik+N*logdetiP-0.5*sum(sum(Ydiff.*(Ydiff/L)/L'));

  logdetiP=-1/2*sum(log(temp1));
  lik=lik+N*logdetiP-0.5*Ydiff^2/temp1;

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

  
  %%%%%%%%

%  Ppre =  Ppre(:,:,t+1);
%  Ppre = A*Pcur(:,:,t)*A'+Q;
%  [L,p] = chol(Ppre);  

  L = Lpre(:,:,t+1);
%  cc = cond(L);

%    if cc>1e15
%      disp('condition number Lpre smoother')
%      keyboard
%    end

%  if p>0
%    disp('smoothing')
%    keyboard
%  end
%  J(:,:,t)=Pcur(:,:,t)*A'/Ppre(:,:,t+1);
%  J(:,:,t)=((Pcur(:,:,t)*A')/L)/L';
  Xfin(:,:,t)=Xcur(:,:,t)+(Xfin(:,:,t+1)-Xcur(:,:,t)*A')*J(:,:,t)';
%  Pfin(:,:,t)=Pcur(:,:,t)+J(:,:,t)*(Pfin(:,:,t+1)-Ppre)*J(:,:,t)';
  Pfin(:,:,t)=Pcur(:,:,t)+J(:,:,t)*(Pfin(:,:,t+1)-L'*L)*J(:,:,t)';  

end


if verbose==1
  fprintf('                                        \r')
end

