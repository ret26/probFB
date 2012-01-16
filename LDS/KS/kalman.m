%
% function  [lik,Xfin,Pfin,Ptsum,YX,A1,A2,A3]=kalmansmooth(A,C,Q,R,x0,P0,Y);
% 
% Kalman Smoother
%
% x_{t}|x_{t-1} ~ Norm(A x_{t-1},Q)
% y_{t}|x_{t} ~ Norm(C x_{t},R) 
% x_1 ~ Norm(x0,P0)  
% 
% INPUTS  
% A = Dynamical Matrix [K,K]
% C = Emission Matrices  [p,K]
% Q = State innovations noise [K,K]
% R = Emission Noise = [p,p]  
% x0 = initial state mean [K,1]
% P0 = initial state covariance [K,K]  
% Y = Data [N,p,T] 
%
% OUTPUTS
% lik = likelihood
% Xfin = posterior means [N,K,T]
% Pfin = posterior covariance [K,K,T]
% Ptsum = \sum_{t=1}^T <x_{k t} x_{k' t}> [K,K]
% YX = \sum_{t=1}^T y_{t}<x_{k t}> [p,K]
% A1 = \sum_{t=2}^T <x_{k t} x_{k' t-1}> [K,K] 
% A2 = \sum_{t=2}^T <x_{k t-1} x_{k' t-1}> [K,K]
% A3 = \sum_{t=2}^T <x_{k t} x_{k' t}> [K,K]

function  [lik,Xfin,Pfin,Ptsum,YX,A1,A2,A3]=kalmansmooth(A,C,Q,R,x0,P0,Y);

[N p T]=size(Y);
K=length(x0);
tiny=exp(-700);
I=eye(K);
const=(2*pi)^(-p/2);
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
Kcur=zeros(K,p);
invP=zeros(p,p);
J=zeros(K,K,T);

%%%%%%%%%%%%%%%
% FORWARD PASS

R=R+(R==0)*tiny;

Xpre=ones(N,1)*x0';
Ppre(:,:,1)=P0;
invR=diag(1./R);

CntInt=T/5; % T / 2*number of progress values displayed

for t=1:T,
  if mod(t-1,CntInt)==0
    disp(['Progress ',num2str(floor(50*t/T)),'%'])
  end
    
    
  if (K<p)
    temp1=rdiv(C,R);
    temp2=temp1*Ppre(:,:,t); % inv(R)*C*Ppre
    temp3=C'*temp2;
    temp4=inv(I+temp3)*temp1';
    invP=invR-temp2*temp4; 
    CP= temp1' - temp3*temp4;    % C'*invP
  else
    temp1=diag(R)+C*Ppre(:,:,t)*C';
    invP=inv(temp1);
    CP=C'*invP;
  end;
  Kcur=Ppre(:,:,t)*CP;
  KC=Kcur*C;
  Ydiff=Y(:,:,t)-Xpre*C';
  Xcur(:,:,t)=Xpre+Ydiff*Kcur'; 
  Pcur(:,:,t)=Ppre(:,:,t)-KC*Ppre(:,:,t);
  if (t<T)
    Xpre=Xcur(:,:,t)*A';
    Ppre(:,:,t+1)=A*Pcur(:,:,t)*A'+Q;
  end;

  % calculate likelihood
  detiP=sqrt(det(invP));
  if (isreal(detiP) & detiP>0)
    lik=lik+N*log(detiP)-0.5*sum(sum(Ydiff.*(Ydiff*invP)));
  else
    problem=1;
  end;
end;  

lik=lik+N*T*log(const);

%%%%%%%%%%%%%%%
% BACKWARD PASS

A1=zeros(K);
A2=zeros(K);
A3=zeros(K);
Ptsum=zeros(K);
YX=zeros(p,K);
  
t=T; 
Xfin(:,:,t)=Xcur(:,:,t);
Pfin(:,:,t)=Pcur(:,:,t); 
Pt=Pfin(:,:,t) + Xfin(:,:,t)'*Xfin(:,:,t)/N; 
A2= -Pt;
Ptsum=Pt;

YX=Y(:,:,t)'*Xfin(:,:,t);

for t=(T-1):-1:1
  if mod(t+1,CntInt)==0
    disp(['Progress ',num2str(50+floor(50*(T-t+1)/T)),'%'])
  end

  J(:,:,t)=Pcur(:,:,t)*A'*inv(Ppre(:,:,t+1));
  Xfin(:,:,t)=Xcur(:,:,t)+(Xfin(:,:,t+1)-Xcur(:,:,t)*A')*J(:,:,t)';
  Pfin(:,:,t)=Pcur(:,:,t)+J(:,:,t)*(Pfin(:,:,t+1)-Ppre(:,:,t+1))*J(:,:,t)';
  Pt=Pfin(:,:,t) + Xfin(:,:,t)'*Xfin(:,:,t)/N; 
  Ptsum=Ptsum+Pt;
  YX=YX+Y(:,:,t)'*Xfin(:,:,t);
end;

A3 = Ptsum-Pt;
A2 = Ptsum+A2;

t=T;  
Pcov=(I-KC)*A*Pcur(:,:,t-1);
A1=A1+Pcov+Xfin(:,:,t)'*Xfin(:,:,t-1)/N;

for t=(T-1):-1:2
  Pcov=(Pcur(:,:,t)+J(:,:,t)*(Pcov-A*Pcur(:,:,t)))*J(:,:,t-1)';
  A1=A1+Pcov+Xfin(:,:,t)'*Xfin(:,:,t-1)/N;
end;    

if problem 
  fprintf('!!!!!!!!! problem  !!!!!!!!!'); problem=0;
  lik = NaN;
end;

% if ( ~all(eig(Ppre(:,:,t+1))>0) | ~all(eig(Pcur(:,:,t))>0)),
% disp('negative variances !');
% keyboard;
% end;

% Hacked by Rich so that the output is the kalman filter
%Xfin=Xcur;
%Pfin=Pcur; 
