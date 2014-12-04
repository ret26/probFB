function varargout = ...
    kalmanVar(A,Ct,Q,Rt,x0,P0,Yt,varargin)

T = length(Yt);

problem=0;
lik=0;

Xcur=cell(T,1);   % P(x_t | y_1 ... y_t)
Xfint=cell(T,1);   % P(x_t | y_1 ... y_T)    given all outputs

Ppre=cell(T,1);
Pcur=cell(T,1);
Pfint=cell(T,1);
J = cell(T,1);

if nargin<=7
    verbose = 0 ;
else
    verbose = varargin{1};
end
if nargin<=8
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

Xpre=x0;

P0tin = P0 + 1e-7*eye(length(P0));
Ppre{1} = P0tin;

CntInt=T/5; % T / 2*number of progress values displayed

for t=1:T
    Rcur = Rt{t};
    Ccur = Ct{t};
    Kt = size(Ccur,2);
    Dt = size(Rcur,1);
    ItK = eye(Kt);
    ItD = eye(Dt);
    Ppret = Ppre{t};
    
    if verbose==1&&mod(t-1,CntInt)==0
        fprintf(['Progress ',num2str(floor(50*t/T)),'%%','\r'])
    end
    
    temp1=Rcur+Ccur*Ppret*Ccur';
    invP=ItD/temp1;
    CP=Ccur'*invP;
    Kcur=Ppret*CP; % Kalman gain
    KC=Kcur*Ccur;
    Ydiff=Yt{t}-Ccur*Xpre;
    Xcur{t}=Xpre+Kcur*Ydiff;
    
    %Pcur{t} = Ppret - KC*Ppret;
    
    % numerical problem with subtraction, use Joseph form
    % for the covariance update
    IsubKC = ItK-KC;
    Pcur{t} = IsubKC*Ppret*IsubKC' + Kcur*Rcur*Kcur';
    
    if (t<T)
        Xpre=A*Xcur{t};
        Ppre{t+1}=A*Pcur{t}*A'+Q;
    end
    
    % calculate likelihood
    if length(varargin) >= 3
        %P = Rcur + Ccur*Ppret*Ccur';
        %[cholP,e2] = chol(P);
        [cholP,e2] = chol(temp1);
        if e2
            problem = 1;
            continue;
        else
            logdetP=sum(log(diag(cholP)));
            lik=lik-Dt/2*log(2*pi)-logdetP-0.5*sum(sum(Ydiff.*(invP*Ydiff)));
        end
        
    end
end

% Kalman Filtering

if KF==1
    % Only Kalman filtering
    Xfint=Xcur;
    Pfint=Pcur;
else
    % Kalman Smoothing
    
    %%%%%%%%%%%%%%%
    % BACKWARD PASS
    
    t=T;
    Xfint{t}=Xcur{t};
    Pfint{t}=Pcur{t};
    
    for t=(T-1):-1:1
        if verbose==1&&mod(t+1,CntInt)==0
            fprintf(['Progress ',num2str(50+floor(50*(T-t+1)/T)),'%%','\r'])
        end
        J{t}=Pcur{t}*A'/Ppre{t+1};
        Xfint{t}=Xcur{t}+J{t}*(Xfint{t+1}-A*Xcur{t});
        Pfint{t}=Pcur{t}+J{t}*(Pfint{t+1}-Ppre{t+1})*J{t}';
    end
end


A4=cell(T-1,1);
Pcov = Pfint{T}*J{T-1}';
A4{T-1}=Pcov;
for t=(T-1):-1:2
    Pcov = Pfint{t}*J{t-1}';
    A4{t-1} = Pcov;
    
end

%% FIND DERIVATIVES
if length(varargin) >= 3
    dA = varargin{3};
    dCt = varargin{4};
    dQ = varargin{5};
    dRt = varargin{6};
    dP0 = varargin{7};
    noVar = size(dCt,2);
    dlik = zeros(noVar,1);
    Qtinv = Q\eye(size(Q));
    
    for j = 1:noVar
        M5{j} = -1/2*sum(sum(Qtinv.*dQ{j}'));
        M6{j} = -Qtinv*dQ{j}*Qtinv;
        M7{j} = Qtinv*dA{j} + M6{j}*A;
        M8{j} = dA{j}'*(Qtinv*A) + A'*M7{j};
    end
    for i = 1:T
        mu1 = Xfint{i};
        sig11 = Pfint{i};
        
        if i>1
            mu2 = Xfint{i-1};
            sig22 = Pfint{i-1};
            sig12 = A4{i-1};
        end
        
        Rinv = Rt{i}\eye(size(Rt{i}));
        RinvY = Rinv*Yt{i};
        RinvC = Rinv*Ct{i};
        RinvCmu1 = RinvC*mu1;
        RinvCsig11 = RinvC*sig11;
        RinvCsig11CtRinv = RinvCsig11*RinvC';
        for j = 1:noVar
            if i==1
                P0inv = P0tin\eye(size(P0tin));
                M15 = -1/2*sum(sum(P0inv.*dP0{j}'));
                M16 = -P0inv*dP0{j}*P0inv;
                dlik(j) = dlik(j) + M15 - 1/2*(sum(sum(M16.*sig11')) + mu1'*M16*mu1);
            else
                Lb2 = -1/2*(sum(sum(M6{j}.*sig11')) + mu1'*M6{j}*mu1);
                Lb4 = -1/2*(sum(sum(M8{j}.*sig22')) + mu2'*M8{j}*mu2);
                Lb31 = mu1'*M7{j}*mu2;
                Lb32 = sum(sum(sig12'.*M7{j}')); % only work when sig12 is square
                Lb3 = Lb31 + Lb32;
                dlik(j) = dlik(j) + M5{j} + Lb2 + Lb3 + Lb4;
            end
            
            if j==noVar % dRt diag and dCt = 0
                Z1 = -1/2*sum(diag(Rinv).*diag(dRt{i,j}'));
                Z2 = 1/2*RinvY'*dRt{i,j}*RinvY;
                Z3 = -RinvY'*dRt{i,j}*RinvCmu1;
                Z4 = 1/2*sum(diag(dRt{i,j}).*diag(RinvCsig11CtRinv));
                Z5 = 1/2*(RinvCmu1'*dRt{i,j}*RinvCmu1);
                dlik(j) = dlik(j) + Z1 + Z2 + Z3 + Z4 + Z5;
            else
                Z1 = -1/2*sum(sum(Rinv.*dRt{i,j}'));
                Z2 = 1/2*RinvY'*dRt{i,j}*RinvY;
                dCmu1 = dCt{i,j}*mu1(:);
                Z3 = -RinvY'*dRt{i,j}*RinvCmu1 + RinvY'*dCmu1;
                
                Bk = RinvCsig11;
                dCktrBk = dCt{i,j}'*Bk;
                Z41 = -sum(diag(dCktrBk));
                Z42 =  1/2*sum(sum(dRt{i,j}.*RinvCsig11CtRinv));
                Z4 = Z41+Z42;
                Z5 = -1/2*(2*dCmu1'*RinvCmu1 - RinvCmu1'*dRt{i,j}*RinvCmu1);
                dlik(j) = dlik(j) + Z1 + Z2 + Z3 + Z4 + Z5;
            end
        end
    end
end

if problem
    fprintf('!!!!!!!!! problem  !!!!!!!!!\r\n');
    %lik = NaN;
end

if verbose==1
    fprintf('                                        \r')
end

if length(varargin) >= 3
    varargout{1} = lik;
    varargout{2} = dlik;
else
    varargout{1} = Xfint;
    varargout{2} = Pfint;
    varargout{3} = A4;
end