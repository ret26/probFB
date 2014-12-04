function varargout = ...
    kalmanVarLocal(Ct,Q,Rt,x0,P0,Yt,varargin)

T = length(Yt);

problem=0;
lik=0;

Xcur=cell(T,1);   % P(x_t | y_1 ... y_t)

Ppre=cell(T,1);
Pcur=cell(T,1);

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
    Ydiff=Yt{t};
    Xcur{t}=Kcur*Ydiff;
    
    %Pcur{t} = Ppret - KC*Ppret;
    
    % numerical problem with subtraction, use Joseph form
    % for the covariance update
    IsubKC = ItK-KC;
    Pcur{t} = IsubKC*Ppret*IsubKC' + Kcur*Rcur*Kcur';
    
    if (t<T)
        Ppre{t+1}=Q;
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


% Only Kalman filtering
Xfint=Xcur;
Pfint=Pcur;

%% FIND DERIVATIVES
if length(varargin) >= 3
    dCt = varargin{3};
    dQ = varargin{4};
    dRt = varargin{5};
    dP0 = varargin{6};
    noVar = size(dCt,2);
    dlik = zeros(noVar,1);
    Qtinv = Q\eye(size(Q));
    for j = 1:noVar
        M5{j} = -1/2*sum(sum(Qtinv.*dQ{j}'));
        M6{j} = -Qtinv*dQ{j}*Qtinv;
    end
    for i = 1:T
        mu1 = Xfint{i};
        sig11 = Pfint{i};

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
                dlik(j) = dlik(j) + M5{j} + Lb2;
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
                dCmu1 = dCt{i,j}*mu1;
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
end