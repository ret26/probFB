function [Ct,Rt,At,Qt,dCt,dRt,dAt,dQt] = getParamsSECosMix(covfunc,hyp,mu,vary,X,tau1,tau2,noBlks,missingInd)
% extract dynamical system weights for a uniformly sampled dataset

K = size(hyp,1); % number of components

sigmaInf = 100;
s = (X(1:tau1));
u = linspace(X(1),X(tau1),tau2)';
uprev = u-X(tau1);
tiny = 1e-7;
Ck = cell(K,1);
Rk = cell(K,1);
Ak = cell(K,1);
Qk = cell(K,1);

% find the transition and emission dynamics for each component
Kuuk = cell(K,1);
Ksuk = cell(K,1);
Kupupk = cell(K,1);
Kuupk = cell(K,1);
%R = [];
for k = 1:K
    hypk = hyp(k,:);
    Kss = feval(covfunc{:},hypk,s)+tiny*eye(tau1);
    Kuu = feval(covfunc{:},hypk,u)+tiny*eye(tau2);
    Ksu = feval(covfunc{:},hypk,s,u);
    Kuuk{k} = Kuu;
    Ksuk{k} = Ksu;
    
    L0 = chol(Kuu);
    Ksu_invL0 = Ksu/L0;
    Ck{k} = Ksu_invL0/L0';
    Rk{k} = Kss - Ksu_invL0*Ksu_invL0';
    
    Kupup = feval(covfunc{:}, hypk, uprev)+tiny*eye(tau2);
    Kuup  = feval(covfunc{:}, hypk, u, uprev);
    Kupupk{k} = Kupup;
    Kuupk{k} = Kuup;
    L1 = chol(Kupup);
    Kuup_invL1 = Kuup/L1;
    Ak{k} = Kuup_invL1/L1';
    Qk{k} = Kuu - Kuup_invL1*Kuup_invL1' + tiny*eye(tau2);
end

T = noBlks;
A1 = zeros(2*K*tau2);
Q1 = zeros(2*K*tau2);
for k = 1:K
    A1(2*(k-1)*tau2+(1:tau2),2*(k-1)*tau2+(1:tau2)) = Ak{k};
    A1((2*k-1)*tau2+(1:tau2),(2*k-1)*tau2+(1:tau2)) = Ak{k};
    Q1(2*(k-1)*tau2+(1:tau2),2*(k-1)*tau2+(1:tau2)) = Qk{k};
    Q1((2*k-1)*tau2+(1:tau2),(2*k-1)*tau2+(1:tau2)) = Qk{k};
end
At = A1;
Qt = Q1;

Ct = cell(T,1);
Rt = cell(T,1);

vec1Mat = cell(T,K);
vec2Mat = cell(T,K);
vec3Mat = cell(T,K);
vec4Mat = cell(T,K);
vec12Mat = cell(T,K);

for t = 1:T
    ind = tau1*(t-1) + (1:tau1);
    Ct{t} = [];
    R1 = zeros(tau1);

    for k = 1:K
        mukInd = mu(k)*ind(:);
        cosVec = cos(2*pi*mukInd);
        sinVec = sin(2*pi*mukInd);
        cosMat11 = cosVec(:,ones(1,tau1));
        sinMat11 = sinVec(:,ones(1,tau1));
        cosMat12 = cosVec(:,ones(1,tau2));
        sinMat12 = sinVec(:,ones(1,tau2));
        
        Ct{t} = [Ct{t} cosMat12.*Ck{k} sinMat12.*Ck{k}];
        temp = (cosMat11.*cosMat11' + sinMat11.*sinMat11');
        R1 = R1 + temp.*Rk{k}; 
        vec1Mat{t,k} = cosMat11;
        vec2Mat{t,k} = sinMat11;
        vec3Mat{t,k} = cosMat12;
        vec4Mat{t,k} = sinMat12;
        vec12Mat{t,k} = temp;
    end
    mInd = missingInd(t,:);
    sigma2 = vary*ones(tau1,1);
    sigma2(mInd) = sigma2(mInd)+sigmaInf; % set large variance for missing indices
    Rt{t} = R1 + diag(sigma2);
end

if nargout>4
    % numel(hyp) should be k*2
    D = numel(hyp);
    noVar = D+length(mu)+1; % SE hypers, freqs and obs noise
    dC = cell(noVar,1);
    dR = cell(noVar,1);
    dA = cell(noVar,1);
    dQ = cell(noVar,1);
    for j = 1:noVar
        if j > D
            dCj = zeros(tau1,tau2);
            dRj = zeros(2*K*tau1);
            dAj = zeros(2*K*tau2);
            dQj = zeros(2*K*tau2);
        else
            c = ceil(j/2); % spectral component c
            paramInd = j - (c-1)*2; % param index
            hypc = hyp(c,:);
            Kssdi = feval(covfunc{:},hypc,s,[],paramInd);
            Kuudi = feval(covfunc{:},hypc,u,[],paramInd);
            Ksudi = feval(covfunc{:},hypc,s,u,paramInd);
            Ka = Kuudi/Kuuk{c};
            Kb = Ksudi/Kuuk{c};
            dCj = Kb - Ck{c}*Ka;
            
            dRj = Kssdi - Kb*Ksuk{c}' + Ck{c}*Ka*Ksuk{c}' - Ck{c}*Ksudi';
            
            Kupupdi = feval(covfunc{:},hypc,uprev,[],paramInd);
            Kuupdi = feval(covfunc{:},hypc,u,uprev,paramInd);
            Ka = Kupupdi/Kupupk{c};
            Kb = Kuupdi/Kupupk{c};
            dA1 = Kb - Ak{c}*Ka;
            dQ1 = Kuudi - Kb*Kuupk{c}' + Ak{c}*Ka*Kuupk{c}' - Ak{c}*Kuupdi';
            
            dAj = zeros(2*K*tau2);
            dQj = zeros(2*K*tau2);
            dAj(2*(c-1)*tau2+(1:tau2),2*(c-1)*tau2+(1:tau2)) = dA1;
            dAj((2*c-1)*tau2+(1:tau2),(2*c-1)*tau2+(1:tau2)) = dA1;
            dQj(2*(c-1)*tau2+(1:tau2),2*(c-1)*tau2+(1:tau2)) = dQ1;
            dQj((2*c-1)*tau2+(1:tau2),(2*c-1)*tau2+(1:tau2)) = dQ1;
        end
        dC{j} = dCj; % tau1 x tau2
        dR{j} = dRj; % 2*K*tau1 x 2*K*tau1
        dA{j} = dAj; % 2*K*tau2 x 2*K*tau2
        dQ{j} = dQj; % 2*K*tau2 x 2*K*tau2
    end
    
    dCt = cell(T,noVar);
    dRt = cell(T,noVar);
    dAt = dA;
    dQt = dQ;
    
    for t = 1:T
        for j = 1:noVar
            if j == noVar % observation noise
                dCt{t,j} = zeros(tau1,2*K*tau2);
                sigma2 = vary*ones(tau1,1);
                dRt{t,j} = diag(sigma2);
            elseif j>D && j<=noVar-1 % frequencies
                c = j-D; % spectral component c
                ind = tau1*(t-1) + (1:tau1);
                pi2mut = 2*pi*mu(c)*ind(:);
                
                dvec1 = -pi2mut.*sin(pi2mut);
                dvec2 = pi2mut.*cos(pi2mut);
                dvec1Mat = dvec1(:,ones(1,tau1));
                dvec2Mat = dvec2(:,ones(1,tau1));
                
                dvec3Mat = dvec1(:,ones(1,tau2));
                dvec4Mat = dvec2(:,ones(1,tau2));
                
                
                dCt{t,j} = [dvec3Mat.*Ck{c} dvec4Mat.*Ck{c}];
                dRt{t,j} = (dvec1Mat.*vec1Mat{t,c}' + dvec2Mat.*vec2Mat{t,c}' + ...
                    vec1Mat{t,c}.*dvec1Mat' + vec2Mat{t,c}.*dvec2Mat').*Rk{c};
            else % other params
                c = ceil(j/2); % spectral component c
                
                dCt{t,j} = [vec3Mat{t,c}.*dC{j} vec4Mat{t,c}.*dC{j}];
                %dRt{t,j} = dR{j}.*(vec1Mat{t,c}.*vec1Mat{t,c}' + vec2Mat{t,c}.*vec2Mat{t,c}');
                dRt{t,j} = dR{j}.*vec12Mat{t,c};
            end
        end
        
    end
    
end
end