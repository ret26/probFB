function [Ct,Rt,Qt,dCt,dRt,dQt] = getParamsSELocal(covfunc,hyp,vary,X,tau1,tau2,noBlks,missingInd)

sigmaInf = 100;
s = (X(1:tau1));
u = linspace(X(1),X(tau1),tau2)';
tiny = 1e-8;
Kss = feval(covfunc{:},hyp,s)+tiny*eye(tau1);
Kuu = feval(covfunc{:},hyp,u)+tiny*eye(tau2);
Ksu = feval(covfunc{:},hyp,s,u);

L0 = chol(Kuu);
Ksu_invL0 = Ksu/L0;
C = Ksu_invL0/L0';
R = Kss - Ksu_invL0*Ksu_invL0' + tiny*eye(tau1);
Q = Kuu + tiny*eye(tau2);
T = noBlks;
Ct = repmat({C},T,1);
Rt = repmat({R},T,1);
%Qt = repmat({Q},T,1);
Qt = Q;
for t = 1:T
    mInd = missingInd(t,:);
    sigma2 = vary*ones(tau1,1);
    sigma2(mInd) = sigma2(mInd)+sigmaInf; % set large variance for missing indices
    Rt{t} = R + diag(sigma2);
end

if nargout>3
    noVar = length(hyp)+1; 
    dC = cell(noVar,1);
    dR = cell(noVar,1);
    dQ = cell(noVar,1);
    for j = 1:noVar
        if j == noVar
            dCj = zeros(size(C));
            dRj = vary*eye(size(R));
            dQj = zeros(length(Q));
        else
            Kssdi = feval(covfunc{:},hyp,s,[],j);
            Kuudi = feval(covfunc{:},hyp,u,[],j);
            Ksudi = feval(covfunc{:},hyp,s,u,j);
            Ka = Kuudi/Kuu;
            Kb = Ksudi/Kuu;

            dCj = Kb - C*Ka;
            dRj = Kssdi - Kb*Ksu' + C*Ka*Ksu' - C*Ksudi';

            dQj = Kuudi;
        end
        dC{j} = dCj;
        dR{j} = dRj;
        dQ{j} = dQj;
    end
    
    dCt = cell(T,noVar);
    dRt = cell(T,noVar);
    %dQt = cell(T,noVar);
    dQt = dQ;
    for j = 1:noVar
        [dQt{1:T,j}] = deal(dQ{j});
        [dCt{1:T,j}] = deal(dC{j});
        [dRt{1:T,j}] = deal(dR{j});
    end

end


end