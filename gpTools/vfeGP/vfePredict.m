function [mf,vf] = vfePredict(covfunc,hyp,x,y,xtest)

tiny = 1e-7;

Kuu = feval(covfunc{:},hyp.cov,hyp.Xu);
Kuu = Kuu + tiny*eye(size(Kuu));
Kfu = feval(covfunc{:},hyp.cov,x,hyp.Xu);
Kuf = Kfu';
Ktu = feval(covfunc{:},hyp.cov,xtest,hyp.Xu);
Ktt = feval(covfunc{:},hyp.cov,xtest,'diag');

sigma2 = exp(2*hyp.lik);
A = Kuu + Kuf*Kfu/sigma2;
LA = chol(A);

Kufy = Kuf*y;
invAKufy = solve_chol(LA,Kufy);
mf = Ktu*invAKufy/sigma2;

Luu = chol(Kuu);
KtuInvLuu = Ktu/Luu;
KtuInvLa = Ktu/LA;
vf = Ktt-sum(KtuInvLuu.*KtuInvLuu,2)+sum(KtuInvLa.*KtuInvLa,2);

end
