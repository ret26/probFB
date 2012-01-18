function test_suite =  test_sampleMPAD
  initTestSuite;

% test 
% [Y,X1,X2,A] = sampleMPAD(Params,Dims)

function test_correctAutoCorrelationGPs

dispFigs = 0;
randn('state',1)
vars = [1.5,0.7];
lens = [5.32,3.4]';
tau = 20;
T = 40000;

Lam1 = [1.1,-0.9;
       1.5,-0.95];
Var1 = [1,0.5]';

Dims.D = 2;
Dims.K = 2;
Dims.T = T;

Params.G = randn(2);
Params.Len2 = lens;
Params.Var2 = vars;
Params.Mu2 = randn(2,1);
Params.Var1 = Var1;
Params.Lam1 = Lam1;
Params.vary = 0.1;

[Y,X1,X2,A] = sampleMPAD(Params,Dims);

autoCor1 = vars(1)*exp(-1/(2*lens(1)^2)*([0:tau-1]).^2);
autoCor2 = vars(2)*exp(-1/(2*lens(2)^2)*([0:tau-1]).^2);
 
autoCorEmp1 = zeros(tau,1);
autoCorEmp2 = zeros(tau,1);

for t=1:tau
  autoCorEmp1(t) = mean(X2(1:T-t+1,1).*X2(t:T,1));
  autoCorEmp2(t) = mean(X2(1:T-t+1,2).*X2(t:T,2));
end

if dispFigs==1
  figure;
  subplot(2,1,1)
  hold on
  plot(autoCor1,'-k','linewidth',2)
  plot(autoCorEmp1,'-r','linewidth',1)

  subplot(2,1,2)
  hold on
  plot(autoCor2,'-k','linewidth',2)
  plot(autoCorEmp2,'-r','linewidth',1)

end

tol = 2e-1;
assertVectorsAlmostEqual(autoCor1,autoCorEmp1','absolute',tol,0)
assertVectorsAlmostEqual(autoCor2,autoCorEmp2','absolute',tol,0)


function test_correctAutoCorrelationCars

dispFigs = 1;
randn('state',1)
vars = [1.5,0.7];
lens = [5.32,3.4]';
tau = 20;
T = 100000;

Lam1 = [1.1,-0.9;
       1.5,-0.95];
Var1 = [1,0.5]';

Dims.D = 2;
Dims.K = 2;
Dims.T = T;

Params.G = randn(2);
Params.Len2 = lens;
Params.Var2 = vars;
Params.Mu2 = randn(2,1);
Params.Var1 = Var1;
Params.Lam1 = Lam1;
Params.vary = 0.1;

[Y,X1,X2,A] = sampleMPAD(Params,Dims);

autoCor1A = getAutoCorARTau(Lam1(1,:)',Var1(1),2);
autoCor2A = getAutoCorARTau(Lam1(2,:)',Var1(2),2);

autoCor1B = [var(X1(:,1)); mean(X1(2:T,1).*X1(1:T-1,1))];
autoCor2B = [var(X1(:,2)); mean(X1(2:T,2).*X1(1:T-1,2))];

tol =3e-1;
assertVectorsAlmostEqual(autoCor1A,autoCor1B,'absolute',tol,0)
assertVectorsAlmostEqual(autoCor2A,autoCor2B,'absolute',tol,0)
