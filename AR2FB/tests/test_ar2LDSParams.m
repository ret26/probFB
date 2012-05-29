function test_suite =  test_ar2LDSParams
  initTestSuite;

% tests function [A,Q,C,R,x0,P0] =  ar2LDSParams(Lam,Var,vary);

  
function test_1DExample

Lam = [1.5,-0.7];
Var = 1;
vary = 0.03;

[A,Q,C,R,x0,P0] =  ar2LDSParams(Lam,Var,vary);

A2 = [Lam;
      1,0];
Q2 = [Var,0;
      0, 0];
C2 = [1,0];
R2 = vary;
x02 = [0;0];

autoCor = getAutoCorARTau(Lam',1,2);

P02 = [autoCor';
      autoCor(2:-1:1)'];

tol = 1e-5;
assertVectorsAlmostEqual(A,A2,'absolute',tol,0)
assertVectorsAlmostEqual(Q,Q2,'absolute',tol,0)
assertVectorsAlmostEqual(R,R2,'absolute',tol,0)
assertVectorsAlmostEqual(C,C2,'absolute',tol,0)
assertVectorsAlmostEqual(x0,x02,'absolute',tol,0)
assertVectorsAlmostEqual(P0,P02,'absolute',tol,0)

function test_2DExample

Lam = [1.5,-0.7;
       1,-0.8];
Var = [1,1];
vary = randn(5,1);

[A,Q,C,R,x0,P0] =  ar2LDSParams(Lam,Var,vary);

A2 = [Lam(1,:),0,0;
      1,0,0,0;
      0,0,Lam(2,:);
      0,0,1,0];

Q2 = [Var(1),0,0,0;
      0, 0,0,0;
      0,0,Var(2),0;
      0,0,0,0];
C2 = [1,0,1,0];
R2 = reshape(vary,[1,1,5]);
x02 = [0;0;0;0];

autoCor1 = getAutoCorARTau(Lam(1,:)',1,2);
autoCor2 = getAutoCorARTau(Lam(2,:)',1,2);

P02 = [autoCor1',0,0;
      autoCor1(2:-1:1)',0,0;
      0,0,autoCor2';
      0,0,autoCor2(2:-1:1)'];

tol = 1e-5;
assertVectorsAlmostEqual(A,A2,'absolute',tol,0)
assertVectorsAlmostEqual(Q,Q2,'absolute',tol,0)
assertVectorsAlmostEqual(R,R2,'absolute',tol,0)
assertVectorsAlmostEqual(C,C2,'absolute',tol,0)
assertVectorsAlmostEqual(x0,x02,'absolute',tol,0)
assertVectorsAlmostEqual(P0,P02,'absolute',tol,0)
