function test_suite = test_AR2FBFFT
  initTestSuite;

% Tests: 
%
% function [X,covX] = AR2FBFFT(y,Lam,Var,vary)
%

function testMultiStep

% Test the function returns the true variables when the observation
% noise is zero and the data are sampled from a single AR(2) process
dispFigs=0;

 T = 100;
 D = 1;
 K = 1;

 Lam = [1.5,-0.9];
 Var = [1];

 vary = 1e-6;
 
x1 = sampleARTau(Lam(1,:),Var(1),T,1);
Y = x1';

[X,covX] = AR2FB(Y,Lam,Var,vary);
X2 = AR2FBFFT(Y,Lam,Var,vary);

if dispFigs==1
figure
hold on
plot(X,'-k','linewidth',2)
plot(X2,'-r')
end

tol = 1e-3;
assertVectorsAlmostEqual(x1(:),X2(:),'absolute',tol,0)
assertVectorsAlmostEqual(X2(:),X(:),'absolute',tol,0)


function testMultiStepMultiChain

% Test the function agrees with the Kalman Smoother 

randn('state',1);

dispFigs=1;

 T = 200;
 D = 1;
 K = 2;

 Lam = [1.5,-0.9;
        1.9,-0.95];
 Var = [1,0.5];

 vary = 1e-1;
 
x1 = sampleARTau(Lam(1,:),Var(1),T,1);
x2 = sampleARTau(Lam(2,:),Var(2),T,1);
Y = (x1+x2)';

X = AR2FB(Y,Lam,Var,vary);
X2 = AR2FBFFT(Y,Lam,Var,vary);

if dispFigs==1
  figure
  hold on
  plot(X','-k','linewidth',2)
  plot(X2','-r')
  error=mean((X(:)-X2(:)).^2)
end



tol = 1e-2;
tEdge = 40;
assertVectorsAlmostEqual(X2(:,tEdge:T-tEdge),X(:,tEdge:T-tEdge),'absolute',tol,0)

