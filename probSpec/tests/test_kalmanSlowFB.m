function test_suite = test_kalmanSlowFB
  initTestSuite;

% Tests: 
%
% function % [lik,Xfin,Pfin,Ptsum,YX,A1,A2,A3]=kalmanSlowFB(lamx,varx,om,vary,Y,verbose,KF);
%
% Might want to test complete, over-complete and under-complete
% models separately.

function testTwoStepComputeByHand

% Two time-step example which can be computed by hand

dispFig = 0;

T = 2;
D = 1;


lamx = [.9;.7];
varx = rand(length(lamx),1);
om = [1/10,1/20]';
vary = rand;
mVar = varx./(1-lamx.^2);

K = length(lamx)*2;
x0 = zeros(K,1);

P0 = zeros(K); 
P0(1,1) = mVar(1); P0(2,2) = mVar(1);
P0(3,3) = mVar(2); P0(4,4) = mVar(2);

R1 = lamx(1)*[cos(om(1)),-sin(om(1));sin(om(1)),cos(om(1))];
R2 = lamx(2)*[cos(om(2)),-sin(om(2));sin(om(2)),cos(om(2))];

A = [R1,zeros(2);zeros(2),R2];
C = [1,0,1,0];
Q = diag([varx(1);varx(1);varx(2);varx(2)]);
R = vary;
Y = randn([1,1,T]);

%%%%%%%%%%%%%%%%%%%%%

Y0 = squeeze(Y(1,:,1))';
Y1 = squeeze(Y(1,:,2))';

pre = [inv(P0)+A'*inv(Q)*A+C'*inv(R)*C,-A'*inv(Q);
       -inv(Q)*A,inv(Q)+C'*inv(R)*C];

muPre = [x0'*inv(P0)+Y0'*inv(R)*C,Y1'*inv(R)*C]';

sigmaPost = inv(pre);
muPost = sigmaPost*muPre;

Xfin1 = reshape(muPost,[1,K,T]);

Pfin1 = zeros(K,K,T);
Pfin1(:,:,1) = sigmaPost(1:K,1:K);
Pfin1(:,:,2) = sigmaPost(K+1:2*K,K+1:2*K);

lik1 = -1/2*log(det(2*pi*P0))-1/2*log(det(2*pi*Q)) ...
       -log(det(2*pi*R))+1/2*log(det(2*pi*sigmaPost)) ...
       -1/2*Y0'*inv(R)*Y0-1/2*Y1'*inv(R)*Y1 ...
       -1/2*x0'*inv(P0)*x0 + 1/2*muPost'*inv(sigmaPost)*muPost;

Ptsum_1 = Pfin1(:,:,1)+Pfin1(:,:,2)+...
	  Xfin1(1,:,1)'*Xfin1(1,:,1)+Xfin1(1,:,2)'*Xfin1(1,:,2);

YX_1 = Y(1,:,1)'*Xfin1(1,:,1)+Y(1,:,2)'*Xfin1(1,:,2);

A1_1 = Xfin1(1,:,2)'*Xfin1(1,:,1)+sigmaPost(K+1:2*K,1:K);
A2_1 = Pfin1(:,:,1)+ Xfin1(1,:,1)'*Xfin1(1,:,1);
A3_1 = Pfin1(:,:,2)+ Xfin1(1,:,2)'*Xfin1(1,:,2);

[lik2,Xfin2,Pfin2,Ptsum_2,YX_2,A1_2,A2_2,A3_2] = kalmanSlowFB(lamx, ...
						  varx,om,vary,Y(:),0);

tol = 1e-5;
assertVectorsAlmostEqual(lik1,lik2,'absolute',tol,0)
assertVectorsAlmostEqual(Xfin1,Xfin2,'absolute',tol,0)
assertVectorsAlmostEqual(Pfin1,Pfin2,'absolute',tol,0)
assertVectorsAlmostEqual(Ptsum_1,Ptsum_2,'absolute',tol,0)
assertVectorsAlmostEqual(YX_1,YX_2,'absolute',tol,0)
assertVectorsAlmostEqual(A1_1,A1_2,'absolute',tol,0)
assertVectorsAlmostEqual(A2_1,A2_2,'absolute',tol,0)
assertVectorsAlmostEqual(A3_1,A3_2,'absolute',tol,0)



function testTwoStepComputeByHandFilering

% Two time-step example which can be computed by hand for a
% filtering example

dispFig = 0;

T = 2;
D = 1;


lamx = [.9;.7];
varx = rand(length(lamx),1);
om = [1/10,1/20]';
vary = rand;
mVar = varx./(1-lamx.^2);

K = length(lamx)*2;
x0 = zeros(K,1);

P0 = zeros(K); 
P0(1,1) = mVar(1); P0(2,2) = mVar(1);
P0(3,3) = mVar(2); P0(4,4) = mVar(2);

R1 = lamx(1)*[cos(om(1)),-sin(om(1));sin(om(1)),cos(om(1))];
R2 = lamx(2)*[cos(om(2)),-sin(om(2));sin(om(2)),cos(om(2))];

A = [R1,zeros(2);zeros(2),R2];
C = [1,0,1,0];
Q = diag([varx(1);varx(1);varx(2);varx(2)]);
R = vary;
Y = randn([1,1,T]);

%%%%%%%%%%%%%%%%%%%%%

Y0 = squeeze(Y(1,:,1))';
Y1 = squeeze(Y(1,:,2))';

pre = [inv(P0)+A'*inv(Q)*A+C'*inv(R)*C,-A'*inv(Q);
       -inv(Q)*A,inv(Q)+C'*inv(R)*C];

muPre = [x0'*inv(P0)+Y0'*inv(R)*C,Y1'*inv(R)*C]';

sigmaPostKS = inv(pre);
muPostKS = sigmaPostKS*muPre;

sigmaPostKF = inv(inv(P0)+C'*inv(R)*C);
muPostKF = sigmaPostKF*(x0'*inv(P0)+Y0'*inv(R)*C)';

Xfin1 = zeros(1,K,T);
Xfin1(1,:,1) = reshape(muPostKF,[1,K,1]);
Xfin1(1,:,2) = reshape(muPostKS(K+1:2*K),[1,K,1]);

Pfin1 = zeros(K,K,T);
Pfin1(:,:,1) = sigmaPostKF;
Pfin1(:,:,2) = sigmaPostKS(K+1:2*K,K+1:2*K);

lik1 = -1/2*log(det(2*pi*P0))-1/2*log(det(2*pi*Q)) ...
       -log(det(2*pi*R))+1/2*log(det(2*pi*sigmaPostKS)) ...
       -1/2*Y0'*inv(R)*Y0-1/2*Y1'*inv(R)*Y1 ...
       -1/2*x0'*inv(P0)*x0 + 1/2*muPostKS'*inv(sigmaPostKS)*muPostKS;
% Not yet implement sufficient statistic computation for Kalman
% filtering:

% Ptsum_1 = Pfin1(:,:,1)+Pfin1(:,:,2)+...
% 	  Xfin1(1,:,1)'*Xfin1(1,:,1)+Xfin1(1,:,2)'*Xfin1(1,:,2);

% YX_1 = Y(1,:,1)'*Xfin1(1,:,1)+Y(1,:,2)'*Xfin1(1,:,2);

% A1_1 = Xfin1(1,:,2)'*Xfin1(1,:,1)+sigmaPost(K+1:2*K,1:K);
% A2_1 = Pfin1(:,:,1)+ Xfin1(1,:,1)'*Xfin1(1,:,1);
% A3_1 = Pfin1(:,:,2)+ Xfin1(1,:,2)'*Xfin1(1,:,2);

%[lik2,Xfin2,Pfin2,Ptsum_2,YX_2,A1_2,A2_2,A3_2] = kalman(A,C,Q,R,x0,P0,Y,0);

[lik2,Xfin2,Pfin2] = kalmanSlowFB(lamx,varx,om,vary,Y(:),0,1);

tol = 1e-5;
assertVectorsAlmostEqual(lik1,lik2,'absolute',tol,0)
assertVectorsAlmostEqual(Xfin1(:),Xfin2(:),'absolute',tol,0)
assertVectorsAlmostEqual(Pfin1(:),Pfin2(:),'absolute',tol,0)
% assertVectorsAlmostEqual(Ptsum_1,Ptsum_2,'absolute',tol,0)
% assertVectorsAlmostEqual(YX_1,YX_2,'absolute',tol,0)
% assertVectorsAlmostEqual(A1_1,A1_2,'absolute',tol,0)
% assertVectorsAlmostEqual(A2_1,A2_2,'absolute',tol,0)
% assertVectorsAlmostEqual(A3_1,A3_2,'absolute',tol,0)


% function testMultiStepRuns

% % Test the function runs when there are multiple time-steps

% T = 100;
% D = 3;
% K = 4;

% x0 = randn(K,1);
% P0 = randn(K); P0 = P0*P0';
% A = randn(K);
% C = randn(D,K);
% Q = randn(K); Q = Q*Q';
% R = randn(D); R = R*R';
% Y = randn([1,D,T]);

% [lik2,Xfin2,Pfin2,Ptsum,YX,A1,A2,A3] = kalman(A,C,Q,R,x0,P0,Y,0);

% %tol = 1e-5;
% %assertVectorsAlmostEqual(lik1,lik2,'absolute',tol,0)
% %assertVectorsAlmostEqual(Xfin1,Xfin2,'absolute',tol,0)
% %assertVectorsAlmostEqual(Pfin1,Pfin2,'absolute',tol,0)

