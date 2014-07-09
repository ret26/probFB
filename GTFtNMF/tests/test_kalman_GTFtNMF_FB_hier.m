function test_suite = test_kalman_GTFtNMF_FB_hier
  initTestSuite;

% Tests: 
%
% [lik,Xfin,Pfin] = kalman_GTFtNMF_FB_hier(y,Amp,lamx,varx,omx,vary);
%



function testTwoStepComputeByHand

% Two time-step example which can be computed by hand

dispFig = 0;

T = 2;
D = 2;


lamx = [.9;.7];
varx = rand(length(lamx),1);
om = [1/10,1/20]';
vary = rand;
mVar = varx./(1-lamx.^2);


Amp = rand(D,T);

K = length(lamx)*2;
x0 = zeros(K,1);

P0 = zeros(K); 
P0(1,1) = Amp(1,1).^2.*mVar(1); P0(2,2) = Amp(1,1).^2.*mVar(1);
P0(3,3) = Amp(2,1).^2.*mVar(2); P0(4,4) = Amp(2,1).^2.*mVar(2);

R1 = Amp(1,2)/Amp(1,1)*lamx(1)*[cos(om(1)),-sin(om(1));sin(om(1)),cos(om(1))];
R2 = Amp(2,2)/Amp(2,1)*lamx(2)*[cos(om(2)),-sin(om(2));sin(om(2)),cos(om(2))];

A = [R1,zeros(2);zeros(2),R2];
C = [1,0,1,0];
Q = diag([Amp(1,2).^2.*varx(1);Amp(1,2).^2.*varx(1); ...
	  Amp(2,2).^2.*varx(2);Amp(2,2).^2.*varx(2)]);
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


[lik2,Xfin2,Pfin2] = kalman_GTFtNMF_FB_hier(Y(:),Amp,lamx,varx,om,vary);
%[lik2,Xfin2,Pfin2] = kalman_mPAD_FB_hier(Params,Y(:),Amp);
%[lik2,Xfin2,Pfin2] = kalmanSlowFB(lamx,varx,om,vary,Y(:));


tol = 1e-5;
assertVectorsAlmostEqual(lik1,lik2,'absolute',tol,0)
assertVectorsAlmostEqual(Xfin1,Xfin2,'absolute',tol,0)
assertVectorsAlmostEqual(Pfin1,Pfin2,'absolute',tol,0)



