function test_suite = test_kalman_GTFtNMF_FB_hier_spec
  initTestSuite;

% Tests: 
%
% function  [lik,Xfin,Pfin,VV,VVlag] = kalman_GTFtNMF_FB_hier_spec(y,Amp,lamx,varx,omx,vary,varargin);
%


function testThreeStepComputeByHand

% Two time-step example which can be computed by hand

dispFig = 0;

T = 3;
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

R11 = Amp(1,2)/Amp(1,1)*lamx(1)*[cos(om(1)),-sin(om(1));sin(om(1)),cos(om(1))];
R21 = Amp(2,2)/Amp(2,1)*lamx(2)*[cos(om(2)),-sin(om(2));sin(om(2)),cos(om(2))];

R12 = Amp(1,3)/Amp(1,2)*lamx(1)*[cos(om(1)),-sin(om(1));sin(om(1)),cos(om(1))];
R22 = Amp(2,3)/Amp(2,2)*lamx(2)*[cos(om(2)),-sin(om(2));sin(om(2)),cos(om(2))];

A1 = [R11,zeros(2);zeros(2),R21];
A2 = [R12,zeros(2);zeros(2),R22];
C = [1,0,1,0];
Q1 = diag([Amp(1,2).^2.*varx(1);Amp(1,2).^2.*varx(1); ...
	  Amp(2,2).^2.*varx(2);Amp(2,2).^2.*varx(2)]);
Q2 = diag([Amp(1,3).^2.*varx(1);Amp(1,3).^2.*varx(1); ...
	  Amp(2,3).^2.*varx(2);Amp(2,3).^2.*varx(2)]);

R = vary;
Y = randn([1,1,T]);

%%%%%%%%%%%%%%%%%%%%%

Y0 = squeeze(Y(1,:,1))';
Y1 = squeeze(Y(1,:,2))';
Y2 = squeeze(Y(1,:,3))';

pre = [inv(P0)+A1'*inv(Q1)*A1+C'*inv(R)*C,-A1'*inv(Q1),zeros(2*D);
       -inv(Q1)*A1,inv(Q1)+A2'*inv(Q2)*A2+C'*inv(R)*C,-A2'*inv(Q2);
        zeros(2*D), -inv(Q2)*A2, inv(Q2)+C'*inv(R)*C];

muPre = [x0'*inv(P0)+Y0'*inv(R)*C,Y1'*inv(R)*C,Y2'*inv(R)*C]';

sigmaPost = inv(pre);
muPost = sigmaPost*muPre;

Xfin1 = reshape(muPost,[1,K,T]);

Pfin1 = zeros(K,K,T);
Pfin1(:,:,1) = sigmaPost(1:K,1:K);
Pfin1(:,:,2) = sigmaPost(K+1:2*K,K+1:2*K);
Pfin1(:,:,3) = sigmaPost(2*K+1:3*K,2*K+1:3*K);

lik1 = -1/2*log(det(2*pi*P0))-1/2*log(det(2*pi*Q1)) -1/2*log(det(2*pi*Q2))...
       -3/2*log(det(2*pi*R))+1/2*log(det(2*pi*sigmaPost)) ...
       -1/2*Y0'*inv(R)*Y0-1/2*Y1'*inv(R)*Y1 ...
       -1/2*Y2'*inv(R)*Y2 ...
       -1/2*x0'*inv(P0)*x0 + 1/2*muPost'*inv(sigmaPost)*muPost;

VV1 = zeros(D,T);
VV1(1,1) = Xfin1(1,1,1).^2+Xfin1(1,2,1).^2+Pfin1(1,1,1)+Pfin1(2,2,1);
VV1(1,2) = Xfin1(1,1,2).^2+Xfin1(1,2,2).^2+Pfin1(1,1,2)+Pfin1(2,2,2);
VV1(1,3) = Xfin1(1,1,3).^2+Xfin1(1,2,3).^2+Pfin1(1,1,3)+Pfin1(2,2,3);
VV1(2,1) = Xfin1(1,3,1).^2+Xfin1(1,4,1).^2+Pfin1(3,3,1)+Pfin1(4,4,1);
VV1(2,2) = Xfin1(1,3,2).^2+Xfin1(1,4,2).^2+Pfin1(3,3,2)+Pfin1(4,4,2);
VV1(2,3) = Xfin1(1,3,3).^2+Xfin1(1,4,3).^2+Pfin1(3,3,3)+Pfin1(4,4,3);

RVV1 = zeros(D,T-1);

R1 = [1,1;-1,1];
R2 = R1;


VVlag11(1,1) = (sigmaPost(1,5)+muPost(1)*muPost(5)) ...
            + (sigmaPost(2,6)+muPost(2)*muPost(6));

VVlag11(2,1) = (sigmaPost(3,7)+muPost(3)*muPost(7)) ...
            + (sigmaPost(4,8)+muPost(4)*muPost(8));

VVlag11(1,2) = (sigmaPost(5,9)+muPost(5)*muPost(9)) ...
            + (sigmaPost(6,10)+muPost(6)*muPost(10));

VVlag11(2,2) = (sigmaPost(7,11)+muPost(7)*muPost(11)) ...
            + (sigmaPost(8,12)+muPost(8)*muPost(12));


VVlag21(1,1) = (sigmaPost(2,5)+muPost(2)*muPost(5)) ...
            - (sigmaPost(1,6)+muPost(1)*muPost(6)); 

VVlag21(2,1) = (sigmaPost(4,7)+muPost(4)*muPost(7)) ...
            - (sigmaPost(3,8)+muPost(3)*muPost(8));
            

VVlag21(1,2) = (sigmaPost(6,9)+muPost(6)*muPost(9)) ...
            - (sigmaPost(5,10)+muPost(5)*muPost(10));

VVlag21(2,2) = (sigmaPost(8,11)+muPost(8)*muPost(11)) ...
            - (sigmaPost(7,12)+muPost(7)*muPost(12));


[lik2,Xfin,Pfin,VV2,VVlag12,VVlag22] = kalman_GTFtNMF_FB_hier_spec(Y(:),Amp,lamx,varx,om,vary);

% % alternative method:
% RVV3 = zeros(D,T-1);
% EXX = sigmaPost+muPost*muPost';

% temp = EXX(1:2,5:6)';
% RVV3(1,1) = R1(:)'*temp(:);

% temp = EXX(3:4,7:8)';
% RVV3(2,1) = R2(:)'*temp(:);

% temp = EXX(5:6,9:10)';
% RVV3(1,2) = R1(:)'*temp(:);

% temp = EXX(7:8,11:12)';
% RVV3(2,2) = R2(:)'*temp(:);

tol = 1e-5;
assertVectorsAlmostEqual(lik1,lik2,'absolute',tol,0)
assertVectorsAlmostEqual(VV1,VV2,'absolute',tol,0)
assertVectorsAlmostEqual(VVlag11,VVlag12,'absolute',tol,0)
assertVectorsAlmostEqual(VVlag21,VVlag22,'absolute',tol,0)



