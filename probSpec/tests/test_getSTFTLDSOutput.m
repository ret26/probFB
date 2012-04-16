function test_suite = test_getSTFTLDSOutput
  initTestSuite;

% Tests: 
%
% function [S,covS] = getSTFTLDSOutput(Xfin,Pfin);
%

function testComputeByHand

% Two time-step example which can be computed by hand

T = 3;
D = 2;

lamx = [.9;.7];
varx = rand(length(lamx),1);
om = [1/10,1/20]';
vary = rand;

y = randn(T,1);

[lik,Xfin,Pfin] = kalmanSlowSTFT(lamx,varx,om,vary,y);

[S1,covS1] = getSTFTLDSOutput(Xfin,Pfin);

S2 = zeros(D,T);
S2(1,:) = squeeze(Xfin(1,1,:)+i*Xfin(1,2,:));
S2(2,:) = squeeze(Xfin(1,3,:)+i*Xfin(1,4,:));

covS2 = [Pfin(1,1,:),Pfin(1,3,:),Pfin(1,2,:),Pfin(1,4,:); 
	 Pfin(3,1,:),Pfin(3,3,:),Pfin(3,2,:),Pfin(3,4,:);
	 Pfin(2,1,:),Pfin(2,3,:),Pfin(2,2,:),Pfin(2,4,:);
	 Pfin(4,1,:),Pfin(4,3,:),Pfin(4,2,:),Pfin(4,4,:)];

tol = 1e-5;
assertVectorsAlmostEqual(covS1,covS2,'absolute',tol,0)
assertVectorsAlmostEqual(S1,S2,'absolute',tol,0)
