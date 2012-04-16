function test_suite = test_probSTFT
  initTestSuite;

% Tests: 
%
% function [S,covS] = probSTFT(y,lamx,varx,om,vary)
%

function testCompareSTFT_And_FB

% Two time-step example which can be computed by hand

T = 3;
D = 2;

lamx = [.9;.7];
varx = rand(length(lamx),1);
om = [1/10,1/20]';
vary = rand;
y = randn(T,1);

[Z1,covZ1] = probFB(y,lamx,varx,om,vary);
[S,covS] = probSTFT(y,lamx,varx,om,vary);

Z2(1,:) = S(1,:).*exp(i*om(1)*[1:T]);
Z2(2,:) = S(2,:).*exp(i*om(2)*[1:T]);

covZ2 = zeros(size(covZ1));

for t=1:T
  c1 = cos(om(1)*t);
  s1 = -sin(om(1)*t);
  c2 = cos(om(2)*t);
  s2 = -sin(om(2)*t);
  R = [c1,0,s1,0;
       0,c2,0,s2;
       -s1,0,c1,0;
       0,-s2,0,c2];
       
  covZ2(:,:,t) = R*covS(:,:,t)*R';
end

tol = 1e-5;
assertVectorsAlmostEqual(Z1,Z2,'absolute',tol,0)
assertVectorsAlmostEqual(covZ1,covZ2,'absolute',tol,0)
