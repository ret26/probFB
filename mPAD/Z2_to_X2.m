function X2 = Z2_to_X2(Z2,Len2,Var2,tol)

% function X2 = Z2_to_X2(Z2,Len,tol)
%
% Derives the function values from the basis function coefficients
% The basis functions are Gaussians of width equal to Len2/sqrt(2) and
% height equal to Var2/Len2*sqrt(pi)
%
% INPUTS
% Z2 = basis-function coefficients [T,K]
% Len2 = GP length scale - relates to length-scales of basis functions [K,1]
% Var2 = GP power - relates to height of basis functions [K,1]
% tol = truncate basis function after tol*basis-function width [K,1]
%
% OUTPUTS
% X2 = function output [T,K]

%%%%%%

% form envelopes from basis functions

[T,K] = size(Z2);
  X2 = zeros(T,K);
  
  for k=1:K
    tau = ceil(Len2(k)/sqrt(2)*tol);
    bh = sqrt(Var2(k)*sqrt(2)/(Len2(k)*sqrt(pi))); % basis function height
                                           % to get correct variance
    bask = bh*exp(-1/(Len2(k)^2)*([-tau:1:tau]').^2);
    
    X2(:,k) =  conv(Z2(:,k),bask,'same');
  end
 